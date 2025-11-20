"""FASTQ alignment and variant calling for QuickVar."""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import os

from .install import ensure_environment, micromamba_run
from .reference import ensure_reference
from .settings import REFERENCE_REGISTRY
from .sra import download_and_convert_bioproject, fasterq_dump, prefetch_sra, query_study_runs

FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
DEFAULT_OUTPUT_NAME = "Results"


@dataclass
class Sample:
    name: str
    reads: List[Path]


def configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="[%(levelname)s] %(message)s")


def normalize_stem(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for suffix in (".fastq", ".fq"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name


def split_read_suffix(stem: str) -> tuple[str, int | None]:
    lowered = stem.lower()
    for marker in ("_r1", "-r1", ".r1", "_1", "-1", ".1"):
        if lowered.endswith(marker):
            return stem[: -len(marker)], 1
    for marker in ("_r2", "-r2", ".r2", "_2", "-2", ".2"):
        if lowered.endswith(marker):
            return stem[: -len(marker)], 2
    return stem, None


def discover_fastqs(input_path: Path) -> List[Path]:
    if input_path.is_file():
        return [input_path]
    if input_path.is_dir():
        return sorted(
            path
            for path in input_path.rglob("*")
            if any(str(path).lower().endswith(suffix) for suffix in FASTQ_SUFFIXES)
        )
    raise FileNotFoundError(f"Input path not found: {input_path}")


def group_samples(fastqs: Iterable[Path]) -> List[Sample]:
    samples: Dict[str, Dict[int, Path]] = {}
    singles: Dict[str, Path] = {}
    for fastq in fastqs:
        stem = normalize_stem(fastq)
        base, read = split_read_suffix(stem)
        if read is None:
            if base in singles:
                raise ValueError(f"Duplicate single-end sample detected: {base}")
            singles[base] = fastq
            continue
        entry = samples.setdefault(base, {})
        if read in entry:
            raise ValueError(f"Duplicate read {read} for sample {base}")
        entry[read] = fastq
    merged: List[Sample] = []
    merged_names: set[str] = set()
    for base, reads in samples.items():
        ordered = [reads[i] for i in sorted(reads)]
        if len(ordered) == 1:
            logging.warning("Sample %s missing mate; processing as single-end", base)
        merged.append(Sample(name=base, reads=ordered))
        merged_names.add(base)
    for base, fastq in singles.items():
        if base in samples:
            # Fallback: treat as separate sample if unpaired read exists.
            logging.warning("Sample %s has both paired and single reads; treating single as separate sample", base)
            suffix = 3
            new_name = f"{base}_single"
            while new_name in merged_names:
                suffix += 1
                new_name = f"{base}_single{suffix}"
            merged.append(Sample(name=new_name, reads=[fastq]))
            merged_names.add(new_name)
        else:
            merged.append(Sample(name=base, reads=[fastq]))
            merged_names.add(base)
    if not merged:
        raise ValueError("No FASTQ files detected")
    return merged


def align_sample(
    sample: Sample,
    reference: Dict[str, Path],
    output_dir: Path,
    threads: int,
    keep_intermediate: bool,
    ploidy: int,
    amplicon: bool,
    deduplicate: bool,
    annotate: bool,
    reference_key: str,
) -> None:
    logging.info("Processing sample %s", sample.name)
    sample_dir = output_dir / sample.name
    sample_dir.mkdir(parents=True, exist_ok=True)
    sam_path = sample_dir / f"{sample.name}.sam"
    bam_path = sample_dir / f"{sample.name}.sorted.bam"
    markdup_bam = sample_dir / f"{sample.name}.sorted.markdup.bam"
    bcf_path = sample_dir / f"{sample.name}.bcf"
    vcf_path = sample_dir / f"{sample.name}.vcf.gz"

    minimap_args = [
        "minimap2",
        "-t",
        str(threads),
        "-ax",
        "sr",
        str(reference["mmi"]),
        *map(str, sample.reads),
        "-o",
        str(sam_path),
    ]
    micromamba_run(minimap_args)

    samtools_sort_args = [
        "samtools",
        "sort",
        "-@",
        str(max(1, threads - 1)),
        "-o",
        str(bam_path),
        str(sam_path),
    ]
    micromamba_run(samtools_sort_args)

    final_bam = bam_path
    if deduplicate:
        micromamba_run(
            [
                "samtools",
                "markdup",
                "-@",
                str(max(1, threads - 1)),
                "-r",
                str(bam_path),
                str(markdup_bam),
            ]
        )
        markdup_bam.replace(bam_path)
        final_bam = bam_path


    micromamba_run(["samtools", "index", str(final_bam)])

    bcftools_mpileup = [
        "bcftools",
        "mpileup",
        "-Ob",
        "-f",
        str(reference["fasta"]),
        "-o",
        str(bcf_path),
        str(final_bam),
    ]
    micromamba_run(bcftools_mpileup)

    bcftools_call = [
        "bcftools",
        "call",
        "-mv",
        "-Oz",
        "--ploidy",
        str(ploidy),
        "-o",
        str(vcf_path),
        str(bcf_path),
    ]
    micromamba_run(bcftools_call)

    micromamba_run(["bcftools", "index", str(vcf_path)])

    # Optional VCF annotation (with consequences)
    if annotate:
        try:
            annotation_bgz = ensure_annotations(reference_key, reference)
            if annotation_bgz:
                annotate_vcf(vcf_path, annotation_bgz, reference)
            else:
                logging.warning("VCF annotation skipped: annotation file not available")
        except Exception as e:
            logging.warning(f"VCF annotation failed (continuing without annotation): {e}")

    if not keep_intermediate:
        if sam_path.exists():
            sam_path.unlink()
        if bcf_path.exists():
            bcf_path.unlink()
        if deduplicate and markdup_bam.exists():
            markdup_bam.unlink()
            markdup_index = markdup_bam.with_suffix(".bam.bai")
            if markdup_index.exists():
                markdup_index.unlink()

    if amplicon:
        generate_amplicon_report(
            sample=sample,
            reference=reference,
            bam_path=final_bam,
            sample_dir=sample_dir,
        )


def parse_pileup_bases(bases: str, ref: str) -> Dict[str, int]:
    counts = defaultdict(int)
    i = 0
    ref_upper = ref.upper()
    while i < len(bases):
        base = bases[i]
        if base == "^":
            i += 2
            continue
        if base == "$":
            i += 1
            continue
        if base in "+-":
            i += 1
            length_digits = []
            while i < len(bases) and bases[i].isdigit():
                length_digits.append(bases[i])
                i += 1
            length = int("".join(length_digits)) if length_digits else 0
            indel_seq = bases[i : i + length].upper()
            i += length
            if indel_seq:
                key = f"{base}{indel_seq}"
                counts[key] += 1
            continue
        if base == "*":
            i += 1
            continue
        if base in {".", ","}:
            counts[ref_upper] += 1
        else:
            counts[base.upper()] += 1
        i += 1
    return counts


def generate_amplicon_report(
    sample: Sample,
    reference: Dict[str, Path],
    bam_path: Path,
    sample_dir: Path,
) -> None:
    logging.info("Generating amplicon summary for sample %s", sample.name)
    igv_depths = load_igv_depths(bam_path)
    result = micromamba_run(
        [
            "samtools",
            "mpileup",
            "-aa",
            "-d",
            "0",
            "-f",
            str(reference["fasta"]),
            str(bam_path),
        ],
        capture_output=True,
        text=True,
    )
    summary_path = sample_dir / f"{sample.name}_amplicon.tsv"
    indel_records: list[dict[str, object]] = []
    position_data: Dict[tuple[str, int], tuple[str, Dict[str, int], int]] = {}
    with open(summary_path, "w", encoding="utf-8") as handle:
        handle.write(
            (
                "chrom\tpos\tref_base\talt_base\talt_count\tdepth\tfrequency\tmutation\t"
                "igv_depth\testimated_coverage\testimated_frequency\n"
            )
        )
        for line in result.stdout.splitlines():
            fields = line.strip().split("\t")
            if len(fields) < 6:
                continue
            chrom, pos_str, ref_base, depth_str, bases, _quals = fields[:6]
            try:
                depth = int(depth_str)
                pos_int = int(pos_str)
            except ValueError:
                continue
            if depth <= 0:
                continue
            counts = parse_pileup_bases(bases, ref_base)
            total_depth = sum(counts.values())
            if total_depth == 0:
                continue
            ref_upper = ref_base.upper()
            igv_depth = igv_depths.get((chrom, pos_int), 0)
            estimated_cov = estimate_neighbor_depth(igv_depths, chrom, pos_int)
            position_data[(chrom, pos_int)] = (ref_upper, counts.copy(), total_depth)
            for base, count in counts.items():
                if base == ref_upper or count == 0:
                    continue
                frequency = count / total_depth if total_depth else 0.0
                estimated_freq = count / estimated_cov if estimated_cov > 0 else 0.0
                if base.startswith("+"):
                    mutation = f"{ref_upper}>+{base[1:]}"
                elif base.startswith("-"):
                    mutation = f"{ref_upper}>-{base[1:]}"
                else:
                    mutation = f"{ref_upper}>{base}"
                handle.write(
                    f"{chrom}\t{pos_int}\t{ref_upper}\t{base}\t{count}\t{total_depth}\t"
                    f"{frequency:.4f}\t{mutation}\t{igv_depth}\t{estimated_cov:.2f}\t"
                    f"{estimated_freq:.4f}\n"
                )
                if base.startswith("+") or base.startswith("-"):
                    indel_records.append(
                        {
                            "chrom": chrom,
                            "pos": pos_int,
                            "ref": ref_upper,
                            "alt": base,
                            "alt_count": count,
                            "total_depth": total_depth,
                            "frequency": frequency,
                            "estimated_frequency": estimated_freq,
                        }
                    )

    if indel_records:
        indel_path = sample_dir / f"{sample.name}_amplicon_indels.tsv"
        with open(indel_path, "w", encoding="utf-8") as indel_file:
            indel_file.write(
                "chrom\twindow_start\twindow_end\tpos\tref_base\talt_base\talt_count\t"
                "total_depth\twt_count_in_10bp\tfrequency\testimated_frequency\n"
            )
            for record in indel_records:
                pos_int = int(record["pos"])
                window_start = pos_int - 5
                window_end = pos_int + 5
                total_depth = int(record["total_depth"])
                alt_count = int(record["alt_count"])
                wt_count = estimate_wildtype_window(position_data, record["chrom"], pos_int, flank=5)
                indel_file.write(
                    f"{record['chrom']}\t{window_start}\t{window_end}\t{pos_int}\t"
                    f"{record['ref']}\t{record['alt']}\t{alt_count}\t{total_depth}\t"
                    f"{wt_count}\t{record['frequency']:.4f}\t"
                    f"{record['estimated_frequency']:.4f}\n"
                )


def estimate_wildtype_window(
    position_data: Dict[tuple[str, int], tuple[str, Dict[str, int], int]],
    chrom: str,
    center_pos: int,
    flank: int = 5,
) -> int:
    ref_counts: list[int] = []
    for offset in range(-flank, flank + 1):
        key = (chrom, center_pos + offset)
        data = position_data.get(key)
        if not data:
            continue
        ref_base, counts, _ = data
        ref_counts.append(counts.get(ref_base, 0))
    if not ref_counts:
        return 0
    return min(ref_counts)


def load_igv_depths(bam_path: Path) -> Dict[tuple[str, int], int]:
    result = micromamba_run(
        [
            "samtools",
            "depth",
            "-aa",
            "-d",
            "0",
            "-Q",
            "0",
            "-q",
            "0",
            str(bam_path),
        ],
        capture_output=True,
        text=True,
    )
    depths: Dict[tuple[str, int], int] = {}
    for line in result.stdout.splitlines():
        parts = line.strip().split("\t")
        if len(parts) != 3:
            continue
        chrom, pos_str, depth_str = parts
        try:
            depths[(chrom, int(pos_str))] = int(depth_str)
        except ValueError:
            continue
    return depths


def estimate_neighbor_depth(
    depths: Dict[tuple[str, int], int],
    chrom: str,
    center_pos: int,
    flank: int = 5,
) -> float:
    values: list[int] = []
    for offset in range(-flank, flank + 1):
        if offset == 0:
            continue
        key = (chrom, center_pos + offset)
        if key in depths:
            depth_value = depths[key]
            if depth_value > 0:
                values.append(depth_value)
    if not values:
        return float(depths.get((chrom, center_pos), 0))
    return sum(values) / len(values)


def parse_gff_attributes(attr_string: str) -> Dict[str, str]:
    """Parse GFF attribute string into a dictionary."""
    attrs: Dict[str, str] = {}
    for pair in attr_string.split(";"):
        if "=" in pair:
            key, value = pair.split("=", 1)
            attrs[key.strip()] = value.strip()
    return attrs


def parse_gff_annotations(gff_path: Path) -> Dict[Tuple[str, int, int], Dict[str, str]]:
    """
    Parse GFF file and extract gene annotations.
    Returns a dictionary mapping (chrom, start, end) -> {gene_id, gene_name, feature_type}.
    """
    annotations: Dict[Tuple[str, int, int], Dict[str, str]] = {}
    
    if not gff_path.exists():
        return annotations
    
    open_func = gzip.open if gff_path.suffix == ".gz" or str(gff_path).endswith(".gz") else open
    mode = "rt" if open_func == open else "rt"
    
    try:
        with open_func(gff_path, mode) as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                
                fields = line.split("\t")
                if len(fields) < 9:
                    continue
                
                chrom = fields[0]
                feature_type = fields[2]
                try:
                    start = int(fields[3])
                    end = int(fields[4])
                except ValueError:
                    continue
                
                # Focus on gene and CDS features
                if feature_type not in ("gene", "CDS", "mRNA"):
                    continue
                
                attrs = parse_gff_attributes(fields[8])
                
                # Extract gene ID and name
                gene_id = attrs.get("ID", attrs.get("gene_id", attrs.get("locus_tag", "")))
                gene_name = attrs.get("Name", attrs.get("gene", attrs.get("gene_name", "")))
                product = attrs.get("product", attrs.get("Note", ""))
                
                # Store annotation for this region
                key = (chrom, start, end)
                if key not in annotations or feature_type == "gene":
                    annotations[key] = {
                        "gene_id": gene_id,
                        "gene_name": gene_name,
                        "feature_type": feature_type,
                        "product": product,
                    }
    except Exception as e:
        logging.warning(f"Failed to parse GFF file {gff_path}: {e}")
        return {}
    
    return annotations


def parse_gff_cds_features(gff_path: Path) -> Dict[str, List[Dict[str, any]]]:
    """
    Parse GFF file and extract CDS features grouped by transcript/gene.
    Returns a dictionary mapping gene_id -> list of CDS features with coordinates, strand, phase.
    """
    cds_features: Dict[str, List[Dict[str, any]]] = {}
    
    if not gff_path.exists():
        return cds_features
    
    open_func = gzip.open if gff_path.suffix == ".gz" or str(gff_path).endswith(".gz") else open
    mode = "rt" if open_func == open else "rt"
    
    try:
        with open_func(gff_path, mode) as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                
                fields = line.split("\t")
                if len(fields) < 9:
                    continue
                
                feature_type = fields[2]
                if feature_type != "CDS":
                    continue
                
                try:
                    chrom = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    strand = fields[6]  # + or -
                    phase = fields[7] if len(fields) > 7 else "0"  # 0, 1, or 2
                except (ValueError, IndexError):
                    continue
                
                attrs = parse_gff_attributes(fields[8])
                parent = attrs.get("Parent", "")
                gene_id = attrs.get("ID", parent)  # Use Parent or ID as gene identifier
                
                # Extract gene ID from parent if available
                if "Parent" in attrs:
                    # Try to get gene ID from parent transcript
                    parent_id = attrs["Parent"]
                    # For now, use parent ID as the key (we can improve this later)
                    gene_id = parent_id.split("-")[0] if "-" in parent_id else parent_id
                
                # Also check for gene-level ID
                gene_id_from_gene = attrs.get("gene_id", "")
                if gene_id_from_gene:
                    gene_id = gene_id_from_gene
                
                # Use transcript ID as key (most specific)
                transcript_id = parent if parent else gene_id
                
                if transcript_id:
                    if transcript_id not in cds_features:
                        cds_features[transcript_id] = []
                    
                    cds_features[transcript_id].append({
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "phase": int(phase) if phase.isdigit() else 0,
                        "gene_id": gene_id,
                    })
    except Exception as e:
        logging.warning(f"Failed to parse CDS features from GFF file {gff_path}: {e}")
        return {}
    
    # Sort CDS features by start position for each transcript
    for transcript_id in cds_features:
        cds_features[transcript_id].sort(key=lambda x: x["start"])
    
    return cds_features


# Standard genetic code
GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Amino acid three-letter to one-letter mapping
AA_3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Stop": "*",
}


_REFERENCE_SEQUENCE_CACHE: Dict[Path, Dict[str, str]] = {}


def _load_reference_sequences(reference_fasta: Path) -> Dict[str, str]:
    """
    Load all reference chromosome sequences into memory once.
    Returns a dict: chrom_id -> sequence (uppercased).
    """
    if reference_fasta in _REFERENCE_SEQUENCE_CACHE:
        return _REFERENCE_SEQUENCE_CACHE[reference_fasta]

    sequences: Dict[str, str] = {}
    current_id: Optional[str] = None
    chunks: List[str] = []

    try:
        with open(reference_fasta, "r", encoding="utf-8") as handle:
            for line in handle:
                if not line:
                    continue
                if line.startswith(">"):
                    # Flush previous record
                    if current_id is not None:
                        sequences[current_id] = "".join(chunks).upper()
                    header = line[1:].strip()
                    current_id = header.split()[0]
                    chunks = []
                else:
                    chunks.append(line.strip())
            # Flush last record
            if current_id is not None:
                sequences[current_id] = "".join(chunks).upper()
    except OSError as e:
        logging.warning("Failed to load reference FASTA %s: %s", reference_fasta, e)
        return {}

    _REFERENCE_SEQUENCE_CACHE[reference_fasta] = sequences
    return sequences


def get_reference_sequence(reference: Dict[str, Path], chrom: str, start: int, end: int) -> str:
    """
    Extract reference sequence for a given region using an in-memory cache.
    Coordinates are 1-based inclusive.
    """
    fasta_path = reference["fasta"]
    sequences = _load_reference_sequences(fasta_path)
    if not sequences:
        return ""

    # Try exact match, then relaxed matching on IDs
    seq = sequences.get(chrom)
    if seq is None:
        for key in sequences:
            if key == chrom or chrom in key or key in chrom:
                seq = sequences[key]
                break
    if seq is None:
        return ""

    # Clamp coordinates to sequence bounds
    start_idx = max(start - 1, 0)
    end_idx = min(end, len(seq))
    if start_idx >= end_idx:
        return ""
    return seq[start_idx:end_idx]


def calculate_codon_position(variant_pos: int, cds_start: int, cds_end: int, strand: str, phase: int = 0) -> Optional[int]:
    """Calculate which codon (1-based) a variant position falls in within a CDS."""
    if not (cds_start <= variant_pos <= cds_end):
        return None
    
    if strand == "+":
        # Forward strand: position is relative to CDS start
        cds_pos = variant_pos - cds_start + 1  # 1-based position in CDS
    else:
        # Reverse strand: position is relative to CDS end (reverse complement)
        cds_pos = cds_end - variant_pos + 1  # 1-based position in CDS
    
    # Account for phase offset
    cds_pos_with_phase = cds_pos - phase
    
    # Calculate codon number (1-based)
    codon_pos = (cds_pos_with_phase - 1) // 3 + 1
    
    return codon_pos


def translate_codon(codon: str) -> str:
    """Translate a DNA codon to amino acid (single letter)."""
    codon_upper = codon.upper()
    if len(codon_upper) != 3:
        return "?"
    return GENETIC_CODE.get(codon_upper, "?")


def calculate_consequence(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    reference: Dict[str, Path],
    cds_index: Dict[str, List[Dict[str, object]]],
    gene_index: Dict[str, List[Dict[str, object]]],
) -> Dict[str, str]:
    """
    Calculate variant consequence for a given variant.
    Returns a dictionary with consequence information.
    """
    consequence_info = {
        "consequence": "intergenic",
        "amino_acid_change": "",
        "codon_position": "",
        "protein_position": "",
    }
    
    # Check if variant overlaps with any CDS using indexed intervals
    overlapping_cds: Optional[Dict[str, object]] = None

    cds_list = cds_index.get(chrom, [])
    if cds_list:
        # Binary search by start coordinate
        starts = [c["start"] for c in cds_list]  # type: ignore[index]
        lo = 0
        hi = len(starts)
        while lo < hi:
            mid = (lo + hi) // 2
            if starts[mid] <= pos:
                lo = mid + 1
            else:
                hi = mid
        # Scan a small window backwards to find overlapping interval(s)
        i = max(0, lo - 5)
        while i < len(cds_list) and cds_list[i]["start"] <= pos:  # type: ignore[index]
            cds = cds_list[i]
            if cds["start"] <= pos <= cds["end"]:  # type: ignore[index]
                overlapping_cds = cds
                break
            i += 1

    if not overlapping_cds:
        # Variant is intergenic or in intron/UTR
        return consequence_info
    
    # Variant is in CDS - calculate consequence
    cds_start = overlapping_cds["start"]
    cds_end = overlapping_cds["end"]
    strand = overlapping_cds["strand"]
    phase = overlapping_cds.get("phase", 0)
    
    # Calculate codon position
    codon_pos = calculate_codon_position(pos, cds_start, cds_end, strand, phase)
    
    if codon_pos is None:
        return consequence_info
    
    # Determine consequence type
    ref_len = len(ref)
    alt_len = len(alt)
    
    if ref_len == alt_len == 1:
        # SNP
        # Get reference codon sequence
        # Calculate which position in CDS the variant is at
        if strand == "+":
            # Forward strand: position in CDS (1-based) = variant_pos - cds_start + 1
            cds_pos = pos - cds_start + 1  # 1-based position in CDS
        else:
            # Reverse strand: position in CDS (1-based) = cds_end - variant_pos + 1
            cds_pos = cds_end - pos + 1  # 1-based position in CDS
        
        # Account for phase offset
        cds_pos_with_phase = cds_pos - phase
        
        # Calculate codon boundaries (1-based)
        codon_start_cds = ((codon_pos - 1) * 3) + 1  # Start of codon in CDS coordinates
        codon_end_cds = codon_start_cds + 2  # End of codon in CDS coordinates
        
        # Convert CDS coordinates to genomic coordinates
        if strand == "+":
            codon_genomic_start = cds_start + codon_start_cds - 1  # 1-based to 0-based
            codon_genomic_end = cds_start + codon_end_cds - 1
        else:
            # Reverse strand: reverse the coordinates
            codon_genomic_end = cds_end - codon_start_cds + 1
            codon_genomic_start = cds_end - codon_end_cds + 1
        
        # Get reference codon sequence
        ref_codon_raw = get_reference_sequence(reference, chrom, codon_genomic_start, codon_genomic_end)
        if not ref_codon_raw or len(ref_codon_raw) != 3:
            consequence_info["consequence"] = "coding_sequence_variant"
            consequence_info["codon_position"] = str(codon_pos)
            return consequence_info
        
        # For reverse strand, reverse complement the codon (get reference sequence is always forward)
        complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
        if strand == "-":
            # Reverse complement: reverse the sequence and complement each base
            ref_codon = "".join(complement.get(b, "N") for b in reversed(ref_codon_raw))
            # For alt, we need to complement it too
            alt_comp = complement.get(alt, alt)
        else:
            ref_codon = ref_codon_raw
            alt_comp = alt
        
        # Calculate position within codon (0-based)
        pos_in_cds = cds_pos_with_phase  # Position in CDS (accounting for phase)
        pos_in_codon = (pos_in_cds - 1) % 3  # 0-based position within codon (0, 1, or 2)
        
        # Create alt codon
        alt_codon = ref_codon[:pos_in_codon] + alt_comp + ref_codon[pos_in_codon + 1:]
        
        ref_aa = translate_codon(ref_codon)
        alt_aa = translate_codon(alt_codon)
        
        # Convert single-letter to three-letter amino acid codes for output
        aa_1_to_3 = {
            "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
            "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
            "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
            "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
            "*": "Ter", "?": "?",
        }
        
        ref_aa_3 = aa_1_to_3.get(ref_aa, ref_aa)
        alt_aa_3 = aa_1_to_3.get(alt_aa, alt_aa)
        protein_pos = codon_pos
        
        if ref_aa == alt_aa:
            consequence_info["consequence"] = "synonymous_variant"
        elif alt_aa == "*":
            consequence_info["consequence"] = "stop_gained"
            consequence_info["amino_acid_change"] = f"p.{ref_aa_3}{protein_pos}Ter"
        elif ref_aa == "*":
            consequence_info["consequence"] = "stop_lost"
            consequence_info["amino_acid_change"] = f"p.Ter{protein_pos}{alt_aa_3}"
        else:
            consequence_info["consequence"] = "missense_variant"
            consequence_info["amino_acid_change"] = f"p.{ref_aa_3}{protein_pos}{alt_aa_3}"
            consequence_info["protein_position"] = str(protein_pos)
        
        consequence_info["codon_position"] = str(codon_pos)
    
    elif ref_len > alt_len:
        # Deletion
        del_len = ref_len - alt_len
        if del_len % 3 == 0:
            consequence_info["consequence"] = "inframe_deletion"
        else:
            consequence_info["consequence"] = "frameshift_variant"
            consequence_info["amino_acid_change"] = f"p.{codon_pos}fs"
    
    elif alt_len > ref_len:
        # Insertion
        ins_len = alt_len - ref_len
        if ins_len % 3 == 0:
            consequence_info["consequence"] = "inframe_insertion"
        else:
            consequence_info["consequence"] = "frameshift_variant"
            consequence_info["amino_acid_change"] = f"p.{codon_pos}fs"
    
    else:
        # Complex variant
        consequence_info["consequence"] = "coding_sequence_variant"
    
    if not consequence_info["amino_acid_change"] and codon_pos:
        consequence_info["codon_position"] = str(codon_pos)
    
    return consequence_info


def build_annotation_tsv(
    gff_path: Path,
    output_tsv: Path,
    reference: Dict[str, Path],
) -> bool:
    """
    Build a TSV annotation file from GFF for bcftools annotate.
    Format: CHROM\tPOS\tGENE_ID\tGENE_NAME\tFEATURE_TYPE
    Returns True if successful, False otherwise.
    """
    try:
        annotations = parse_gff_annotations(gff_path)
        if not annotations:
            logging.warning("No annotations found in GFF file")
            return False
        
        # Read reference FASTA to get chromosome names
        chrom_names: set[str] = set()
        fai_path = reference["fasta"].with_suffix(".fai")
        if fai_path.exists():
            with open(fai_path) as f:
                for line in f:
                    chrom_names.add(line.split("\t")[0])
        
        # Build position-based annotation lookup
        position_annotations: Dict[Tuple[str, int], Dict[str, str]] = {}
        for (chrom, start, end), info in annotations.items():
            if chrom_names and chrom not in chrom_names:
                # Try to match chromosome name variations
                matched = False
                for ref_chrom in chrom_names:
                    if chrom in ref_chrom or ref_chrom in chrom:
                        chrom = ref_chrom
                        matched = True
                        break
                if not matched and chrom_names:
                    continue
            
            # Store annotation for each position in the range
            for pos in range(start, end + 1):
                key = (chrom, pos)
                # Prefer gene-level annotations over CDS/mRNA
                if key not in position_annotations or info["feature_type"] == "gene":
                    position_annotations[key] = info
        
        if not position_annotations:
            logging.warning("No position annotations could be matched to reference chromosomes")
            return False
        
        # Write TSV file
        with open(output_tsv, "w") as f:
            f.write("CHROM\tPOS\tGENE_ID\tGENE_NAME\tFEATURE_TYPE\tPRODUCT\n")
            for (chrom, pos), info in sorted(position_annotations.items()):
                f.write(
                    f"{chrom}\t{pos}\t{info['gene_id']}\t{info['gene_name']}\t"
                    f"{info['feature_type']}\t{info['product']}\n"
                )
        
        logging.info(f"Built annotation TSV with {len(position_annotations)} positions")
        return True
    except Exception as e:
        logging.warning(f"Failed to build annotation TSV: {e}")
        return False


def index_annotation_tsv(tsv_path: Path) -> bool:
    """
    Index annotation TSV with tabix for bcftools annotate.
    Returns True if successful, False otherwise.
    """
    try:
        # Sort the TSV file (required for tabix) - skip header
        sorted_tsv = tsv_path.with_suffix(".sorted.tsv")
        with open(tsv_path) as infile, open(sorted_tsv, "w") as outfile:
            header = infile.readline()
            outfile.write(header)
            # Sort remaining lines
            lines = sorted(infile, key=lambda x: (x.split("\t")[0], int(x.split("\t")[1])))
            outfile.writelines(lines)
        
        # Replace original with sorted
        sorted_tsv.replace(tsv_path)
        
        # Compress with bgzip (part of htslib, should be available via bcftools)
        # bgzip compresses in place and renames to .gz, so we'll rename to .bgz after
        bgz_path = tsv_path.with_suffix(".tsv.bgz")
        result = micromamba_run(
            ["bgzip", str(tsv_path)],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            logging.warning(f"bgzip failed: {result.stderr or 'unknown error'}")
            return False
        
        # bgzip creates .gz file, rename to .bgz for clarity
        gz_path = tsv_path.with_suffix(".tsv.gz")
        if gz_path.exists():
            gz_path.rename(bgz_path)
        elif not bgz_path.exists():
            logging.warning("bgzip did not create expected output file")
            return False
        
        # Index with tabix (part of htslib)
        # For TSV format: -s 1 (sequence/chrom), -b 2 (start/pos), -e 2 (end/pos)
        # Skip header with -S 1
        result = micromamba_run(
            ["tabix", "-S", "1", "-s", "1", "-b", "2", "-e", "2", str(bgz_path)],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            error_msg = result.stderr or result.stdout or "unknown error"
            logging.warning(f"tabix failed: {error_msg}")
            # Don't delete the bgz file - it might be useful for debugging
            return False
        
        # Remove uncompressed TSV (already compressed)
        if tsv_path.exists():
            tsv_path.unlink()
        
        logging.info(f"Indexed annotation file: {bgz_path}")
        return True
    except Exception as e:
        import traceback
        logging.warning(f"Failed to index annotation TSV: {e}")
        logging.debug(f"Traceback: {traceback.format_exc()}")
        return False


def annotate_vcf(
    vcf_path: Path,
    annotation_bgz: Path,
    reference: Dict[str, Path],
) -> bool:
    """
    Annotate VCF file with gene information and consequences using bcftools annotate.
    Adds consequences (missense, frameshift, etc.) calculated from GFF and reference.
    Returns True if successful, False otherwise.
    """
    import tempfile
    
    try:
        if not annotation_bgz.exists():
            logging.warning(f"Annotation file not found: {annotation_bgz}")
            return False
        
        annotated_vcf = vcf_path.with_suffix(".annotated.vcf.gz")
        
        # Create temporary header file for INFO field definitions
        header_lines = [
            "##INFO=<ID=GENE_ID,Number=1,Type=String,Description=\"Gene ID\">",
            "##INFO=<ID=GENE_NAME,Number=1,Type=String,Description=\"Gene name\">",
            "##INFO=<ID=FEATURE_TYPE,Number=1,Type=String,Description=\"Feature type (gene/CDS/mRNA)\">",
            "##INFO=<ID=PRODUCT,Number=1,Type=String,Description=\"Gene product/description\">",
            "##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description=\"Variant consequence (missense_variant, frameshift_variant, synonymous_variant, etc.)\">",
            "##INFO=<ID=AMINO_ACID_CHANGE,Number=1,Type=String,Description=\"Amino acid change (e.g., p.Arg123Lys)\">",
            "##INFO=<ID=CODON_POSITION,Number=1,Type=Integer,Description=\"Codon position in CDS\">",
            "##INFO=<ID=PROTEIN_POSITION,Number=1,Type=Integer,Description=\"Protein position\">",
        ]
        
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as header_file:
            header_file.write("\n".join(header_lines) + "\n")
            header_path = Path(header_file.name)
        
        try:
            # First: annotate with gene information using bcftools
            result = micromamba_run(
                [
                    "bcftools",
                    "annotate",
                    "-a",
                    str(annotation_bgz),
                    "-c",
                    "CHROM,POS,GENE_ID,GENE_NAME,FEATURE_TYPE,PRODUCT",
                    "-h",
                    str(header_path),
                    "-o",
                    str(annotated_vcf),
                    str(vcf_path),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            
            if result.returncode != 0:
                error_msg = result.stderr or result.stdout or "unknown error"
                logging.warning(f"bcftools annotate failed: {error_msg}")
                return False
            
            # Second: add consequences by parsing VCF and calculating from GFF
            # Try to find GFF file from reference registry
            gff_path = None
            for ref_key in REFERENCE_REGISTRY:
                ref_metadata = REFERENCE_REGISTRY[ref_key]
                if "local_gff" in ref_metadata:
                    gff_candidate = ref_metadata["local_gff"]
                    if gff_candidate:
                        gff_candidate_path = Path(gff_candidate)
                        if gff_candidate_path.exists():
                            gff_path = gff_candidate_path
                            break
            
            if gff_path and gff_path.exists():
                logging.info("Calculating variant consequences from GFF...")
                add_consequences_to_vcf(annotated_vcf, reference, gff_path)
            else:
                logging.debug("GFF file not found - skipping consequence calculation")
            
            # Replace original VCF with annotated version
            annotated_vcf.replace(vcf_path)
            
            # Re-index
            micromamba_run(["bcftools", "index", str(vcf_path)], check=False)
            
            header_path.unlink()
            return True
        except Exception as e:
            logging.warning(f"Failed to annotate VCF: {e}")
            if header_path.exists():
                header_path.unlink()
            return False
    except Exception as e:
        logging.warning(f"VCF annotation failed: {e}")
        return False


def add_consequences_to_vcf(
    vcf_path: Path,
    reference: Dict[str, Path],
    gff_path: Path,
) -> None:
    """
    Add consequence annotations to VCF file by parsing variants and calculating consequences.
    Updates INFO fields: CONSEQUENCE, AMINO_ACID_CHANGE, CODON_POSITION, PROTEIN_POSITION
    """
    import tempfile
    import gzip
    
    try:
        # Parse CDS features and gene annotations once
        cds_features = parse_gff_cds_features(gff_path)
        gene_annotations = parse_gff_annotations(gff_path)
        
        if not cds_features:
            logging.debug("No CDS features found for consequence calculation")
            return

        # Build per-chromosome interval indexes for fast overlap queries
        cds_index: Dict[str, List[Dict[str, object]]] = {}
        for transcript_id, cds_list in cds_features.items():
            for cds in cds_list:
                chrom = cds["chrom"]  # type: ignore[index]
                cds_entry = {
                    "chrom": chrom,
                    "start": cds["start"],  # type: ignore[index]
                    "end": cds["end"],      # type: ignore[index]
                    "strand": cds["strand"],  # type: ignore[index]
                    "phase": cds.get("phase", 0),
                    "gene_id": cds.get("gene_id", transcript_id),
                }
                cds_index.setdefault(chrom, []).append(cds_entry)
        # Sort intervals by start for each chromosome
        for chrom in cds_index:
            cds_index[chrom].sort(key=lambda x: x["start"])  # type: ignore[index]

        gene_index: Dict[str, List[Dict[str, object]]] = {}
        for (gff_chrom, start, end), info in gene_annotations.items():
            entry = {
                "start": start,
                "end": end,
                "gene_id": info.get("gene_id", ""),
                "gene_name": info.get("gene_name", ""),
                "feature_type": info.get("feature_type", ""),
                "product": info.get("product", ""),
            }
            gene_index.setdefault(gff_chrom, []).append(entry)
        for chrom in gene_index:
            gene_index[chrom].sort(key=lambda x: x["start"])  # type: ignore[index]
        
        # Read VCF file
        open_func = gzip.open if str(vcf_path).endswith(".gz") else open
        mode = "rt" if open_func == open else "rt"
        
        # Create temporary output VCF
        temp_vcf = vcf_path.with_suffix(".temp.vcf.gz")
        
        with open_func(vcf_path, mode) as vcf_in:
            with gzip.open(temp_vcf, "wt") as vcf_out:
                for line in vcf_in:
                    if line.startswith("##"):
                        # Header line - pass through
                        vcf_out.write(line)
                    elif line.startswith("#CHROM"):
                        # Column header - pass through
                        vcf_out.write(line)
                    else:
                        # Variant line - calculate and add consequences
                        fields = line.strip().split("\t")
                        if len(fields) < 5:
                            vcf_out.write(line)
                            continue
                        
                        chrom = fields[0]
                        pos_str = fields[1]
                        ref = fields[3]
                        alt = fields[4]
                        info = fields[7] if len(fields) > 7 else "."
                        
                        try:
                            pos = int(pos_str)
                        except ValueError:
                            vcf_out.write(line)
                            continue
                        
                        # Calculate consequence using indexed annotations
                        consequence_info = calculate_consequence(
                            chrom, pos, ref, alt, reference, cds_index, gene_index
                        )
                        
                        # Add consequence fields to INFO
                        if consequence_info["consequence"] and consequence_info["consequence"] != "intergenic":
                            info_fields = info.split(";")
                            
                            # Remove existing consequence fields if present
                            info_fields = [f for f in info_fields if not f.startswith(("CONSEQUENCE=", "AMINO_ACID_CHANGE=", "CODON_POSITION=", "PROTEIN_POSITION="))]
                            
                            # Add new consequence fields
                            info_fields.append(f"CONSEQUENCE={consequence_info['consequence']}")
                            if consequence_info["amino_acid_change"]:
                                info_fields.append(f"AMINO_ACID_CHANGE={consequence_info['amino_acid_change']}")
                            if consequence_info["codon_position"]:
                                info_fields.append(f"CODON_POSITION={consequence_info['codon_position']}")
                            if consequence_info["protein_position"]:
                                info_fields.append(f"PROTEIN_POSITION={consequence_info['protein_position']}")
                            
                            fields[7] = ";".join(info_fields)
                        
                        vcf_out.write("\t".join(fields) + "\n")
        
        # Replace original with annotated version
        temp_vcf.replace(vcf_path)
        
    except Exception as e:
        logging.warning(f"Failed to add consequences to VCF: {e}")


def annotate_vcf_basic(
    vcf_path: Path,
    annotation_bgz: Path,
    reference: Dict[str, Path],
) -> bool:
    """
    Basic VCF annotation with gene information only (no consequences).
    Returns True if successful, False otherwise.
    """
    import tempfile
    
    try:
        if not annotation_bgz.exists():
            logging.warning(f"Annotation file not found: {annotation_bgz}")
            return False
        
        annotated_vcf = vcf_path.with_suffix(".annotated.vcf.gz")
        
        # Create temporary header file for INFO field definitions
        header_lines = [
            "##INFO=<ID=GENE_ID,Number=1,Type=String,Description=\"Gene ID\">",
            "##INFO=<ID=GENE_NAME,Number=1,Type=String,Description=\"Gene name\">",
            "##INFO=<ID=FEATURE_TYPE,Number=1,Type=String,Description=\"Feature type (gene/CDS/mRNA)\">",
            "##INFO=<ID=PRODUCT,Number=1,Type=String,Description=\"Gene product/description\">",
        ]
        
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as header_file:
            header_file.write("\n".join(header_lines) + "\n")
            header_path = Path(header_file.name)
        
        try:
            result = micromamba_run(
                [
                    "bcftools",
                    "annotate",
                    "-a",
                    str(annotation_bgz),
                    "-c",
                    "CHROM,POS,GENE_ID,GENE_NAME,FEATURE_TYPE,PRODUCT",
                    "-h",
                    str(header_path),
                    "-o",
                    str(annotated_vcf),
                    str(vcf_path),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            
            if result.returncode != 0:
                error_msg = result.stderr or result.stdout or "unknown error"
                logging.warning(f"bcftools annotate failed: {error_msg}")
                return False
            
            # Replace original VCF with annotated version
            annotated_vcf.replace(vcf_path)
            
            # Re-index
            micromamba_run(["bcftools", "index", str(vcf_path)], check=False)
            
            logging.info(f"Successfully annotated VCF: {vcf_path}")
            return True
        finally:
            # Clean up temporary header file
            if header_path.exists():
                header_path.unlink()
    except Exception as e:
        logging.warning(f"Failed to annotate VCF: {e}")
        return False


def ensure_annotations(
    reference_key: str,
    reference: Dict[str, Path],
) -> Optional[Path]:
    """
    Ensure annotation files are available for the given reference.
    Returns path to annotation BGZ file if successful, None otherwise.
    """
    if reference_key not in REFERENCE_REGISTRY:
        return None
    
    metadata = REFERENCE_REGISTRY[reference_key]
    gff_path = metadata.get("local_gff")
    
    if not gff_path or not gff_path.exists():
        logging.debug(f"No GFF file available for reference {reference_key}")
        return None
    
    # Check cache for existing annotation file
    from .settings import REFERENCE_DIR
    
    annotation_bgz = REFERENCE_DIR / f"{reference_key}_annotations.tsv.bgz"
    annotation_tsv = REFERENCE_DIR / f"{reference_key}_annotations.tsv"
    
    tbi_path = annotation_bgz.with_name(annotation_bgz.name + ".tbi")
    if annotation_bgz.exists() and tbi_path.exists():
        return annotation_bgz
    
    # Build annotation TSV
    if not build_annotation_tsv(gff_path, annotation_tsv, reference):
        return None
    
    # Index annotation TSV
    if not index_annotation_tsv(annotation_tsv):
        return None
    
    return annotation_bgz


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Align FASTQ reads and call variants")
    parser.add_argument(
        "--input",
        help="Path to FASTQ file or directory containing FASTQs (required if --bioproject not used)",
    )
    parser.add_argument(
        "--bioproject",
        help="NCBI BioProject ID (e.g., PRJNA123456) to download and process SRA files",
    )
    parser.add_argument("--output", default=DEFAULT_OUTPUT_NAME, help="Directory to write results (default: Results)")
    parser.add_argument("--threads", type=int, default=0, help="Number of CPU threads (default: auto)")
    parser.add_argument("--ploidy", type=int, default=1, help="Organism ploidy for variant calling (default: haploid)")
    parser.add_argument(
        "--reference",
        choices=sorted(REFERENCE_REGISTRY.keys()),
        default="c_glabrata",
        help="Reference genome to use (default: c_glabrata)",
    )
    parser.add_argument("--amplicon", action="store_true", help="Produce per-position mutation frequency table")
    parser.add_argument("--deduplicate", action="store_true", help="Remove PCR duplicates using samtools markdup")
    parser.add_argument(
        "--annotate",
        action="store_true",
        help="Annotate VCF with gene information from GFF (requires GFF file for reference)",
    )
    parser.add_argument(
        "--skip-prefetch",
        action="store_true",
        help="Skip prefetch step when downloading SRA files (fasterq-dump will download if needed)",
    )
    parser.add_argument("--keep-intermediate", action="store_true", help="Retain SAM and BCF intermediates")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
    parser.add_argument("--force-reference", action="store_true", help="Redownload and reindex reference genome")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    configure_logging(args.verbose)

    # Validate arguments
    if not args.input and not args.bioproject:
        logging.error("Either --input or --bioproject must be provided")
        return 1
    
    if args.input and args.bioproject:
        logging.error("Cannot use both --input and --bioproject. Use one or the other.")
        return 1

    output_dir = Path(args.output).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    threads = args.threads or max(1, (os.cpu_count() or 2) - 1)

    ensure_environment()
    reference = ensure_reference(reference_key=args.reference, force=args.force_reference)

    # Handle BioProject download - process each SRA sequentially
    if args.bioproject:
        logging.info(f"Processing BioProject {args.bioproject}...")
        try:
            # Query BioProject and generate/update metadata table with study IDs
            # Returns study IDs, metadata file, and set of already processed run IDs
            study_ids, metadata_file, processed_runs = download_and_convert_bioproject(
                bioproject_id=args.bioproject,
                output_dir=output_dir,
                threads=threads,
                skip_prefetch=args.skip_prefetch,
            )
            
            if not study_ids:
                logging.error(f"No study IDs found for BioProject {args.bioproject}")
                return 1
            
            logging.info(f"Metadata table: {metadata_file}")
            logging.info(f"Processing {len(study_ids)} study/studies sequentially...")
            logging.info(f"Skipping {len(processed_runs)} already processed run(s)")
            
            # Process each study sequentially - query runs on-demand
            for study_id in study_ids:
                logging.info(f"Querying runs from study {study_id}...")
                
                # Query this study's runs on-demand (only if not already in metadata table)
                run_ids = query_study_runs(study_id, metadata_file=metadata_file)
                
                if not run_ids:
                    logging.warning(f"No runs found for study {study_id}, skipping...")
                    continue
                
                # Filter out already processed runs
                unprocessed_runs = [r for r in run_ids if r not in processed_runs]
                
                if not unprocessed_runs:
                    logging.info(f"All {len(run_ids)} run(s) from study {study_id} already processed, skipping...")
                    continue
                
                if len(unprocessed_runs) < len(run_ids):
                    logging.info(f"Study {study_id}: {len(unprocessed_runs)} unprocessed, {len(run_ids) - len(unprocessed_runs)} already processed")
                
                logging.info(f"Found {len(unprocessed_runs)} unprocessed run(s) in study {study_id}")
                
                # Process each unprocessed SRA run from this study
                for run_id in unprocessed_runs:
                    logging.info(f"Processing SRA run {run_id}...")
                    
                    # Create output subfolder for this SRA run
                    sra_output_dir = output_dir / run_id
                    sra_output_dir.mkdir(parents=True, exist_ok=True)
                    
                    # Create temporary directory for FASTQ files for this run
                    fastq_temp_dir = sra_output_dir / "fastq"
                    
                    try:
                        # Optionally prefetch (for better caching)
                        sra_file = None
                        if not args.skip_prefetch:
                            sra_file = prefetch_sra(run_id)
                        
                        # Convert to FASTQ
                        fastq_files = fasterq_dump(
                            run_id, 
                            fastq_temp_dir, 
                            sra_file=sra_file, 
                            threads=threads
                        )
                        
                        if not fastq_files:
                            logging.warning(f"No FASTQ files generated for {run_id}, skipping...")
                            continue
                        
                        # Group FASTQ files into samples
                        fastqs = discover_fastqs(fastq_temp_dir)
                        samples = group_samples(fastqs)
                        
                        # Process each sample from this SRA run
                        for sample in samples:
                            # Use the SRA run ID as the sample name prefix
                            original_name = sample.name
                            sample.name = f"{run_id}_{original_name}"
                            
                            align_sample(
                                sample,
                                reference,
                                sra_output_dir,  # Results go in SRA-specific subfolder
                                threads,
                                args.keep_intermediate,
                                args.ploidy,
                                args.amplicon,
                                args.deduplicate,
                                args.annotate,
                                args.reference,
                            )
                        
                    except Exception as e:
                        logging.error(f"Failed to process SRA run {run_id}: {e}")
                        logging.error(f"Continuing with next SRA run...")
                        continue
            
            logging.info(f"Completed processing BioProject {args.bioproject}")
            return 0
            
        except Exception as e:
            logging.error(f"Failed to process BioProject {args.bioproject}: {e}")
            return 1
    else:
        # Standard processing of local FASTQ files
        input_path = Path(args.input).expanduser().resolve()
    fastqs = discover_fastqs(input_path)
    samples = group_samples(fastqs)

    for sample in samples:
            align_sample(
                sample,
                reference,
                output_dir,
                threads,
                args.keep_intermediate,
                args.ploidy,
                args.amplicon,
                args.deduplicate,
                args.annotate,
                args.reference,
            )

    logging.info("Completed processing %d sample(s)", len(samples))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
