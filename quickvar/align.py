"""FASTQ alignment and variant calling for QuickVar."""

from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List

import os

from .install import ensure_environment, micromamba_run
from .reference import ensure_reference
from .settings import REFERENCE_REGISTRY

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


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Align FASTQ reads and call variants")
    parser.add_argument("--input", required=True, help="Path to FASTQ file or directory containing FASTQs")
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
    parser.add_argument("--keep-intermediate", action="store_true", help="Retain SAM and BCF intermediates")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
    parser.add_argument("--force-reference", action="store_true", help="Redownload and reindex reference genome")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    configure_logging(args.verbose)

    input_path = Path(args.input).expanduser().resolve()
    output_dir = Path(args.output).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    threads = args.threads or max(1, (os.cpu_count() or 2) - 1)

    ensure_environment()
    reference = ensure_reference(reference_key=args.reference, force=args.force_reference)

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
        )

    logging.info("Completed processing %d sample(s)", len(samples))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
