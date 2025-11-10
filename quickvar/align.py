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
) -> None:
    logging.info("Processing sample %s", sample.name)
    sample_dir = output_dir / sample.name
    sample_dir.mkdir(parents=True, exist_ok=True)
    sam_path = sample_dir / f"{sample.name}.sam"
    bam_path = sample_dir / f"{sample.name}.sorted.bam"
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

    micromamba_run(["samtools", "index", str(bam_path)])

    bcftools_mpileup = [
        "bcftools",
        "mpileup",
        "-Ob",
        "-f",
        str(reference["fasta"]),
        "-o",
        str(bcf_path),
        str(bam_path),
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

    if amplicon:
        generate_amplicon_report(
            sample=sample,
            reference=reference,
            bam_path=bam_path,
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
    with open(summary_path, "w", encoding="utf-8") as handle:
        handle.write(
            "chrom\tpos\tref_base\talt_base\talt_count\tdepth\tfrequency\tmutation\n"
        )
        for line in result.stdout.splitlines():
            fields = line.strip().split("\t")
            if len(fields) < 6:
                continue
            chrom, pos, ref_base, depth_str, bases, _quals = fields[:6]
            try:
                depth = int(depth_str)
            except ValueError:
                continue
            if depth <= 0:
                continue
            counts = parse_pileup_bases(bases, ref_base)
            total_depth = sum(counts.values())
            if total_depth == 0:
                continue
            ref_upper = ref_base.upper()
            for base, count in counts.items():
                if base == ref_upper or count == 0:
                    continue
                frequency = count / total_depth if total_depth else 0.0
                if base.startswith("+"):
                    mutation = f"{ref_upper}>+{base[1:]}"
                elif base.startswith("-"):
                    mutation = f"{ref_upper}>-{base[1:]}"
                else:
                    mutation = f"{ref_upper}>{base}"
                handle.write(
                    f"{chrom}\t{pos}\t{ref_upper}\t{base}\t{count}\t{total_depth}\t{frequency:.4f}\t{mutation}\n"
                )


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Align FASTQ reads and call variants")
    parser.add_argument("--input", required=True, help="Path to FASTQ file or directory containing FASTQs")
    parser.add_argument("--output", default=DEFAULT_OUTPUT_NAME, help="Directory to write results (default: Results)")
    parser.add_argument("--threads", type=int, default=0, help="Number of CPU threads (default: auto)")
    parser.add_argument("--ploidy", type=int, default=1, help="Organism ploidy for variant calling (default: haploid)")
    parser.add_argument("--amplicon", action="store_true", help="Produce per-position mutation frequency table")
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
    reference = ensure_reference(force=args.force_reference)

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
        )

    logging.info("Completed processing %d sample(s)", len(samples))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
