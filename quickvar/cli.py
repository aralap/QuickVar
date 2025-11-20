"""Command-line interface entry point for QuickVar."""

from __future__ import annotations

import argparse
import sys

from . import __version__
from . import align as align_module
from . import install as install_module
from .settings import REFERENCE_REGISTRY


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="quickvar", description="QuickVar CLI")
    parser.add_argument("--version", action="store_true", help="Show QuickVar version and exit")
    subparsers = parser.add_subparsers(dest="command")

    install_parser = subparsers.add_parser("install", help="Install pipeline dependencies")
    install_parser.add_argument("--force", action="store_true", help="Recreate the environment")
    install_parser.add_argument("--remove", action="store_true", help="Remove the environment")
    install_parser.add_argument("--purge", action="store_true", help="Remove environment and micromamba")
    install_parser.add_argument("--show-path", action="store_true", help="Print environment path")

    align_parser = subparsers.add_parser("align", help="Align FASTQ reads and call variants")
    align_parser.add_argument("--input", help="Path to FASTQ file or directory (required if --bioproject not used)")
    align_parser.add_argument(
        "--bioproject",
        help="NCBI BioProject ID (e.g., PRJNA123456) to download and process SRA files",
    )
    align_parser.add_argument("--output", default=align_module.DEFAULT_OUTPUT_NAME, help="Output directory")
    align_parser.add_argument("--threads", type=int, default=0, help="Number of CPU threads")
    align_parser.add_argument("--ploidy", type=int, default=1, help="Organism ploidy (default: 1)")
    align_parser.add_argument(
        "--reference",
        choices=sorted(REFERENCE_REGISTRY.keys()),
        default="c_glabrata",
        help="Reference genome identifier (default: c_glabrata)",
    )
    align_parser.add_argument("--amplicon", action="store_true", help="Generate per-position amplicon summary")
    align_parser.add_argument("--deduplicate", action="store_true", help="Remove PCR duplicates with samtools markdup")
    align_parser.add_argument(
        "--annotate",
        action="store_true",
        help="Annotate VCF with gene information from GFF (requires GFF file for reference)",
    )
    align_parser.add_argument(
        "--skip-prefetch",
        action="store_true",
        help="Skip prefetch step when downloading SRA files (fasterq-dump will download if needed)",
    )
    align_parser.add_argument("--keep-intermediate", action="store_true", help="Keep SAM and BCF files")
    align_parser.add_argument("--verbose", action="store_true", help="Verbose logging")
    align_parser.add_argument("--force-reference", action="store_true", help="Redownload reference data")

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.version:
        print(__version__)
        return 0

    if args.command == "install":
        install_args = []
        if args.force:
            install_args.append("--force")
        if args.remove:
            install_args.append("--remove")
        if args.purge:
            install_args.append("--purge")
        if args.show_path:
            install_args.append("--show-path")
        return install_module.main(install_args)

    if args.command == "align":
        align_args = [
            "--output",
            args.output,
            "--threads",
            str(args.threads),
            "--ploidy",
            str(args.ploidy),
        ]
        if args.input:
            align_args.extend(["--input", args.input])
        if args.bioproject:
            align_args.extend(["--bioproject", args.bioproject])
        if args.amplicon:
            align_args.append("--amplicon")
        if args.deduplicate:
            align_args.append("--deduplicate")
        if args.annotate:
            align_args.append("--annotate")
        if args.skip_prefetch:
            align_args.append("--skip-prefetch")
        align_args.extend(["--reference", args.reference])
        if args.keep_intermediate:
            align_args.append("--keep-intermediate")
        if args.verbose:
            align_args.append("--verbose")
        if args.force_reference:
            align_args.append("--force-reference")
        return align_module.main(align_args)

    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
