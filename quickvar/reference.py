"""Reference genome management for *Candida glabrata*."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Dict

from .install import ensure_environment, micromamba_run
from .settings import REFERENCE_DIR
from .utils import gunzip_file, stream_download

REFERENCE_FASTA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/545/"
    "GCF_000002545.3_ASM254v2/GCF_000002545.3_ASM254v2_genomic.fna.gz"
)
REFERENCE_FASTA_NAME = "c_glabrata_cbs138.fna"
REFERENCE_FASTA_GZ = REFERENCE_DIR / f"{REFERENCE_FASTA_NAME}.gz"
REFERENCE_FASTA = REFERENCE_DIR / REFERENCE_FASTA_NAME
REFERENCE_MMI = REFERENCE_DIR / "c_glabrata_cbs138.mmi"
REFERENCE_DICT = REFERENCE_DIR / "c_glabrata_cbs138.dict"

LOCAL_REFERENCE_DIR = Path(__file__).resolve().parents[1] / "ref_genome"
LOCAL_REFERENCE_FASTA_GZ = LOCAL_REFERENCE_DIR / "C_glabrata_CBS138_current_chromosomes.fasta.gz"


def ensure_reference(force: bool = False) -> Dict[str, Path]:
    """Ensure the reference genome and indexes are available."""
    ensure_environment()
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    if force:
        if REFERENCE_FASTA_GZ.exists():
            REFERENCE_FASTA_GZ.unlink()
        if REFERENCE_FASTA.exists():
            REFERENCE_FASTA.unlink()
        if REFERENCE_MMI.exists():
            REFERENCE_MMI.unlink()
    if force or not REFERENCE_FASTA.exists():
        if LOCAL_REFERENCE_FASTA_GZ.exists():
            shutil.copyfile(LOCAL_REFERENCE_FASTA_GZ, REFERENCE_FASTA_GZ)
        else:
            stream_download(REFERENCE_FASTA_URL, REFERENCE_FASTA_GZ, desc="C. glabrata reference")
        gunzip_file(REFERENCE_FASTA_GZ, REFERENCE_FASTA)
    if force or not REFERENCE_MMI.exists():
        micromamba_run(["minimap2", "-d", str(REFERENCE_MMI), str(REFERENCE_FASTA)])
    fai_path = REFERENCE_FASTA.with_suffix(".fai")
    if force or not fai_path.exists():
        micromamba_run(["samtools", "faidx", str(REFERENCE_FASTA)])
    return {
        "fasta": REFERENCE_FASTA,
        "mmi": REFERENCE_MMI,
        "fai": fai_path,
    }
