"""Reference genome management for supported Candida species."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Dict

from .install import ensure_environment, micromamba_run
from .settings import REFERENCE_DIR, REFERENCE_REGISTRY
from .utils import gunzip_file, stream_download

LOCAL_REFERENCE_DIR = Path(__file__).resolve().parents[1] / "ref_genome"


def ensure_reference(reference_key: str = "c_glabrata", force: bool = False) -> Dict[str, Path]:
    """Ensure the reference genome and indexes are available."""
    if reference_key not in REFERENCE_REGISTRY:
        raise ValueError(
            f"Unsupported reference '{reference_key}'. "
            f"Available options: {', '.join(sorted(REFERENCE_REGISTRY))}"
        )
    metadata = REFERENCE_REGISTRY[reference_key]
    local_gz = metadata.get("local_gz")
    fasta_name = metadata["fasta_name"]
    mmi_name = metadata["mmi_name"]

    fasta_gz_path = REFERENCE_DIR / f"{fasta_name}.gz"
    fasta_path = REFERENCE_DIR / fasta_name
    mmi_path = REFERENCE_DIR / mmi_name

    ensure_environment()
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    if force:
        for path in (fasta_gz_path, fasta_path, mmi_path):
            if path.exists():
                path.unlink()
    if force or not fasta_path.exists():
        if local_gz and local_gz.exists():
            shutil.copyfile(local_gz, fasta_gz_path)
        else:
            stream_download(metadata["url"], fasta_gz_path, desc=f"{reference_key} reference")
        gunzip_file(fasta_gz_path, fasta_path)
    if force or not mmi_path.exists():
        micromamba_run(["minimap2", "-d", str(mmi_path), str(fasta_path)])
    fai_path = fasta_path.with_suffix(".fai")
    if force or not fai_path.exists():
        micromamba_run(["samtools", "faidx", str(fasta_path)])
    return {
        "fasta": fasta_path,
        "mmi": mmi_path,
        "fai": fai_path,
    }
