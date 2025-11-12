"""Shared configuration values for QuickVar."""

from pathlib import Path

CACHE_DIR = Path.home() / ".quickvar"
MAMBA_ROOT_PREFIX = CACHE_DIR / "micromamba"
ENV_NAME = "quickvar"
REFERENCE_DIR = CACHE_DIR / "reference"
DOWNLOAD_CHUNK_SIZE = 1_048_576  # 1 MB

LOCAL_REFERENCE_DIR = Path(__file__).resolve().parents[1] / "ref_genome"

REFERENCE_REGISTRY = {
    "c_glabrata": {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/545/"
            "GCF_000002545.3_ASM254v2/GCF_000002545.3_ASM254v2_genomic.fna.gz"
        ),
        "fasta_name": "c_glabrata_cbs138.fna",
        "mmi_name": "c_glabrata_cbs138.mmi",
        "local_gz": LOCAL_REFERENCE_DIR / "C_glabrata_CBS138_current_chromosomes.fasta.gz",
    },
    "c_auris": {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/775/015/"
            "GCF_002775015.1_Cand_auris_B11221_V1/GCF_002775015.1_Cand_auris_B11221_V1_genomic.fna.gz"
        ),
        "fasta_name": "c_auris_b11221.fna",
        "mmi_name": "c_auris_b11221.mmi",
        "local_gz": LOCAL_REFERENCE_DIR / "GCF_002775015.1_Cand_auris_B11221_V1_genomic.fna.gz",
    },
}
