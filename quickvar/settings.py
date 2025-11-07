"""Shared configuration values for QuickVar."""

from pathlib import Path

CACHE_DIR = Path.home() / ".quickvar"
MAMBA_ROOT_PREFIX = CACHE_DIR / "micromamba"
ENV_NAME = "quickvar"
REFERENCE_DIR = CACHE_DIR / "reference"
DOWNLOAD_CHUNK_SIZE = 1_048_576  # 1 MB
