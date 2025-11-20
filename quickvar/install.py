"""Environment bootstrap logic for QuickVar."""

from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path
from typing import Iterable

import urllib.request
from urllib.parse import urlparse

from . import __version__
from .settings import CACHE_DIR, ENV_NAME, MAMBA_ROOT_PREFIX
from .utils import ensure_executable, stream_download

MICROMAMBA_API = "https://micro.mamba.pm/api/micromamba"
MICROMAMBA_DIR = CACHE_DIR / "micromamba-bin"
MICROMAMBA_BINARY_NAME = "micromamba.exe" if platform.system().lower().startswith("win") else "micromamba"
QUICKVAR_PACKAGES = [
    "python=3.11",
    "minimap2>=2.24",
    "samtools>=1.19",
    "bcftools>=1.19",
    "sra-tools>=3.0",
]


def ensure_cache_dirs() -> None:
    """Create cache directories used by QuickVar."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    MAMBA_ROOT_PREFIX.mkdir(parents=True, exist_ok=True)


def detect_platform_tag() -> str:
    """Return the micromamba platform tag for the current system."""
    system = platform.system().lower()
    machine = platform.machine().lower()
    if system == "darwin":
        arch = "arm64" if "arm" in machine or "aarch64" in machine else "64"
        return f"osx-{arch}"
    if system == "linux":
        if any(token in machine for token in ("aarch64", "arm64")):
            return "linux-aarch64"
        return "linux-64"
    if system == "windows":
        return "win-64"
    raise RuntimeError(f"Unsupported platform: {system} {machine}")


def resolve_suffix_from_url(url: str) -> str:
    """Determine archive suffix from a URL, ignoring query parameters."""
    parsed = urlparse(url)
    suffixes = Path(parsed.path).suffixes
    if suffixes[-2:] == [".tar", ".bz2"]:
        suffix = ".tar.bz2"
    elif suffixes[-1:] == [".zip"]:
        suffix = ".zip"
    else:
        suffix = suffixes[-1] if suffixes else ".tar.bz2"
    return suffix


def download_micromamba() -> Path:
    """Download and extract micromamba, returning the binary path."""
    ensure_cache_dirs()
    if MICROMAMBA_DIR.exists():
        candidate = MICROMAMBA_DIR / MICROMAMBA_BINARY_NAME
        if candidate.exists():
            return ensure_executable(candidate)
    platform_tag = detect_platform_tag()
    url = f"{MICROMAMBA_API}/{platform_tag}/latest"
    response = urllib.request.urlopen(url)
    final_url = response.geturl()
    suffix = resolve_suffix_from_url(final_url)
    archive_path = CACHE_DIR / f"micromamba{suffix}"
    if archive_path.exists():
        archive_path.unlink()
    stream_download(url, archive_path, desc="Micromamba", response=response)
    response.close()
    if MICROMAMBA_DIR.exists():
        shutil.rmtree(MICROMAMBA_DIR)
    MICROMAMBA_DIR.mkdir(parents=True, exist_ok=True)
    if suffix.endswith(".tar.bz2"):
        with tarfile.open(archive_path, "r:bz2") as tar:
            tar.extractall(MICROMAMBA_DIR)
    elif suffix.endswith(".zip"):
        with zipfile.ZipFile(archive_path, "r") as zf:
            zf.extractall(MICROMAMBA_DIR)
    else:
        raise RuntimeError(f"Unsupported micromamba archive format: {suffix}")
    # micromamba is typically under bin/
    for candidate in MICROMAMBA_DIR.rglob(MICROMAMBA_BINARY_NAME):
        return ensure_executable(candidate)
    raise RuntimeError("Micromamba binary not found after extraction")


def get_micromamba_binary() -> Path:
    """Get the micromamba executable path, downloading it if necessary."""
    candidate = MICROMAMBA_DIR / MICROMAMBA_BINARY_NAME
    if candidate.exists():
        return ensure_executable(candidate)
    return download_micromamba()


def micromamba_env_exists() -> bool:
    """Check whether the QuickVar environment already exists."""
    env_dir = MAMBA_ROOT_PREFIX / "envs" / ENV_NAME
    return env_dir.exists() and any(env_dir.iterdir())


def run_micromamba(
    args: Iterable[str],
    *,
    check: bool = True,
    capture_output: bool = False,
    text: bool = False,
) -> subprocess.CompletedProcess:
    """Execute micromamba with the provided arguments."""
    binary = get_micromamba_binary()
    env = os.environ.copy()
    env["MAMBA_ROOT_PREFIX"] = str(MAMBA_ROOT_PREFIX)
    command = [str(binary), *args]
    return subprocess.run(
        command,
        check=check,
        env=env,
        capture_output=capture_output,
        text=text,
    )


def create_environment(force: bool = False) -> None:
    """Create the QuickVar micromamba environment."""
    if force and micromamba_env_exists():
        run_micromamba(["env", "remove", "-y", "-n", ENV_NAME], check=False)
    if micromamba_env_exists():
        return
    run_micromamba(
        [
            "create",
            "-y",
            "-n",
            ENV_NAME,
            "-c",
            "conda-forge",
            "-c",
            "bioconda",
            *QUICKVAR_PACKAGES,
        ]
    )


def remove_environment() -> None:
    """Remove the QuickVar environment and downloaded binaries."""
    if micromamba_env_exists():
        run_micromamba(["env", "remove", "-y", "-n", ENV_NAME], check=False)
    # Preserve micromamba binary unless --purge


def ensure_environment(force: bool = False) -> Path:
    """Ensure the QuickVar environment exists and return its path."""
    ensure_cache_dirs()
    create_environment(force=force)
    return MAMBA_ROOT_PREFIX / "envs" / ENV_NAME


def micromamba_run(
    command: Iterable[str],
    *,
    check: bool = True,
    capture_output: bool = False,
    text: bool = False,
) -> subprocess.CompletedProcess:
    """Run a command within the QuickVar environment."""
    run_args = ["run", "-n", ENV_NAME, *command]
    return run_micromamba(
        run_args,
        check=check,
        capture_output=capture_output,
        text=text,
    )


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Install QuickVar dependencies")
    parser.add_argument("--force", action="store_true", help="Recreate the environment")
    parser.add_argument("--remove", action="store_true", help="Remove the environment")
    parser.add_argument("--purge", action="store_true", help="Remove the environment and micromamba binary")
    parser.add_argument("--show-path", action="store_true", help="Print the environment path")
    parser.add_argument("--version", action="store_true", help="Print QuickVar version")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    if args.version:
        print(__version__)
        return 0
    if args.remove or args.purge:
        remove_environment()
        if args.purge and MICROMAMBA_DIR.exists():
            shutil.rmtree(MICROMAMBA_DIR)
        return 0
    env_path = ensure_environment(force=args.force)
    
    if args.show_path:
        print(env_path)
    else:
        print(f"QuickVar environment available at {env_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
