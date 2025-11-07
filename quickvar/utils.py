"""Utility helpers for QuickVar."""

from __future__ import annotations

import gzip
import shutil
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import IO, Optional

from .settings import DOWNLOAD_CHUNK_SIZE


class DownloadError(RuntimeError):
    """Raised when a download fails."""


def _format_progress_downloaded(downloaded: int, total: int | None, elapsed: float) -> str:
    mb_downloaded = downloaded / 1_000_000
    speed = mb_downloaded / elapsed if elapsed > 0 else 0.0
    if total:
        mb_total = total / 1_000_000
        percent = downloaded / total * 100
        return f"{percent:5.1f}% ({mb_downloaded:.1f}/{mb_total:.1f} MB) {speed:.2f} MB/s"
    return f"{mb_downloaded:.1f} MB {speed:.2f} MB/s"


def stream_download(
    url: str,
    destination: Path,
    *,
    desc: Optional[str] = None,
    response: Optional[IO[bytes]] = None,
) -> Path:
    """Download *url* to *destination* with a simple progress indicator."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    close_response = False
    try:
        if response is None:
            response = urllib.request.urlopen(url)
            close_response = True
        total_header = response.getheader("Content-Length") if hasattr(response, "getheader") else None
        total = int(total_header) if total_header else None
        tmp_path = destination.with_suffix(destination.suffix + ".partial")
        downloaded = 0
        start_time = time.time()
        label = desc or destination.name
        with open(tmp_path, "wb") as handle:
            while True:
                chunk = response.read(DOWNLOAD_CHUNK_SIZE)
                if not chunk:
                    break
                handle.write(chunk)
                downloaded += len(chunk)
                status = _format_progress_downloaded(downloaded, total, time.time() - start_time)
                sys.stdout.write(f"\r{label}: {status}")
                sys.stdout.flush()
        sys.stdout.write("\n")
        tmp_path.replace(destination)
        return destination
    except (urllib.error.HTTPError, urllib.error.URLError) as exc:  # pragma: no cover - network failures
        raise DownloadError(f"Failed to download {url}: {exc}") from exc
    finally:
        if close_response and response is not None:
            response.close()


def gunzip_file(src: Path, dest: Path) -> Path:
    """Decompress a GZip file from *src* into *dest*."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(src, "rb") as gz_in, open(dest, "wb") as out:
        shutil.copyfileobj(gz_in, out)
    return dest


def ensure_executable(path: Path) -> Path:
    """Make sure *path* is executable on POSIX systems."""
    if path.exists():
        try:
            current_mode = path.stat().st_mode
            path.chmod(current_mode | 0o111)
        except OSError:
            pass
    return path
