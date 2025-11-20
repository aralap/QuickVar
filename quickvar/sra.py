"""SRA download and conversion functionality for QuickVar."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import List, Optional

from .install import ensure_environment, micromamba_run
from .settings import CACHE_DIR

SRA_CACHE_DIR = CACHE_DIR / "sra"


def ensure_sra_cache() -> Path:
    """Ensure SRA cache directory exists and set environment variable."""
    SRA_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    # Set NCBI SRA cache directory
    os.environ["NCBI_SRA_CACHE"] = str(SRA_CACHE_DIR)
    return SRA_CACHE_DIR


def query_bioproject(bioproject_id: str) -> List[str]:
    """
    Query BioProject ID to get list of SRA run IDs using NCBI Entrez E-utilities API.
    
    Args:
        bioproject_id: NCBI BioProject ID (e.g., "PRJNA123456")
    
    Returns:
        List of SRA run IDs (e.g., ["SRR123456", "SRR123457"])
    """
    import urllib.request
    import urllib.parse
    import xml.etree.ElementTree as ET
    import time
    
    logging.info(f"Querying BioProject {bioproject_id} for SRA runs...")
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    
    try:
        # Step 1: Search for the BioProject in BioProject database
        search_params = {
            "db": "bioproject",
            "term": bioproject_id,
            "retmode": "xml",
        }
        search_url = f"{base_url}esearch.fcgi?{urllib.parse.urlencode(search_params)}"
        
        with urllib.request.urlopen(search_url) as response:
            search_xml = response.read().decode()
            root = ET.fromstring(search_xml)
            bioproject_ids = [id_elem.text for id_elem in root.findall(".//Id")]
            
            if not bioproject_ids:
                logging.warning(f"BioProject {bioproject_id} not found in NCBI")
                return []
        
        # Small delay to respect NCBI rate limits
        time.sleep(0.34)
        
        # Step 2: Link BioProject to SRA studies (SRP)
        link_params = {
            "dbfrom": "bioproject",
            "db": "sra",
            "id": bioproject_ids[0],
            "retmode": "xml",
        }
        link_url = f"{base_url}elink.fcgi?{urllib.parse.urlencode(link_params)}"
        
        with urllib.request.urlopen(link_url) as response:
            link_xml = response.read().decode()
            link_root = ET.fromstring(link_xml)
            sra_study_ids = [id_elem.text for id_elem in link_root.findall(".//Id")]
            
            if not sra_study_ids:
                logging.warning(f"No SRA studies found for BioProject {bioproject_id}")
                return []
        
        logging.info(f"Found {len(sra_study_ids)} SRA study/studies for BioProject {bioproject_id}")
        
        # Step 3: Get runs from each SRP study using efetch
        # First get the study record ID, then fetch detailed XML containing run accessions
        all_run_ids = []
        for srp in sra_study_ids:
            time.sleep(0.34)  # Rate limiting
            
            try:
                # First, search for the study to get its internal record ID
                study_search_params = {
                    "db": "sra",
                    "term": srp,
                    "retmode": "xml",
                }
                study_search_url = f"{base_url}esearch.fcgi?{urllib.parse.urlencode(study_search_params)}"
                
                with urllib.request.urlopen(study_search_url) as response:
                    study_search_xml = response.read().decode()
                    study_search_root = ET.fromstring(study_search_xml)
                    study_record_ids = [id_elem.text for id_elem in study_search_root.findall(".//Id")]
                    
                    if not study_record_ids:
                        logging.debug(f"Could not find record ID for study {srp}")
                        continue
                
                time.sleep(0.34)  # Rate limiting
                
                # Use efetch to get detailed record XML, which includes run accessions
                fetch_params = {
                    "db": "sra",
                    "id": study_record_ids[0],
                    "retmode": "xml",
                }
                fetch_url = f"{base_url}efetch.fcgi?{urllib.parse.urlencode(fetch_params)}"
                
                with urllib.request.urlopen(fetch_url) as response:
                    fetch_xml = response.read().decode()
                    fetch_root = ET.fromstring(fetch_xml)
                    
                    # Extract run accessions from RUN elements
                    # RUN elements have accession attribute with SRR* values
                    runs = fetch_root.findall(".//RUN")
                    for run in runs:
                        accession = run.get("accession")
                        if accession and accession.startswith("SRR"):
                            all_run_ids.append(accession)
                    
                    logging.debug(f"Found {len([r for r in all_run_ids if r.startswith('SRR') and r not in [x for x in all_run_ids[:-1] if x.startswith('SRR')]])} runs from study {srp}")
                        
            except Exception as e:
                logging.debug(f"Failed to get runs from study {srp}: {e}")
                continue
        
        if not all_run_ids:
            logging.warning(f"No SRA runs found for BioProject {bioproject_id}")
            return []
        
        # Step 4: Filter and clean run IDs (keep only SRR* accessions)
        run_ids = [str(run_id).strip() for run_id in all_run_ids if run_id and str(run_id).startswith("SRR")]
        run_ids = list(set(run_ids))  # Remove duplicates
        
        logging.info(f"Found {len(run_ids)} unique SRA run(s) for BioProject {bioproject_id}")
        return run_ids
        
    except urllib.error.HTTPError as e:
        logging.error(f"HTTP error querying BioProject {bioproject_id}: {e}")
        raise
    except ET.ParseError as e:
        logging.error(f"XML parse error querying BioProject {bioproject_id}: {e}")
        raise
    except Exception as e:
        logging.error(f"Failed to query BioProject {bioproject_id}: {e}")
        raise


def prefetch_sra(run_id: str, output_dir: Optional[Path] = None) -> Optional[Path]:
    """
    Download SRA file using prefetch.
    
    Args:
        run_id: SRA run ID (e.g., "SRR123456")
        output_dir: Optional directory to save SRA file (default: SRA cache)
    
    Returns:
        Path to downloaded SRA file, or None if failed
    """
    ensure_sra_cache()
    ensure_environment()
    
    if output_dir is None:
        output_dir = SRA_CACHE_DIR
    
    logging.info(f"Downloading SRA file for {run_id}...")
    
    try:
        result = micromamba_run(
            [
                "prefetch",
                "--output-directory",
                str(output_dir),
                run_id,
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        
        if result.returncode != 0:
            error_msg = result.stderr or result.stdout or "unknown error"
            logging.warning(f"prefetch failed for {run_id}: {error_msg}")
            return None
        
        # prefetch creates a directory with the run_id, containing the .sra file
        sra_dir = output_dir / run_id
        sra_file = sra_dir / f"{run_id}.sra"
        
        if sra_file.exists():
            logging.info(f"Downloaded SRA file: {sra_file}")
            return sra_file
        else:
            # Sometimes prefetch creates the file directly
            direct_file = output_dir / f"{run_id}.sra"
            if direct_file.exists():
                logging.info(f"Downloaded SRA file: {direct_file}")
                return direct_file
            else:
                logging.warning(f"SRA file not found after prefetch for {run_id}")
                return None
                
    except Exception as e:
        logging.warning(f"Failed to prefetch {run_id}: {e}")
        return None


def fasterq_dump(
    run_id: str,
    output_dir: Path,
    sra_file: Optional[Path] = None,
    threads: int = 1,
) -> List[Path]:
    """
    Convert SRA file to FASTQ using fasterq-dump.
    
    Args:
        run_id: SRA run ID (e.g., "SRR123456")
        output_dir: Directory to save FASTQ files
        sra_file: Optional path to SRA file (if None, will use cache or download)
        threads: Number of threads for fasterq-dump
    
    Returns:
        List of paths to generated FASTQ files
    """
    ensure_environment()
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logging.info(f"Converting {run_id} to FASTQ...")
    
    # If SRA file not provided, try to find it in cache or download
    if sra_file is None or not sra_file.exists():
        ensure_sra_cache()
        # Check cache first - prefetch creates a directory structure
        sra_dir = SRA_CACHE_DIR / run_id
        if sra_dir.exists():
            # Look for .sra file in the directory
            sra_files = list(sra_dir.glob("*.sra"))
            if sra_files:
                sra_file = sra_files[0]
            else:
                sra_file = sra_dir / f"{run_id}.sra"
        else:
            sra_file = SRA_CACHE_DIR / f"{run_id}.sra"
        
        if not sra_file.exists():
            # Download it
            sra_file = prefetch_sra(run_id)
            if sra_file is None:
                logging.warning(f"Could not find or download SRA file for {run_id}, fasterq-dump will attempt to download")
                sra_file = None  # Let fasterq-dump handle it
    
    try:
        # fasterq-dump arguments
        args = [
            "fasterq-dump",
            "--outdir",
            str(output_dir),
            "--threads",
            str(threads),
            "--skip-technical",
            "--split-files",  # Split paired-end reads
        ]
        
        # Use SRA file path if available, otherwise use run_id (fasterq-dump will download it)
        if sra_file and sra_file.exists():
            args.append(str(sra_file))
        else:
            # fasterq-dump can download directly if we provide the run_id
            args.append(run_id)
        
        result = micromamba_run(
            args,
            check=False,
            capture_output=True,
            text=True,
        )
        
        if result.returncode != 0:
            error_msg = result.stderr or result.stdout or "unknown error"
            logging.warning(f"fasterq-dump failed for {run_id}: {error_msg}")
            return []
        
        # Find generated FASTQ files
        fastq_files = []
        for pattern in [f"{run_id}_1.fastq", f"{run_id}_2.fastq", f"{run_id}.fastq"]:
            fastq_path = output_dir / pattern
            if fastq_path.exists():
                fastq_files.append(fastq_path)
        
        if not fastq_files:
            # Try without underscore
            for pattern in [f"{run_id}.fastq", f"{run_id}_1.fastq", f"{run_id}_2.fastq"]:
                fastq_path = output_dir / pattern
                if fastq_path.exists():
                    fastq_files.append(fastq_path)
        
        if fastq_files:
            logging.info(f"Generated {len(fastq_files)} FASTQ file(s) for {run_id}")
        else:
            logging.warning(f"No FASTQ files found after fasterq-dump for {run_id}")
        
        return sorted(fastq_files)
        
    except Exception as e:
        logging.warning(f"Failed to convert {run_id} to FASTQ: {e}")
        return []


def download_and_convert_bioproject(
    bioproject_id: str,
    output_dir: Path,
    threads: int = 1,
    skip_prefetch: bool = False,
) -> List[Path]:
    """
    Download and convert all SRA runs from a BioProject to FASTQ.
    
    Args:
        bioproject_id: NCBI BioProject ID (e.g., "PRJNA123456")
        output_dir: Directory to save FASTQ files
        threads: Number of threads for fasterq-dump
        skip_prefetch: If True, skip prefetch step (fasterq-dump will download if needed)
    
    Returns:
        List of paths to all generated FASTQ files
    """
    # Query BioProject for SRA run IDs
    run_ids = query_bioproject(bioproject_id)
    
    if not run_ids:
        logging.warning(f"No SRA runs found for BioProject {bioproject_id}")
        return []
    
    all_fastq_files = []
    
    for run_id in run_ids:
        logging.info(f"Processing {run_id}...")
        
        # Optionally prefetch (for better caching)
        sra_file = None
        if not skip_prefetch:
            sra_file = prefetch_sra(run_id)
        
        # Convert to FASTQ
        fastq_files = fasterq_dump(run_id, output_dir, sra_file=sra_file, threads=threads)
        all_fastq_files.extend(fastq_files)
    
    logging.info(f"Downloaded and converted {len(all_fastq_files)} FASTQ file(s) from BioProject {bioproject_id}")
    return all_fastq_files

