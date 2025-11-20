"""SRA download and conversion functionality for QuickVar."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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
        # Stream output to show native progress from prefetch
        result = micromamba_run(
            [
                "prefetch",
                "--output-directory",
                str(output_dir),
                run_id,
            ],
            check=False,
            capture_output=False,  # Stream output to show progress
            text=True,
        )
        
        if result.returncode != 0:
            logging.warning(f"prefetch failed for {run_id} (exit code {result.returncode})")
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
    Checks for existing FASTQ files first and skips conversion if found.
    
    Args:
        run_id: SRA run ID (e.g., "SRR123456")
        output_dir: Directory to save FASTQ files
        sra_file: Optional path to SRA file (if None, will use cache or download)
        threads: Number of threads for fasterq-dump
    
    Returns:
        List of paths to generated FASTQ files (existing or newly created)
    """
    ensure_environment()
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check for existing FASTQ files first
    existing_fastq_files = []
    for pattern in [f"{run_id}_1.fastq", f"{run_id}_2.fastq", f"{run_id}.fastq"]:
        fastq_path = output_dir / pattern
        if fastq_path.exists():
            existing_fastq_files.append(fastq_path)
    
    # Also check patterns without underscore
    if not existing_fastq_files:
        for pattern in [f"{run_id}.fastq", f"{run_id}_1.fastq", f"{run_id}_2.fastq"]:
            fastq_path = output_dir / pattern
            if fastq_path.exists():
                existing_fastq_files.append(fastq_path)
    
    if existing_fastq_files:
        logging.info(f"Found existing FASTQ file(s) for {run_id}, skipping conversion")
        return sorted(existing_fastq_files)
    
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
        
        # Stream output to show native progress from fasterq-dump
        result = micromamba_run(
            args,
            check=False,
            capture_output=False,  # Stream output to show progress
            text=True,
        )
        
        if result.returncode != 0:
            logging.warning(f"fasterq-dump failed for {run_id} (exit code {result.returncode})")
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


def query_bioproject_with_metadata(
    bioproject_id: str, 
    metadata_file: Optional[Path] = None
) -> tuple[List[str], List[Dict[str, str]]]:
    """
    Query BioProject ID to get list of SRA run IDs and their metadata.
    
    Args:
        bioproject_id: NCBI BioProject ID (e.g., "PRJNA123456")
    
    Returns:
        Tuple of (list of SRA run IDs, list of metadata dictionaries)
    """
    import urllib.request
    import urllib.parse
    import xml.etree.ElementTree as ET
    import time
    
    logging.info(f"Querying BioProject {bioproject_id} for SRA runs with metadata...")
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
                return [], []
        
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
                return [], []
        
        logging.info(f"Found {len(sra_study_ids)} SRA study/studies for BioProject {bioproject_id}")
        
        # Save the study IDs from XML response as metadata table FIRST (before processing individual studies)
        if metadata_file:
            metadata_file.parent.mkdir(parents=True, exist_ok=True)
            header_fields = ["study_accession", "bioproject_id"]
            with open(metadata_file, "w", encoding="utf-8") as f:
                f.write("\t".join(header_fields) + "\n")
                for study_id in sra_study_ids:
                    f.write(f"{study_id}\t{bioproject_id}\n")
            logging.info(f"Saved SRA study metadata table: {metadata_file} ({len(sra_study_ids)} studies)")
        
        # No bulk querying of runs - they will be queried on-demand during processing
        # Return empty run_ids and metadata - actual run processing happens later
        return [], []
        
    except urllib.error.HTTPError as e:
        logging.error(f"HTTP error querying BioProject {bioproject_id}: {e}")
        raise
    except ET.ParseError as e:
        logging.error(f"XML parse error querying BioProject {bioproject_id}: {e}")
        raise
    except Exception as e:
        logging.error(f"Failed to query BioProject {bioproject_id}: {e}")
        raise


def save_sra_metadata(metadata: List[Dict[str, str]], output_file: Path) -> None:
    """Save SRA metadata to a TSV file. Creates file even if metadata is empty (header only)."""
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    header_fields = ["run_accession", "study_accession", "study_title", "total_spots", "total_bases", "size"]
    
    with open(output_file, "w", encoding="utf-8") as f:
        # Write header
        f.write("\t".join(header_fields) + "\n")
        
        # Write data
        if metadata:
            for meta in metadata:
                row = [str(meta.get(field, "")) for field in header_fields]
                f.write("\t".join(row) + "\n")
    
    if metadata:
        logging.info(f"Saved SRA metadata to {output_file} ({len(metadata)} entries)")
    else:
        logging.info(f"Created empty metadata table: {output_file}")


def read_metadata_table(metadata_file: Path) -> tuple[List[str], set[str]]:
    """Read study IDs and run IDs from metadata table.
    
    Returns:
        Tuple of (list of study IDs, set of run IDs)
    """
    study_ids = []
    run_ids = set()
    
    if not metadata_file.exists():
        return study_ids, run_ids
    
    with open(metadata_file, "r", encoding="utf-8") as f:
        lines = f.readlines()
        if len(lines) < 2:  # Header + at least one data row
            return study_ids, run_ids
        
        # Skip header, read data
        for line in lines[1:]:
            parts = line.strip().split("\t")
            if not parts:
                continue
            
            # Check if this is a study row (study_accession starts with SRP or is a number)
            # or a run row (run_accession starts with SRR)
            if parts[0].startswith("SRR"):
                run_ids.add(parts[0])
            elif parts[0] and (parts[0].startswith("SRP") or parts[0].isdigit()):
                study_ids.append(parts[0])
    
    return study_ids, run_ids


def query_study_runs(study_id: str, metadata_file: Optional[Path] = None) -> List[str]:
    """Query a single study to get its run IDs. Optionally update metadata table."""
    import urllib.request
    import urllib.parse
    import xml.etree.ElementTree as ET
    import time
    
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    run_ids = []
    
    try:
        # Search for the study to get its internal record ID
        study_search_params = {
            "db": "sra",
            "term": study_id,
            "retmode": "xml",
        }
        study_search_url = f"{base_url}esearch.fcgi?{urllib.parse.urlencode(study_search_params)}"
        
        time.sleep(0.34)  # Rate limiting
        
        with urllib.request.urlopen(study_search_url) as response:
            study_search_xml = response.read().decode()
            study_search_root = ET.fromstring(study_search_xml)
            study_record_ids = [id_elem.text for id_elem in study_search_root.findall(".//Id")]
            
            if not study_record_ids:
                logging.debug(f"Could not find record ID for study {study_id}")
                return []
        
        time.sleep(0.34)  # Rate limiting
        
        # Use efetch to get detailed record XML
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
            runs = fetch_root.findall(".//RUN")
            for run in runs:
                accession = run.get("accession")
                if accession and accession.startswith("SRR"):
                    run_ids.append(accession)
            
            # Update metadata table if provided
            if metadata_file and run_ids:
                # Check if runs are already in table
                existing_runs = set()
                if metadata_file.exists():
                    with open(metadata_file, "r", encoding="utf-8") as f:
                        for line in f.readlines()[1:]:  # Skip header
                            parts = line.strip().split("\t")
                            if len(parts) > 0 and parts[0].startswith("SRR"):
                                existing_runs.add(parts[0])
                
                # Append new runs to table
                new_runs = [r for r in run_ids if r not in existing_runs]
                if new_runs:
                    header_fields = ["run_accession", "study_accession", "bioproject_id"]
                    with open(metadata_file, "a", encoding="utf-8") as f:
                        for run_id in new_runs:
                            f.write(f"{run_id}\t{study_id}\t\n")  # bioproject_id empty for now
                    logging.debug(f"Added {len(new_runs)} runs from study {study_id} to metadata table")
        
        return run_ids
        
    except Exception as e:
        logging.debug(f"Failed to query study {study_id}: {e}")
        return []


def check_processed_runs(output_dir: Path) -> set[str]:
    """Check which runs have already been processed by looking for output directories.
    
    A run is considered processed if its output directory exists and contains results.
    """
    processed_runs = set()
    
    if not output_dir.exists():
        return processed_runs
    
    for run_dir in output_dir.iterdir():
        if not run_dir.is_dir():
            continue
        
        run_id = run_dir.name
        
        # Check if this directory contains processed results (VCF file or subdirectories with results)
        # Look for VCF files or sample directories
        vcf_files = list(run_dir.rglob("*.vcf.gz"))
        sample_dirs = [d for d in run_dir.iterdir() if d.is_dir() and d.name.startswith(run_id)]
        
        if vcf_files or sample_dirs:
            processed_runs.add(run_id)
    
    return processed_runs


def download_and_convert_bioproject(
    bioproject_id: str,
    output_dir: Path,
    threads: int = 1,
    skip_prefetch: bool = False,
) -> tuple[List[str], Path, set[str]]:
    """
    Query BioProject, generate/update metadata table with study IDs, 
    then return study IDs for sequential processing (runs queried on-demand).
    
    Supports resuming interrupted runs:
    - Checks for existing metadata table
    - Queries BioProject again to get fresh study IDs
    - Compares to find missing/new studies
    - Returns only unprocessed studies
    
    Args:
        bioproject_id: NCBI BioProject ID (e.g., "PRJNA123456")
        output_dir: Directory to save metadata table
        threads: Number of threads for fasterq-dump (not used in query phase)
        skip_prefetch: If True, skip prefetch step (not used in query phase)
    
    Returns:
        Tuple of (list of study IDs to process, path to metadata table file, set of already processed run IDs)
    """
    # Step 1: Create metadata file path
    metadata_file = output_dir / f"{bioproject_id}_metadata.tsv"
    
    # Step 2: Check existing metadata and processed runs
    existing_study_ids, existing_run_ids = read_metadata_table(metadata_file)
    processed_runs = check_processed_runs(output_dir)
    
    if metadata_file.exists():
        logging.info(f"Found existing metadata table: {metadata_file}")
        logging.info(f"Found {len(existing_study_ids)} study/studies and {len(existing_run_ids)} run(s) in existing metadata")
        logging.info(f"Found {len(processed_runs)} already processed run(s)")
    
    # Step 3: Query BioProject again to get fresh study IDs from XML response
    logging.info(f"Querying BioProject {bioproject_id} for fresh metadata...")
    _, _ = query_bioproject_with_metadata(bioproject_id, metadata_file=metadata_file)
    
    # Step 4: Read updated study IDs from metadata table (after fresh query)
    all_study_ids, all_run_ids = read_metadata_table(metadata_file)
    
    if not all_study_ids:
        logging.warning(f"No study IDs found in metadata table for BioProject {bioproject_id}")
        return [], metadata_file, processed_runs
    
    # Step 5: Determine which studies to process
    # We process all studies, but skip runs that are already processed
    # Studies without runs will be queried on-demand during processing
    new_studies = [s for s in all_study_ids if s not in existing_study_ids]
    if new_studies:
        logging.info(f"Found {len(new_studies)} new study/studies in BioProject")
    
    if existing_study_ids:
        logging.info(f"Found {len(existing_study_ids)} existing study/studies in metadata table")
    
    logging.info(f"Total: {len(all_study_ids)} study/studies, {len(all_run_ids)} run(s) in metadata")
    logging.info(f"Will process all studies, skipping {len(processed_runs)} already processed run(s)")
    
    return all_study_ids, metadata_file, processed_runs

