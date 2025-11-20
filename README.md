# QuickVar

QuickVar provides a cross-platform, one-command workflow to align *Candida glabrata* sequencing reads and perform variant calling. The pipeline automatically downloads the reference genome, installs the required bioinformatics tools, and produces alignment and variant files with minimal user input.

## Features
- Automatic installation of alignment and variant-calling tools via Micromamba.
- Automatic download and indexing of the *Candida glabrata* CBS138 reference genome.
- Supports *Candida glabrata* (default) and *Candida auris* references (switch via `--reference`).
- Supports single-end and paired-end FASTQ files (optionally compressed with gzip).
- **NEW:** Download and process SRA files directly from NCBI BioProjects using `--bioproject`.
- Generates sorted BAM files with indexes and VCF outputs per sample.
- Cross-platform support (macOS, Linux, Windows) with identical commands.
- Variant calling defaults to haploid (`--ploidy 1`) but can be configured per run.
- Amplicon mode (`--amplicon`) emits per-position mutation frequency summaries.
- Optional PCR duplicate removal via `--deduplicate` (powered by `samtools markdup`).
- **NEW:** VCF annotation with gene information from GFF files using `--annotate`.

## Prerequisites
- Python 3.10 or newer.
- Internet connection for downloading the reference genome and Micromamba binaries.

Micromamba is installed automatically into a user-level cache (`~/.quickvar`). No system-wide changes are required.

## Quick Start

### Optional: Install the CLI
```bash
python -m pip install -e .
```
Installing in editable mode exposes the `quickvar`, `quickvar-install`, and `quickvar-align` commands on your PATH. You can also run the modules directly (as shown below) without installing.

**Note:** All dependencies are automatically installed in the Micromamba environment when you run `python -m quickvar.install`. BioProject/SRA functionality uses Python's standard library only (no extra dependencies needed).

### 1. Install QuickVar Dependencies

#### Linux / macOS
```bash
python -m quickvar.install
```
This downloads Micromamba (if needed) and creates the `quickvar` environment containing `minimap2`, `samtools`, `bcftools`, and `sra-tools`. BioProject queries use NCBI Entrez API (no additional dependencies needed).

#### Windows (WSL2)
1. Open **PowerShell as Administrator** and enable Ubuntu on WSL2:
   ```powershell
   wsl --install -d Ubuntu
   ```
   Reboot if prompted and complete the first-run setup (username/password).
2. Launch the **Ubuntu** terminal and run QuickVar from there:
   ```bash
   cd /mnt/c/path/to/QuickVar   # adjust the repo path
   python -m quickvar.install
   ```
   Windows drives are available inside WSL under `/mnt/<drive-letter>/...`, and you can open WSL paths from Explorer via `\\wsl$\Ubuntu\`.

### 2. Run the Pipeline

#### Option A: Using Local FASTQ Files
```bash
python -m quickvar.align --input /path/to/fastqs --output /path/to/results
```
- `--input` can point to a single FASTQ file or a directory containing one or more FASTQ files.
- When not provided, `--output` defaults to a `Results` directory in the current working directory.

#### Option B: Download from NCBI BioProject
```bash
python -m quickvar.align --bioproject PRJNA123456 --output /path/to/results
```
- Downloads all SRA runs from the specified BioProject, converts them to FASTQ, and processes them.
- Use `--skip-prefetch` to skip the prefetch step (faster, but less caching).

Each sample results in:
- `sample.sorted.bam` and `sample.sorted.bam.bai`
- `sample.vcf.gz` and `sample.vcf.gz.tbi`
- (Optional) Annotated VCF with gene information if `--annotate` is used

### Sample Detection
- Paired-end files should include `_R1`/`_R2`, `.R1`/`.R2`, or `_1`/`_2` in their names.
- Single-end files that do not follow pairing conventions are processed individually.

### Cleaning Up
To remove the QuickVar Micromamba environment:
```bash
python -m quickvar.install --remove
```

## Hands-On Tutorial

The example below walks through a complete run using public test data.

1. **Prepare working directory**  
   ```bash
   mkdir -p ~/quickvar-demo && cd ~/quickvar-demo
   git clone https://github.com/<your-org>/quickvar.git
   cd quickvar
   ```

2. **Install dependencies**  
   ```bash
   python -m quickvar.install
   ```
   This downloads Micromamba (if necessary) and builds the `quickvar` environment with `minimap2`, `samtools`, `bcftools`, and `sra-tools`. BioProject queries use NCBI Entrez API via Python's standard library (no additional dependencies needed).  
   *Windows users:* run this command inside an Ubuntu WSL2 session (see Quick Start step 1).

3. **Run alignment and variant calling**  
   ```bash
   python -m quickvar.align \
     --input test_data/amplicon/glabrata_amplicon.fastq.gz \
     --output DemoResults \
     --amplicon \
     --annotate
   ```
   Progress prints to the terminal. Results for each sample land inside `DemoResults/<sample>/`. 
   - Add `--deduplicate` if you want duplicate reads removed before variant calling.
   - The `--amplicon` flag adds a per-position summary TSV in addition to the alignment/variant files.
   - The `--annotate` flag adds gene annotations to the VCF file (requires GFF file for the reference).

4. **(Optional) Use your own FASTQs or download from NCBI**  
   - **Local FASTQs:** Point `--input` at your FASTQ file or a directory containing multiple FASTQs. Paired-end files are paired automatically when they follow `_R1`/`_R2` (or similar) naming.
   - **NCBI BioProject:** Use `--bioproject PRJNA123456` to automatically download and process all SRA runs from a BioProject. The pipeline will download SRA files, convert them to FASTQ, and process them automatically.

5. **Inspect outputs**  
   ```bash
   ls DemoResults/sample/
   samtools flagstat DemoResults/sample/sample.sorted.bam
   bcftools view DemoResults/sample/sample.vcf.gz | head
   column -t DemoResults/sample/sample_amplicon.tsv | head
   ```
   If `--annotate` was used, check the VCF INFO column for gene annotations:
   ```bash
   bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/GENE_ID\t%INFO/GENE_NAME\n' DemoResults/sample/sample.vcf.gz | head
   ```

6. **Clean up (optional)**  
   ```bash
   python -m quickvar.install --remove
   rm -rf ~/.quickvar
   ```
   This removes the Micromamba environment and cached reference data.

### Built-in Test Dataset

The repository ships with a tiny synthetic amplicon dataset in `test_data/amplicon/glabrata_amplicon.fastq.gz`.

Run it end-to-end with:
```bash
python -m quickvar.align \
  --input test_data/amplicon/glabrata_amplicon.fastq.gz \
  --output AmpliconResults \
  --ploidy 1 \
  --amplicon \
  --threads 1
```

Append `--deduplicate` if you want the example to remove duplicate reads before summarising/variant calling.

Expected variant coordinates (contig + 1-based position) are listed in `test_data/amplicon/variants.tsv`. The dataset includes two SNPs plus an example insertion (`+AT`) and deletion (`-T`), so with the haploid default you should see homozygous alternate calls at those loci. When `--amplicon` is provided, the pipeline also writes `<sample>_amplicon.tsv` containing per-position alternate counts and frequencies (including insertions/deletions) using all reads at each site (no pileup depth cap), plus `igv_depth` (raw IGV-style coverage), `estimated_coverage` (mean flanking coverage within ±5 bp excluding the focal base), and `estimated_frequency` (alternate counts over the estimated coverage). A secondary `<sample>_amplicon_indels.tsv` highlights 10 bp windows around indels with wild-type read counts (`wt_count_in_10bp`).

### Using the *Candida auris* reference

Supply `--reference c_auris` when installing/aligning to use the bundled *Candida auris* B11221 reference:

```bash
python -m quickvar.align \
  --input your_reads.fastq.gz \
  --output AurisResults \
  --reference c_auris \
  --amplicon \
  --annotate
```

### Downloading and Processing SRA Files from NCBI BioProjects

QuickVar can automatically download and process SRA files from NCBI BioProjects:

```bash
python -m quickvar.align \
  --bioproject PRJNA123456 \
  --output BioProjectResults \
  --reference c_glabrata \
  --amplicon \
  --annotate
```

This will:
1. Query the BioProject to find all associated SRA runs
2. Download SRA files (cached in `~/.quickvar/sra/`)
3. Convert SRA files to FASTQ format
4. Process each sample through the alignment and variant calling pipeline

**Options:**
- `--skip-prefetch`: Skip the prefetch step (faster, but fasterq-dump will download if needed)
- SRA files are cached, so re-running with the same BioProject will reuse cached files

**Note:** Some SRA runs may require dbGaP authorization (controlled access). The pipeline will skip these with a warning and continue processing other runs.

### VCF Annotation with Gene Information

Add gene annotations to your VCF files using the bundled GFF files:

```bash
python -m quickvar.align \
  --input your_reads.fastq.gz \
  --output AnnotatedResults \
  --annotate
```

The `--annotate` flag adds the following INFO fields to variants:
- `GENE_ID`: Gene identifier from the GFF file
- `GENE_NAME`: Gene name (if available)
- `FEATURE_TYPE`: Type of feature (gene/CDS/mRNA)
- `PRODUCT`: Gene product/description (if available)

Annotation files are built once and cached in `~/.quickvar/reference/` for reuse. If annotation fails (e.g., no GFF file available), the pipeline continues with a warning and produces an unannotated VCF.

## Development
- `pyproject.toml` configures QuickVar as a Python package with console entry points.
- Unit tests (coming soon) can be run with `pytest` within the QuickVar environment.

## License
QuickVar is distributed under the MIT License. See `LICENSE` for details.
