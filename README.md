# QuickVar

QuickVar provides a cross-platform, one-command workflow to align *Candida glabrata* sequencing reads and perform variant calling. The pipeline automatically downloads the reference genome, installs the required bioinformatics tools, and produces alignment and variant files with minimal user input.

## Features
- Automatic installation of alignment and variant-calling tools via Micromamba.
- Automatic download and indexing of the *Candida glabrata* CBS138 reference genome.
- Supports *Candida glabrata* (default) and *Candida auris* references (switch via `--reference`).
- Supports single-end and paired-end FASTQ files (optionally compressed with gzip).
- Generates sorted BAM files with indexes and VCF outputs per sample.
- Cross-platform support (macOS, Linux, Windows) with identical commands.
- Variant calling defaults to haploid (`--ploidy 1`) but can be configured per run.
- Amplicon mode (`--amplicon`) emits per-position mutation frequency summaries.
- Optional PCR duplicate removal via `--deduplicate` (powered by `samtools markdup`).

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

### 1. Install QuickVar Dependencies

#### Linux / macOS
```bash
python -m quickvar.install
```
This downloads Micromamba (if needed) and creates the `quickvar` environment containing `minimap2`, `samtools`, and `bcftools`.

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
```bash
python -m quickvar.align --input /path/to/fastqs --output /path/to/results
```
- `--input` can point to a single FASTQ file or a directory containing one or more FASTQ files.
- When not provided, `--output` defaults to a `Results` directory in the current working directory.

Each sample results in:
- `sample.sorted.bam` and `sample.sorted.bam.bai`
- `sample.vcf.gz` and `sample.vcf.gz.tbi`

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
   This downloads Micromamba (if necessary) and builds the `quickvar` environment with `minimap2`, `samtools`, and `bcftools`.  
   *Windows users:* run this command inside an Ubuntu WSL2 session (see Quick Start step 1).

3. **Run alignment and variant calling**  
   ```bash
   python -m quickvar.align \
     --input test_data/amplicon/glabrata_amplicon.fastq.gz \
     --output DemoResults \
     --amplicon
   ```
   Progress prints to the terminal. Results for each sample land inside `DemoResults/<sample>/`. Add `--deduplicate` if you want duplicate reads removed before variant calling. The `--amplicon` flag adds a per-position summary TSV in addition to the alignment/variant files.

4. **(Optional) Use your own FASTQs**  
   Point `--input` at your FASTQ file or a directory containing multiple FASTQs. Paired-end files are paired automatically when they follow `_R1`/`_R2` (or similar) naming.

5. **Inspect outputs**  
   ```bash
   ls DemoResults/sample/
   samtools flagstat DemoResults/sample/sample.sorted.bam
   bcftools view DemoResults/sample/sample.vcf.gz | head
   column -t DemoResults/sample/sample_amplicon.tsv | head
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
  --amplicon
```

## Development
- `pyproject.toml` configures QuickVar as a Python package with console entry points.
- Unit tests (coming soon) can be run with `pytest` within the QuickVar environment.

## License
QuickVar is distributed under the MIT License. See `LICENSE` for details.
