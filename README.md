## Installation and Running

### 1) Environment setup

Create and activate a conda environment (Python 3.9 required):

```bash
conda create -n cowbat_hybrid python=3.9 -y
conda activate cowbat_hybrid
```

Or activate existing environment:

```bash
conda activate /mnt/nas2/virtual_environments/cowbat_hybrid
```

Install required system packages:

```bash
sudo apt-get update
sudo apt-get install -y build-essential pkg-config libcurl4-openssl-dev libssl-dev zlib1g-dev mummer ncbi-epcr
```

Clone repository and install Python requirements:

```bash
git clone https://github.com/OLC-LOC-Bioinformatics/CowbatHybrid.git
cd CowbatHybrid
SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True pip install -r requirements.txt
pip install -e .
```

Install bioinformatics dependencies (conda-forge + bioconda channels required):

```bash
mamba install -y flye racon bbmap blast ectyper mob_suite clark nanoplot prodigal filtlong porechop sistr_cmd mash pilon \
  unicycler=0.5.1 dnaapler=1.4.0 mmseqs2=13.45111
```

> Note: with Unicycler `0.5.1`, do **not** use `--no_correct` (unsupported in this version).

Verify pipeline dependencies:

```bash
python -c "from cowbathybrid.dependency_checks import check_dependencies; print(check_dependencies())"
```

Expected output:

```text
True
```

---

### 2) Input format

Provide a CSV with headers:

- `MinION`
- `Illumina_R1`
- `Illumina_R2`
- `OutName`

Use absolute paths for FASTQ files.

---

### 3) Run the pipeline

Check options:

```bash
cowbat-hybrid-assembly.py --help
```

Example command:

```bash
cowbat-hybrid-assembly.py \
  -i /data/pipelines/seqs/test.csv \
  -r /mnt/nas2/databases/assemblydatabases/0.5.0.18 \
  -t 46 \
  -o /data/pipelines/test_MIN0157 \
  --min-read-length 6000 \
  --filtlong-keep-percent 95 \
  --run-dnaapler
```

---

### 4) Tuning notes (long-read filtering)

- `--min-read-length 6000` is a good starting point for high-coverage datasets.
- `--filtlong-keep-percent 95` keeps most reads while applying quality filtering.
- Tune per sample:
  - Increase min length for stricter filtering.
  - Lower min length for lower coverage samples.
  - Adjust keep-percent if quality is poor or assembly is fragmented.

---

### 5) Output

Main outputs are written to your chosen output directory, including:

- `BestAssemblies/` (final assembly FASTA files)
- `reports/` (typing and summary reports; e.g., `combinedMetadata.csv`)

### Additional runtime notes

1. **Pilon memory limit**
   By default, some Bioconda `pilon` wrappers set max JVM memory to `-Xmx1g`, which may be too low.

   - Find pilon wrapper:
     ```bash
     which pilon
     ```
   - Edit the wrapper and change:
     ```python
     default_jvm_mem_opts = ['-Xms512m', '-Xmx1g']
     ```
     to e.g.:
     ```python
     default_jvm_mem_opts = ['-Xms512m', '-Xmx12g']
     ```

2. **`famap` / `fahash` not found**
   If dependency checks fail for `famap`/`fahash`, locate them and add their directory to `PATH`.

   Example:
   ```bash
   find "$CONDA_PREFIX" -type f \( -name famap -o -name fahash \) 2>/dev/null
   ```

   If found, add parent directory to `PATH`, e.g.:
   ```bash
   export PATH="/path/to/dir/with/famap_and_fahash:$PATH"
   ```

   > Older instructions referenced a Python 3.5 site-packages path.  
   > For this pipeline use Python 3.9, so always use paths from your active environment (`$CONDA_PREFIX`).
