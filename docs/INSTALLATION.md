# Installation

MotifSTaR is a standalone Python program. It can run on a personal computer or a Linux high-performance-computing (HPC) system when Python and the required packages are available. The current release was developed and validated on Linux with Python 3.10/3.11.

## Recommended method: Miniconda

Miniconda provides Python, Conda, and isolated environments without requiring administrator access. Always obtain the installer from the [official Miniconda page](https://www.anaconda.com/docs/getting-started/miniconda/install) and select the file for your operating system and processor.

### Linux

This example is for Linux x86-64. Use the official page for another architecture.

```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Read the licence shown by the installer, choose an installation directory you own, and allow the installer to initialize Conda. Close and reopen the shell, then check:

```bash
conda --version
```

If Conda was installed but is not yet active in the current shell:

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
```

Replace that path if you chose another installation directory.

### Windows

1. Download the Windows Miniconda installer from the official page.
2. Run the installer for the current user.
3. Open **Miniconda Prompt** from the Start menu.
4. Change to the cloned MotifSTaR directory and run the environment commands below.

### macOS

Download the installer that matches Apple silicon (`arm64`) or Intel (`x86_64`). Open Terminal after installation and confirm that `conda --version` works.

## Create the environment

In the MotifSTaR repository directory:

```bash
conda env create --file environment.yml
conda activate motifstar
python motifstar.py --version
```

Expected result:

```text
motifstar.py 1.4.7
```

The environment includes:

- Python 3.10;
- Matplotlib for figures;
- openpyxl for XLSX input and assisted-review workbooks;
- Pillow and ReportLab for graphics support and fallback rendering.

Update the environment after pulling a newer release:

```bash
conda env update --name motifstar --file environment.yml --prune
```

Remove it if it is no longer needed:

```bash
conda deactivate
conda env remove --name motifstar
```

Conda's environment-management documentation is available at [docs.conda.io](https://docs.conda.io/projects/conda/en/stable/user-guide/tasks/manage-environments.html).

## Alternative: Python virtual environment and pip

If Python 3.10 or later is already installed:

### Linux or macOS

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python motifstar.py --version
```

### Windows PowerShell

```powershell
py -3.10 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python motifstar.py --version
```

## HPC installation

Follow the software and storage policies of the HPC facility. A user-owned Conda installation can be created in a permitted home or project directory. On systems that provide Conda as a module, load that module instead of installing a second copy.

For a scheduler job, initialize Conda before activation. A typical PBS job contains:

```bash
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate motifstar
cd /path/to/MotifSTaR
python -u motifstar.py --help
```

Set `--workers` to the number of CPU cores allocated to the job. MotifSTaR also reads `PBS_NCPUS` and `SLURM_CPUS_PER_TASK` when `--workers` is omitted; otherwise the automatic default is at most four workers.

## Verify the installation

Run the version, syntax, and regression checks:

```bash
python motifstar.py --version
python -m py_compile motifstar.py
python -m unittest discover -s tests -p "test_*.py" -v
```

Then perform an input audit before a production run:

```bash
python -u motifstar.py \
  --svg-root examples/input/svgs \
  --catalog config/variant_catalog_without_offtargets.GRCh38.json \
  --annotation config/gene_strands_companion_locus.xlsx \
  --output examples/output/input_check \
  --component-mode primary \
  --validate-only
```

## Common installation problems

### `conda: command not found`

Open a new terminal or source the Conda initialization script from the chosen installation directory.

### `ModuleNotFoundError`

Activate the environment and install/update it again:

```bash
conda activate motifstar
conda env update --name motifstar --file environment.yml --prune
```

### XLSX file is not created

The main CSV outputs can be written without openpyxl, but the assisted-review workbook requires it. Confirm:

```bash
python -c "import openpyxl; print(openpyxl.__version__)"
```

### Figures are not created

Confirm that Matplotlib is available:

```bash
python -c "import matplotlib; print(matplotlib.__version__)"
```

For a headless Linux system, MotifSTaR uses a non-interactive plotting backend; no desktop display is required.
