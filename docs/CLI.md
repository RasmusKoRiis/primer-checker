# Run Primer Checker on your own computer

The command-line tool checks consensus FASTA sequences against your primer
database and writes a CSV table and an optional, self-contained HTML report.
It uses the same Python analysis engine as the website. Once the software is
installed, analysis runs locally without uploading sequences or using Vercel.

This guide covers macOS (Intel and Apple Silicon), Linux, and Windows through
WSL. Commands below use a macOS/Linux terminal, including the Ubuntu terminal
inside WSL. The website's [English/Norsk CLI walkthrough](https://primercheck.rasmuskriis.no/#documentation-cli)
provides the setup steps in both languages.

## 1. Install Conda

If `conda --version` already works, keep your existing installation. Otherwise,
install [Miniforge](https://github.com/conda-forge/miniforge#requirements-and-installers)
for your operating system and processor, allow it to initialize your shell,
then close and reopen the terminal. Choose the macOS arm64 installer for an
Apple Silicon Mac, or x86_64 for an Intel Mac.

On Windows, first install [WSL with Ubuntu](https://learn.microsoft.com/en-us/windows/wsl/install).
Open Ubuntu and install the **Linux** version of Miniforge inside it. Run all
remaining commands there: Bioconda does not provide native Windows packages.

Check that Conda is available:

```bash
conda --version
```

The CLI needs Python and NCBI BLAST+. The environment below installs both;
you do not need Node.js, npm, FastAPI, or the web application's requirements files.

## 2. Download the complete source code

[Download the current source ZIP](https://github.com/RasmusKRiis/primer-checker/archive/refs/heads/feat/webapp-vercel.zip),
extract it, and open a terminal in the extracted folder. For example, if you
extracted it in Downloads on a Mac:

```bash
cd ~/Downloads/primer-checker-feat-webapp-vercel
```

Use the actual path if your browser chose a different folder name. On WSL,
prefer a folder inside your Linux home directory. Keep the whole folder:
`primer_checker.py` also needs `primer_analysis.py`, `primer_report.py`,
`report_text/`, and the database/example files.

If you already use Git, the alternative is:

```bash
git clone --branch feat/webapp-vercel --single-branch https://github.com/RasmusKRiis/primer-checker.git
cd primer-checker
```

The current website/CLI work is on `feat/webapp-vercel` in this fork while the
upstream pull request is under review. The links above select that branch;
the upstream default branch may have older behavior and files.

## 3. Create and activate the environment

The repository includes [`public/environment.yml`](../public/environment.yml).
You can also [download environment.yml directly](https://primercheck.rasmuskriis.no/environment.yml).
The file installs Python 3.12 and BLAST+ 2.16.0 in an environment named
`primer-checker`, using conda-forge and Bioconda without the defaults channel.

From the extracted repository folder, run:

```bash
CONDA_CHANNEL_PRIORITY=strict conda env create --file public/environment.yml
conda activate primer-checker
```

The first command downloads the dependencies, so it requires internet access
and may take several minutes. Strict channel priority applies to this command
without changing your global Conda settings. If using only the separately
downloaded file, replace `public/environment.yml` with its actual path, such
as `~/Downloads/environment.yml`. The environment file installs the runtime;
you still need the source folder from step 2.

Check the installed tools:

```bash
python --version
blastn -version
python primer_checker.py --help
```

Expect Python 3.12.x, BLAST 2.16.0+, and the Primer Checker options. The environment
sets `BLASTN_PATH=blastn` so the CLI selects Conda's BLAST, including if this
checkout contains a Linux deployment bundle in `bin/`.

## 4. Run the bundled test example

The bundled database and sequences are **invented software test data**, not a
validated primer set. Start with this small example to verify your installation:

```bash
mkdir -p results
python primer_checker.py \
  --primers primer_db/dummy_primers.json \
  --virus Demo-virus \
  --assay-type pcr \
  --fasta public/example.fasta \
  --output results/demo.csv \
  --html-report results/demo.html
```

Expected output: **2 sequence records × 3 primers = 6 comparisons**, all six
with hits. One comparison has a mismatch, `9:T>A`, for `DUMMY_F`; the reverse
primer hits have zero mismatches. A dummy-data notice in the terminal is expected.

Open `results/demo.html` in your browser by double-clicking it. On macOS, you
can also run `open results/demo.html`; on Linux, `xdg-open results/demo.html`.
In WSL, use `explorer.exe results` and open the HTML file from that folder.
The report works offline and includes English/Norsk controls, alignment views,
filters, and its own Help and documentation tab. The CSV contains the full
comparison table and can be opened in spreadsheet software.

## 5. Use your own sequences and database

Create a JSON database with the website's **Build a database** control and
download it, or supply an existing supported database. The [database guide](../primer_db/README.md)
describes the schema. Enter all oligos **5′ to 3′ as ordered**, including reverse
primers; do not reverse-complement them before entry.

Save the database as `my-primers.json` and the sequences as `my-sequences.fasta`
in your project folder, or substitute their actual paths. Quote paths that
contain spaces. Validate the database first:

```bash
python primer_checker.py --primers my-primers.json --validate-primers
```

Then run the analysis, replacing `YOUR_DATABASE_ORGANISM` with an organism label
present in your JSON:

```bash
python primer_checker.py \
  --primers my-primers.json \
  --virus "YOUR_DATABASE_ORGANISM" \
  --assay-type pcr \
  --fasta my-sequences.fasta \
  --output results/my-analysis.csv \
  --html-report results/my-analysis.html
```

Create the `results` directory first if you skipped the example. Existing files
at these output paths will be replaced; choose new names to retain earlier runs.

Useful options:

| Option                              | Meaning                                                                                               |
| ----------------------------------- | ----------------------------------------------------------------------------------------------------- |
| `--primers`                         | Path to your primer JSON database. Relative BED/FASTA assets are resolved from the database's folder. |
| `--virus`                           | Organism label from that database.                                                                    |
| `--assay-type pcr`, `ngs`, or `all` | Restrict to one technology, or include both. The default is `all`.                                    |
| `--assay-id`                        | Select an exact `scheme_id` or `panel_id` from the database.                                          |
| `--fasta file1.fasta file2.fasta`   | Analyze several FASTA files in one invocation.                                                        |
| `--output`                          | CSV destination.                                                                                      |
| `--html-report`                     | Optional interactive HTML destination.                                                                |
| `--help`                            | All available arguments.                                                                              |

FASTA records start with `>` followed by a unique ID, then sequence lines.
Only the first space-delimited token is the ID. Use nucleotide IUPAC symbols
such as A, C, G, T, N, R, and Y. FASTQ files are not FASTA files.
See the [FASTA guide](https://primercheck.rasmuskriis.no/#documentation-fasta)
for examples that also work with the website.

## 6. Influenza segments and subtypes

Use the type or exact subtype tag defined in **your database**. For example,
if it contains type A primers with an H5N1 tag:

```bash
python primer_checker.py \
  --primers my-primers.json \
  --virus Influenza \
  --flu-type H5N1 \
  --assay-type pcr \
  --fasta influenza.fasta \
  --output results/influenza.csv \
  --html-report results/influenza.html
```

This is an input template, not a claim that the dummy database contains that
assay. `--flu-type A` selects all type A primers in the supplied database;
an exact subtype selects its tagged primers plus untagged primers of that type.

Match segment labels in the JSON to FASTA IDs such as `>01-PB2|sample_001` or
`>06-NA|sample_001`. A segment token consists of one or two digits, a hyphen,
and the database segment label, before any whitespace. The number does not
determine the segment. Labels are case-insensitive, and supported segments
and subtypes come from the database.

## 7. Larger batches and results

The CLI does not enforce the website's 3 MB / 200-record limits or its
240-second request deadline. Runtime and memory still depend on your data and
primer count. BLAST has a shared maximum of 5,000 target sequences per primer
search, so split very large multi-record inputs; do not assume unlimited
reporting just because the web limits are absent.

For several files using the same organism/assay, pass them after `--fasta`.
There is also a folder wrapper that infers targets from filenames. Inspect
its routing before launching a mixed batch:

```bash
python scripts/run_primer_checker_batch.py \
  --input-folder my-fasta-folder \
  --primers my-primers.json \
  --dry-run \
  --fail-on-unclassified
```

Check every inferred target. Then remove `--dry-run` and add
`--output results/batch.csv --html-report results/batch.html`. Add `--recursive`
to include subfolders. The wrapper uses all matching assays; use the main CLI
if you need `--assay-type` or `--assay-id`. It skips unclassified files unless
`--fail-on-unclassified` is set. See the [README batch examples](../README.md#cli-usage).

The CSV and HTML report describe sequence compatibility, not measured assay
performance. Review no-hit records separately from zero-mismatch hits.
See the [analysis and report guide](https://primercheck.rasmuskriis.no/#documentation-method)
and the HTML report's built-in documentation for interpretation.

## 8. Return to the tool and record your versions

For each new terminal session:

```bash
conda activate primer-checker
cd /path/to/primer-checker
```

Use your actual source directory. After your work, `conda deactivate` leaves
the environment. To update a Git checkout without discarding local edits, use
`git pull --ff-only`; ZIP users can extract a fresh download into a new folder.

Save the code version, input FASTA files, primer database and any referenced
assets, exact command, and output reports. These commands record the runtime:

```bash
python --version
blastn -version
conda list --explicit > results/conda-explicit.txt
```

If you used Git, also save `git rev-parse HEAD`. An explicit Conda export pins
all resolved package builds **for the same operating system and architecture**;
the downloadable YAML is a portable specification, not a complete lockfile.
The hosted website currently bundles BLAST 2.15.0; this environment pins 2.16.0
to support native Apple Silicon as well as Intel/Linux. Search settings and
Python logic are shared, but BLAST-version differences can affect alignments.
Record the version when comparing results. The CLI's CSV/HTML do not include
the website's separate provenance JSON download.

## Troubleshooting

If Python or BLAST shows a different version after activation, check
`command -v python` and `command -v blastn`: both should be inside the active
environment. An existing virtual environment or custom shell PATH may take
precedence. Reopen a terminal, deactivate other environments, then activate
`primer-checker`. For an explicit selection, run `"$CONDA_PREFIX/bin/python"`
instead of `python` and set `export BLASTN_PATH="$CONDA_PREFIX/bin/blastn"`.

| Symptom                                                                     | What to do                                                                                                                                     |
| --------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| `conda: command not found`                                                  | Finish Miniforge installation and shell initialization, then reopen the terminal. On WSL, install it inside Ubuntu.                            |
| `conda activate` asks for initialization                                    | Run `conda init zsh` for macOS's default shell or `conda init bash` for Ubuntu, then reopen the terminal. Use the shell you actually run.      |
| `PackagesNotFoundError` on Windows                                          | Use WSL and Linux Conda. The Bioconda environment cannot be installed by native Windows Conda.                                                 |
| The environment already exists                                              | Run `conda activate primer-checker`. To create a separate copy, add `--name primer-checker-new` to the create command and activate that name.  |
| `can't open file primer_checker.py`                                         | Change into the extracted source folder; the environment does not install this script as a command.                                            |
| `BLAST executable unavailable`, wrong architecture, or shared-library error | Activate the environment and check `blastn -version`. Run `export BLASTN_PATH="$CONDA_PREFIX/bin/blastn"` to select its executable explicitly. |
| Missing `report_text/english.json` or a Python module                       | Download/extract the entire repository, not only `primer_checker.py`.                                                                          |
| No matching primers, unknown organism, or no results for a segment          | Validate the JSON, check organism/assay/subtype labels, and compare segment tags in the FASTA IDs with the database.                           |
| A path is not found                                                         | Check the working directory and quote paths containing spaces. Keep BED/FASTA assets beside the database as its relative paths specify.        |
| No BLAST hits                                                               | Inspect the selected primers and sequences. A no-hit result is not a zero-mismatch result. Check terminal messages for BLAST failures.         |
| The process is slow or runs out of memory                                   | Try a smaller FASTA batch or one assay at a time. Large HTML reports also need browser memory.                                                 |

Package installation follows the [Bioconda channel guidance](https://bioconda.github.io/)
and Conda's [environment-file workflow](https://docs.conda.io/projects/conda/en/stable/commands/env/create.html).
