"use client";

import { useLanguage } from "./language";
import { Download } from "./pixel-icons";

const source = "https://github.com/RasmusKRiis/primer-checker";
const branch = "feat/webapp-vercel";

export default function CliSetup() {
  const { t } = useLanguage();
  const code = (label: string, command: string) => (
    <pre tabIndex={0} aria-label={t(label)}>
      <code>{command}</code>
    </pre>
  );
  return (
    <section
      className="card documentation-section format-guide cli-setup"
      id="documentation-cli"
    >
      <h2>{t("Run locally with the CLI")}</h2>
      <p>
        {t(
          "Run the same Python analysis engine on your own computer and create CSV and HTML reports. After installation, analysis works offline and your sequences stay on your computer. The website's upload and batch limits do not apply; available memory and processing power still matter.",
        )}
      </p>
      <div className="format-options">
        <a className="button secondary" href="/environment.yml" download>
          <Download size={16} />
          {t("Download Conda environment")}
        </a>
        <a
          className="button secondary"
          href={`${source}/archive/refs/heads/${branch}.zip`}
        >
          <Download size={16} />
          {t("Download source code (ZIP)")}
        </a>
        <a
          className="button secondary"
          href={`${source}/blob/${branch}/docs/CLI.md`}
          target="_blank"
          rel="noreferrer"
        >
          {t("Full CLI tutorial (English)")}
        </a>
      </div>
      <h3>{t("1. Install Conda")}</h3>
      <p>
        {t(
          "If conda --version already works, use your existing installation. Otherwise, install Miniforge for your computer, allow shell initialization, then reopen the terminal. Choose arm64 for Apple Silicon or x86_64 for Intel Macs.",
        )}
      </p>
      <p>
        <a href="https://github.com/conda-forge/miniforge#requirements-and-installers">
          {t("Miniforge installers")}
        </a>
        {" · "}
        <a href="https://learn.microsoft.com/en-us/windows/wsl/install">
          {t("Windows: install WSL / Ubuntu first")}
        </a>
      </p>
      <p>
        {t(
          "On Windows, install Linux Miniforge inside Ubuntu / WSL and run the commands there. This Bioconda environment supports macOS and Linux, not native Windows.",
        )}
      </p>
      <h3>{t("2. Download and extract the source code")}</h3>
      <p>
        {t(
          "Download the source ZIP above, extract it, and open a terminal in that folder. Keep the complete folder, including report_text. The download selects the current feature branch while the upstream pull request is under review.",
        )}
      </p>
      {code(
        "Open the source folder (example)",
        "cd ~/Downloads/primer-checker-feat-webapp-vercel",
      )}
      <p>{t("Adjust this path to wherever you extracted the ZIP.")}</p>
      <h3>{t("3. Create and activate the environment")}</h3>
      <p>
        {t(
          "The source ZIP includes public/environment.yml. It installs Python 3.12 and BLAST 2.16.0 into an environment named primer-checker. No website packages are needed. Creating the environment requires internet access.",
        )}
      </p>
      {code(
        "Create the Conda environment",
        "CONDA_CHANNEL_PRIORITY=strict conda env create --file public/environment.yml\nconda activate primer-checker\npython --version\nblastn -version\npython primer_checker.py --help",
      )}
      <p>
        {t(
          "If you downloaded environment.yml separately, use its actual path after --file. That file installs the runtime; you still need the source code. The environment selects its own BLAST executable.",
        )}
      </p>
      <h3>{t("4. Check the installation with dummy data")}</h3>
      {code(
        "Run the CLI dummy example",
        "mkdir -p results\npython primer_checker.py \\\n  --primers primer_db/dummy_primers.json \\\n  --virus Demo-virus \\\n  --assay-type pcr \\\n  --fasta public/example.fasta \\\n  --output results/demo.csv \\\n  --html-report results/demo.html",
      )}
      <p>
        {t(
          "Expect 2 records, 3 primers, 6 hits, and one mismatch (9:T>A for DUMMY_F). Open results/demo.html in your browser. These are invented test sequences, not a validated primer set.",
        )}
      </p>
      <h3>{t("5. Analyze your own inputs")}</h3>
      <p>
        {t(
          "Build and download a JSON database in the website, or use your own supported database. Save it as my-primers.json and your FASTA as my-sequences.fasta, or substitute their paths below. Replace YOUR_DATABASE_ORGANISM with the organism label in your database. Keep all primers in their original 5′ to 3′ orientation, including reverse primers.",
        )}
      </p>
      {code(
        "Validate and analyze your own inputs",
        'python primer_checker.py --primers my-primers.json --validate-primers\n\npython primer_checker.py \\\n  --primers my-primers.json \\\n  --virus "YOUR_DATABASE_ORGANISM" \\\n  --assay-type pcr \\\n  --fasta my-sequences.fasta \\\n  --output results/my-analysis.csv \\\n  --html-report results/my-analysis.html',
      )}
      <p>
        {t(
          "Create the results folder if you skipped the example. Use new output filenames to keep previous reports. Use --assay-type ngs for sequencing panels, --assay-id for one scheme or panel, or list multiple files after --fasta. For influenza, add --flu-type with a database type/subtype and use matching segment tags in FASTA headers.",
        )}
      </p>
      <details>
        <summary>{t("Returning to the CLI and recording versions")}</summary>
        <p>
          {t(
            "In each new terminal, activate primer-checker and change into the source folder. Save your command, inputs, database, code version, and reports. The tutorial explains environment exports, batch folders, and influenza examples.",
          )}
        </p>
        {code(
          "Activate the environment again",
          "conda activate primer-checker\ncd /path/to/primer-checker\nconda list --explicit > results/conda-explicit.txt",
        )}
        <p>
          {t(
            "The website uses BLAST 2.15.0; this portable environment uses 2.16.0 for native Apple Silicon support. Version differences can affect alignments. Keep a record of the BLAST version when comparing runs. An explicit Conda export applies to the same operating system and architecture.",
          )}
        </p>
      </details>
      <details>
        <summary>{t("CLI setup troubleshooting")}</summary>
        <ul>
          <li>
            {t(
              "Conda not found: finish Miniforge installation and shell initialization, then reopen the terminal. On Windows, work inside Ubuntu / WSL.",
            )}
          </li>
          <li>
            {t(
              "Environment already exists: activate primer-checker, or create another copy with --name primer-checker-new.",
            )}
          </li>
          <li>
            {t(
              "Script or report text not found: extract the entire repository and run commands from the folder containing primer_checker.py.",
            )}
          </li>
          <li>
            {t(
              "BLAST not found or wrong architecture: activate the environment and check blastn -version. The full tutorial explains how to select the environment's executable explicitly.",
            )}
          </li>
          <li>
            {t(
              "No matching primers: validate the database and check organism, assay, subtype, and FASTA segment labels. For slow runs, split inputs into smaller batches; BLAST reports at most 5,000 target sequences per primer search.",
            )}
          </li>
        </ul>
      </details>
    </section>
  );
}
