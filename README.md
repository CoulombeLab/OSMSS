# One Stop MS Shop (OSMSS)

**Statistical analysis and interpretation of quantitative proteomics data for non-coders.**

OSMSS is a no-code R Shiny application for differential abundance analysis of mass spectrometry proteomics data. It reads search results exported from Proteome Discoverer, MaxQuant, FragPipe, DIA-NN or Spectronaut. A short wizard records your experimental design and selects the analysis route. OSMSS checks whether the data are already normalized, then runs established statistical methods: proDA, limma, SAINTexpress and Fisher's exact test. Every analysis produces a **case-specific methods text** that records the counts at each filter, the checks, the decisions made and why, and the software versions used.

<!-- TODO before release: add the repository URL, release version, Zenodo DOI and the bioRxiv/JPR citation (see "Citing OSMSS"). -->

---

## Contents

- [What OSMSS does](#what-osmss-does)
- [System requirements](#system-requirements)
- [Installation](#installation)
- [Running OSMSS](#running-osmss)
- [Supported inputs](#supported-inputs)
- [Spectronaut export schemes (.rs files)](#spectronaut-export-schemes-rs-files)
- [Reference data you download yourself](#reference-data-you-download-yourself)
- [SAINTexpress (AP-MS)](#saintexpress-ap-ms)
- [Saved Experiments](#saved-experiments)
- [Data sent to external services](#data-sent-to-external-services)
- [Running OSMSS on a server](#running-osmss-on-a-server)
- [Troubleshooting](#troubleshooting)
- [Developer](#developer)
- [Citing OSMSS](#citing-osmss)
- [License](#license)

---

## What OSMSS does

| Experiment | Statistical route |
|---|---|
| Label-free DDA, group comparison | proDA on protein intensities (global-median or invariant-protein-panel normalization) |
| Affinity purification (AP-MS / IP) | SAINTexpress (intensity mode) against negative controls; proDA with the invariant-protein panel for fraction or cross-condition controls |
| Label-free DIA, group comparison | limma (eBayes, intensity trend, robust) with kNN imputation of sparse missing values, plus a presence/absence (Fisher) test |
| TMT, single plex or multi-plex | limma; IRS bridging or plex covariate for multi-plex designs |
| Phosphoproteomics, DIA or TMT | limma on site-level data; absolute, protein-adjusted and parent-protein views |
| No replicates (N=1) | Rank fold change and trend analysis |

**Interpretation tools:**
- volcano plots
- heatmaps
- PCA
- trend clusters (fuzzy c-means)
- protein coverage
- cross-experiment comparison
- GSEA (MSigDB)
- Kinase-Site Analysis (PhosR kinase–substrate scoring)
- residue-level prediction (The Kinase Library, Scansite, 14-3-3-Pred, AlphaFold structures, conservation)
- a local FragPipe PTM database

---

## System requirements

- **R ≥ 4.4.** R 4.6.1 was used for the manuscript, and R ≥ 4.4 is needed for safe loading of shared files (see [Saved Experiments](#saved-experiments)).
- **Windows, macOS or Linux.** SAINTexpress, used for affinity-purification (AP-MS) scoring, is installed separately (see [SAINTexpress](#saintexpress-ap-ms)).
- **Memory:** more than 8 GB of RAM is recommended. OSMSS starts three background R workers for long tasks, and each uses 200–400 MB.
- **Internet access** to install packages and to use the online annotation tools (UniProt, AlphaFold DB, Ensembl, Scansite, 14-3-3-Pred, MSigDB). Statistical analyses run offline.
- **A modern web browser.** OSMSS opens in your default browser.

**Before installing on Linux:** install the system libraries several packages need. On Ubuntu/Debian:

```bash
sudo apt install build-essential libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev
```

**Before installing on macOS:** if R asks to compile packages from source, install the Xcode command-line tools first:

```bash
xcode-select --install
```

---

## Installation

### Windows: the setup wizard

1. Install [R for Windows](https://cran.r-project.org/bin/windows/base/) (4.4 or later) if you do not have it.
2. Download OSMSS (on GitHub: **Code → Download ZIP**, or a release ZIP) and unzip it.
3. Double-click **`Install OSMSS (Windows).cmd`**. If Windows warns that the file came from the internet, choose **More info → Run anyway**.

The wizard asks where to install OSMSS and checks for R. It then offers two optional steps, each of which you can skip and do later:
- **SAINTexpress, for AP-MS:** point it to the downloaded `SAINTexpress_v3.6.3__2018-03-09.tar.gz`.
- **The Kinase Library files, for kinase prediction:** point it to the downloaded `.tsv` files or ZIP.

It installs the R packages OSMSS needs, which takes 10–30 minutes and needs internet. Packages you already have are **not updated**, and the last page lists them and any that failed. Finally it creates an **OSMSS** shortcut on the desktop. The shortcut opens OSMSS in your browser; closing its console window stops OSMSS. The package-install log is kept in `logs/` in the install folder.

### macOS and Linux: the setup script

1. Install R 4.4 or later from CRAN ([macOS](https://cran.r-project.org/bin/macosx/), [Linux](https://cran.r-project.org/bin/linux/)).
2. Download OSMSS and unzip it.
3. Start Setup:
   - **macOS:** double-click **`Install OSMSS (macOS).command`**. It opens in Terminal. If macOS blocks it, right-click it, choose **Open**, then **Open** again.
   - **Linux:** in a terminal, run `bash "Install OSMSS (Linux).sh"` from the OSMSS folder.

Setup runs the same steps as on Windows, in the terminal, using your system's file picker where one is available.

The differences:
- **Linux build tools:** on Linux, Setup checks the compilers and libraries R needs to build packages from source, and offers to install them with `apt` or `dnf` (this needs your password). Package installation takes 20–60 minutes on Linux.
- **SAINTexpress:**
  - On Linux, Setup uses the ready-made Linux program and compiles from source (`make`) only if that program does not run.
  - On macOS it always compiles. That takes about 10 minutes and needs Apple's Command Line Tools, which Setup offers to install.
- **Shortcut:**
  - macOS: an **OSMSS** app on the Desktop, which opens OSMSS in a Terminal window.
  - Linux: an **OSMSS** entry in the applications menu and on the desktop. Some desktops ask you to choose **Allow Launching** the first time.
- **Starting by hand:** you can also run `start-osmss.sh` in the install folder.

On any system, to install by hand, follow steps 1–3 below.

### 1. Get the code

Clone the repository, or download and unzip a release:

```bash
git clone <REPOSITORY_URL> OSMSS
```

### 2. Install all R packages (one command)

Paste this into an R console. It installs every CRAN and Bioconductor package OSMSS uses; the first install can take 15–30 minutes.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org"); BiocManager::install(c("shiny", "shinyjs", "bslib", "DT", "plotly", "sortable", "htmltools", "jsonlite", "xml2", "httr", "base64enc", "dplyr", "tidyr", "readr", "stringr", "glue", "scales", "data.table", "future", "promises", "ggplot2", "ggvenn", "pheatmap", "e1071", "MASS", "statmod", "ruv", "bio3d", "NGLVieweR", "msigdbr", "rstudioapi", "limma", "proDA", "QFeatures", "MsCoreUtils", "SummarizedExperiment", "S4Vectors", "impute", "PhosR", "Biostrings", "msa", "drawProteins", "clusterProfiler"), update = FALSE, ask = FALSE)
```

- **What it does:**
  - `BiocManager` picks the Bioconductor release that matches your R version and installs CRAN and Bioconductor packages together.
  - `msigdbr` downloads the MSigDB gene sets for GSEA the first time you use them.
- **Missing packages:** if an optional package fails to install, OSMSS still starts. The module that needs it shows a message naming the missing package.
- **To reproduce the manuscript's package versions exactly,** restore the lockfile instead. `renv.lock` records the release environment: R 4.6.1, Bioconductor 3.23 and 246 packages. From an R session in the OSMSS folder:

  ```r
  install.packages("renv")
  renv::activate()   # once: gives this folder its own package library
  renv::restore()    # installs the exact versions listed in renv.lock
  ```

  OSMSS then uses those versions whenever R is started in this folder.

### 3. Optional: install SAINTexpress

This is needed only for AP-MS scoring with SAINTexpress; see [SAINTexpress (AP-MS)](#saintexpress-ap-ms).

---

## Running OSMSS

**From RStudio:** open `app.R` and click **Run App**.

**From an R console** (working directory set to the OSMSS folder):

```r
shiny::runApp(launch.browser = TRUE)
```

**From a terminal:**

```bash
Rscript app.R
```

OSMSS opens in your web browser. The **Home** page lists the analysis sections. Each section starts with an upload wizard that asks about your experiment, then shows the route it selected and why.

---

## Supported inputs

Export the files below from your search software. The wizard detects the file type and tells you if something is missing. The upload limit is **1 GB per file**.

| Experiment | Software | Files to upload |
|---|---|---|
| Label-free DDA, and AP-MS | Proteome Discoverer | Protein or peptide export (CSV) with raw, non-normalized `Abundance: F#:` columns. PD-normalized columns are ignored. |
| | MaxQuant | `proteinGroups.txt`, then `evidence.txt`. Raw `Intensity` columns are used; `LFQ intensity` only if no raw columns exist. |
| | FragPipe | `combined_protein.tsv` and the run's `fragpipe.workflow`. The workflow file tells OSMSS whether IonQuant already normalized the data. |
| Label-free DIA | Spectronaut | Pivot Report (2D) at protein-group level (directly, or after export through Perseus), or a report made with **`OSMSS_Global.rs`** |
| | DIA-NN / FragPipe-DIA | `report.pg_matrix.tsv` (tick "Generate matrices" / `--matrices` in DIA-NN 1.x) |
| TMT (single or multi-plex) | Proteome Discoverer | Nested protein/peptide export (CSV), one per plex |
| | Other software | Flat peptide table with one reporter-intensity column per channel |
| Phospho, DIA | Spectronaut | Phospho (Site) Pivot Report, made with **`OSMSS_PTM.rs`**, plus a global protein pivot made with **`OSMSS_Global.rs`** |
| Phospho, TMT | Proteome Discoverer | Phospho-enriched and global exports for each plex |
| No replicates (N=1) | Proteome Discoverer | CSV with individual `Abundance` columns |
| FragPipe PTM database | FragPipe | `psm.tsv` and `peptide.tsv` |
| GSEA | Any | A pre-ranked `.rnk` file, or `.gct` + `.cls` files |
| Kinase-Site Analysis | OSMSS | The PhosR-format CSV exported from the PTM DE results |

**Current limitations:**
- PTM analysis supports **phosphorylation** only.
- Label-free DDA phospho input is not supported.
- MaxQuant TMT output is not read directly; export a flat peptide table instead.

---

## Spectronaut export schemes (.rs files)

Spectronaut's long-format reports can exceed the 1 GB upload limit. Two report schemes in this repository export exactly the columns OSMSS needs:

| File | Use it for |
|---|---|
| `OSMSS_Global.rs` | Protein-group pivot (`PG.Quantity` per run) for the DIA route, and the global file of the phospho DIA route |
| `OSMSS_PTM.rs` | Site-level phospho pivot (`PTM.Quantity` per run) for the phospho DIA route |

**To import a scheme:**
1. In Spectronaut, open the **Report** perspective.
2. Import the `.rs` file as a report scheme.
3. Select it and export the report for your experiment.

The schemes were created with **Spectronaut 20.5**. Other Spectronaut versions may rename report columns; if OSMSS rejects a report, check that the `PG.Quantity` / `PTM.Quantity` columns and the identifier columns listed in the upload wizard are present.

---

## Reference data you download yourself

Two resources used by the residue-prediction and Kinase-Site tools are licensed **CC BY-NC-SA 3.0** (non-commercial). OSMSS cannot redistribute them, so they are **not included** in this repository. When a module needs them, OSMSS shows a dialog with the download links and creates the destination folders for you.

### The Kinase Library

Download from the [Kinase Library GitHub repository](https://github.com/TheKinaseLibrary/kinase-library), folder `src/kinase_library/databases/matrices`:

```
data/matrices/ser_thr/ser_thr_all_log2_matrices.tsv
data/matrices/ser_thr/st_favorability.tsv
data/matrices/tyrosine/tyrosine_all_log2_matrices.tsv
```

Restart OSMSS after adding these files.

### PhosphoSitePlus

Download from the [PhosphoSitePlus downloads page](https://www.phosphosite.org/staticDownloads) (free registration required), and unzip:

```
data/Phosphorylation_site_dataset/
data/Kinase_Substrate_Dataset/
```

### Downloaded automatically

These need no action:
- **PhosR** reference data (kinase motifs and substrates) is extracted from the installed PhosR package on first use.
- **MSigDB gene sets** for GSEA are downloaded through `msigdbr` and cached in `data/GSEA/`.
- **UniProt reference proteomes** (human UP000005640, mouse UP000000589; about 8 MB each, CC BY 4.0) are downloaded from UniProt the first time a PhosR or KSEA export needs sequence windows, and kept in `data/KSEA/`. To work offline, place `UP000005640_9606.fasta` (or `.fasta.gz`) there yourself.

---

## SAINTexpress (AP-MS)

OSMSS scores affinity-purification data with **SAINTexpress 3.6.3** (intensity mode, `SAINTexpress-int`), by Hyungwon Choi, Damian Fermin and colleagues. SAINTexpress is a separate program distributed by its authors under the GNU GPL v3 or later, and **it is not included in this repository**. Install it once, as follows.

**1. Download.** Get `SAINTexpress_v3.6.3__2018-03-09` from the [SAINT project page on SourceForge](https://sourceforge.net/projects/saint-apms/files/) and unpack it.

**2. Place the folder.** Move the unpacked `SAINTexpress_v3.6.3__2018-03-09` folder into the `SAINT/` folder of OSMSS. Keep the folder name unchanged, and create `bin/` inside it if it does not exist.

**3. Put the binary in `bin/`:**

| Platform | What to do |
|---|---|
| Windows | Copy `Precompiled_binaries/Windows64/SAINTexpress-int.exe` into `bin/` |
| Linux | Copy `Precompiled_binaries/Linux64/SAINTexpress-int` into `bin/`, then make it executable (see below). Alternatively, run `make -j` in the folder, which builds the binaries into `bin/`. |
| macOS | No precompiled binary is provided upstream. Compile from source by running `make -j` in the folder (needs the Xcode command-line tools); the binaries are written to `bin/`. Compiling on macOS has not been tested by the OSMSS authors. Do **not** copy the Linux binary, as it will not run on macOS. |

Make the binary executable on Linux:

```bash
chmod +x SAINT/SAINTexpress_v3.6.3__2018-03-09/bin/SAINTexpress-int
```

**4. Check the result.** You should end up with:

```
OSMSS/
└── SAINT/
    └── SAINTexpress_v3.6.3__2018-03-09/
        └── bin/
            └── SAINTexpress-int        (Windows: SAINTexpress-int.exe)
```

**5. Restart OSMSS.** It checks for the binary at start-up. Only `SAINTexpress-int` is used; `SAINTexpress-spc` (spectral counts) is not needed.

**If SAINTexpress is not installed:**
- Every other part of OSMSS works.
- Running the SAINT engine shows "SAINTexpress binary not found".
- The other AP-MS options (proDA with negative, fraction or cross-condition controls) do not need SAINTexpress.

**Citing:** if you use the SAINT route, cite SAINTexpress (Teo et al., *J Proteomics* 2014). The generated methods text includes this citation.

---

## Saved Experiments

Under **Saved Experiments**, OSMSS saves a finished analysis as a `.zip` file so you or a colleague can reopen it, for example from a PRIDE deposition.

> **Security note:** Saved Experiments contain R data (`.rds`) files. Only open `.zip` files from sources you trust, and use R ≥ 4.4, which fixes a known vulnerability in loading R data files (CVE-2024-27322).

Reloading is currently supported for the protein routes and for the DIA-PTM and TMT-PTM routes.

---

## Data sent to external services

**No quantitative MS data leave your computer.** Only these annotation tools contact online services:

| Service | What is sent | Used by |
|---|---|---|
| UniProt | Protein accessions (ID mapping asks before sending your dataset's accessions); a single accession for annotations | Protein coverage, residue prediction, ID mapping |
| AlphaFold DB | A single UniProt accession | Structure view, solvent accessibility |
| Ensembl REST | A gene symbol, then the Ensembl IDs of its orthologs | Conservation across species |
| Scansite 4 | The protein sequence | Residue prediction |
| 14-3-3-Pred | The protein sequence | Residue prediction |
| MSigDB (via msigdbr) | Nothing user-specific; gene-set download only | GSEA |
| UniProt FTP | Nothing user-specific; one-time reference-proteome download | PhosR and KSEA exports |

---

## Running OSMSS on a server

OSMSS runs unchanged under [Shiny Server](https://posit.co/products/open-source/shiny-server/) or Posit Connect. To set it up:
- Copy the repository into the server's application directory.
- Install the packages for the user that runs Shiny Server.
- For AP-MS scoring, install SAINTexpress (see [SAINTexpress](#saintexpress-ap-ms)).

**Server notes:**
- Uploads go to the R session's temporary directory and are kept separate per session.
- Each worker process needs 200–400 MB; three workers are started per R process.
- Proxies in front of Shiny Server must allow 1 GB uploads.

---

## Troubleshooting

- **Logs:**
  - OSMSS writes a daily log to `logs/osmss_YYYY-MM-DD.log`, deleted automatically after two days.
  - The log records every module load, upload, analysis step and error.
  - Attach it when you report a problem.
- **"Package … is not installed":** run the installation command again, or install the named package with `BiocManager::install("<name>")`.
- **"SAINTexpress binary not found":** check the folder layout in [SAINTexpress](#saintexpress-ap-ms) step 4, then restart OSMSS.
- **A module says "Reference data needed":** see [Reference data you download yourself](#reference-data-you-download-yourself).
- **A Spectronaut file is refused as too large:** export with `OSMSS_Global.rs` or as a Pivot Report.
- **Bugs and questions:** open an issue on the GitHub repository.

---

## Developer

**Erez Cohen, PhD**\
Research Investigator, Coulombe Lab\
Department of Cell & Developmental Biology\
University of Michigan Medical School, Ann Arbor, MI\
[erezco@umich.edu](mailto:erezco@umich.edu)

2026

---

## Citing OSMSS

If you use OSMSS, please cite:

> Cohen, E.; Coulombe, P. A. One Stop MS Shop (OSMSS): statistical analysis and interpretation of quantitative proteomics data for non-coders. *bioRxiv* **2026**. <!-- TODO: DOI -->

Please also cite the methods OSMSS ran for you. The generated methods text for each analysis lists them with versions:
- proDA
- limma
- SAINTexpress
- QFeatures/MsCoreUtils
- PhosR
- clusterProfiler/MSigDB
- The Kinase Library
- Scansite
- 14-3-3-Pred

Software archive: <!-- TODO: Zenodo DOI for the release -->

---

## License

OSMSS is free software: you can redistribute it and/or modify it under the terms of the **GNU General Public License, version 3 or (at your option) any later version** (SPDX: `GPL-3.0-or-later`); see `LICENSE`. OSMSS is distributed WITHOUT ANY WARRANTY; see the license for details.

Third-party components keep their own licenses:

| Component | License | Distributed here? |
|---|---|---|
| Kinase family table (`data/Kinase_Families.tsv`) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); source articles © The Author(s) 2023 and 2024. Adapted from Johnson JL et al. (2023) *Nature* 613:759–766, Supplementary Table 1 ([doi:10.1038/s41586-022-05575-3](https://doi.org/10.1038/s41586-022-05575-3)); Yaron-Barir TM et al. (2024) *Nature* 629:1174–1181, Supplementary Table 1 ([doi:10.1038/s41586-024-07407-y](https://doi.org/10.1038/s41586-024-07407-y)); and the UniProt Consortium's `pkinfam.txt` and UniProtKB protein families ([UniProt](https://www.uniprot.org), CC BY 4.0). Changes: the three sources were merged into one group/family table, and group names were harmonized to the Manning (2002) convention; `data/Kinase_Families_build.R` lists every change and rebuilds the table from the original files. | Yes |
| SAINTexpress 3.6.3 | GPL-3.0-or-later | No; installed by the user |
| The Kinase Library matrices | CC BY-NC-SA 3.0 | No; downloaded by the user |
| PhosphoSitePlus datasets | CC BY-NC-SA 3.0 | No; downloaded by the user |
| PhosR reference data | From the PhosR package (includes PhosphoSitePlus-derived substrates) | No; extracted at run time |
| MSigDB gene sets | MSigDB terms; KEGG-derived collections have their own terms | No; downloaded at run time |
| R packages | Their respective licenses | No; installed from CRAN/Bioconductor |
