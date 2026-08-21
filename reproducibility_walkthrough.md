# Reproducibility Walkthrough

This document provides a guide for reproducing the results reported in **Mansour L. *et al.* (2026). Spectral Normative Modeling of Brain Structure.** It is intended as the primary entry point for anyone wishing to verify, reproduce, or build on the manuscript's findings.

The walkthrough is organised into five sections:

1. **Prerequisites**: system tools and access requirements before installing the Python environment.
2. **Setting up the environment**: installing dependencies and verifying the installation.
3. **Runtime expectations**: approximate runtimes and hardware requirements for each part of the pipeline.
4. **Figure generation walkthrough**: a figure-by-figure map linking each manuscript result to the notebook (and where possible, the specific cell) that produces it.
5. **If something doesn't work**: common issues and where to find help.

Throughout this document, references to manuscript figures use the form *Fig. N* for main-text figures and *Ext. Data Fig. E.N* for extended data figures.

---

## 1. Prerequisites

### System tools

In addition to the Python environment described in the next section, the following system tools must be installed.

**LaTeX** is required for figure rendering via matplotlib's `text.usetex` setting. The figures use Computer Modern fonts and standard scientific LaTeX packages.

- *Ubuntu/Debian:* `sudo apt install texlive-latex-extra texlive-fonts-recommended cm-super dvipng`
- *macOS:* Install [MacTeX](https://www.tug.org/mactex/) (full distribution), or [BasicTeX](https://www.tug.org/mactex/morepackages.html) followed by `sudo tlmgr install cm-super type1cm dvipng standalone`.
- *Windows:* Install [MiKTeX](https://miktex.org/), which auto-installs missing packages on first use.

To verify the installation works with matplotlib:

```python
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.figure(); plt.text(0.5, 0.5, r'$\alpha + \beta = \gamma$'); plt.show()
```

If this raises a `RuntimeError` mentioning `latex` or `dvipng`, the LaTeX installation is incomplete.

**git-annex** is required by [datalad](https://www.datalad.org/) for fetching the public datasets provided by the [Reproducible Brain Charts platform](https://reprobrainchart.github.io/) used in the early pipeline notebooks (stages 01–02). It is included as part of the conda environment described below, so most users do not need to install it separately. If you encounter issues with the conda-provided version, it can also be installed via the system package manager:

- *Ubuntu/Debian:* `sudo apt install git-annex`
- *macOS:* `brew install git-annex`
- *Windows:* See [the git-annex Windows installation guide](https://git-annex.branchable.com/install/Windows/).

### Data access

The manuscript integrates over 30 publicly available neuroimaging datasets (see manuscript Section S.1.1, "Healthy Lifespan Imaging Data"), plus the MACC clinical cohort used in Sections 2.5–2.6. Each dataset is subject to its own access and usage agreements. While some datasets are openly available, others require approval through data access committees or the execution of material transfer or data-use agreements. Researchers interested in accessing these datasets should consult the corresponding data providers for specific access procedures.

Note that **no raw imaging data is redistributed through this repository**. Researchers wishing to re-execute the complete pipeline from scratch should obtain each dataset through its conventional access route. Where applicable, the figure-generating notebooks operate on outputs of the trained model — provided as part of the [pretrained release](https://github.com/sina-mansour/normative_brain_charts/releases/latest) — and can therefore be reproduced without access to the underlying imaging data. Reproducing the full training pipeline from raw data is supported but requires substantial compute resources (see Section 3).

Even for users who do not obtain the underlying data, the notebooks in this repository serve as a transparent record of the computational pipeline behind every manuscript figure. Each analysis step, from preprocessing through model fitting, evaluation, and figure generation, is documented in executable form alongside its corresponding outputs, allowing the methodology to be audited and adapted without rerunning the pipeline end-to-end.

---

## 2. Setting up the environment

Two environment files are provided at the repository root.

**`environment.yml`** is the recommended file for most users. It lists direct dependencies pinned to the versions used to produce the manuscript. To install:

```bash
conda env create -f environment.yml
conda activate snm_env
```

The installation takes less than 5 minutes on a typical laptop with a reasonable internet connection.

Note: **`environment-full.yml`** contains the complete dependency state of the environment used to produce the manuscript, including all transitive dependencies. This file is provided for forensic reproducibility and reflects a Linux installation; it will not resolve cleanly on macOS or Windows. Use `environment.yml` for cross-platform installation.

### Verifying the installation

Once the environment is active, the following imports should all succeed:

```python
import numpy, scipy, pandas, sklearn, matplotlib, seaborn
import nibabel, nilearn
import pymc, pytensor, arviz, xarray, h5py
import spectranorm, neuromaps
import Connectome_Spatial_Smoothing, cerebro
import gdist, trimesh, diptest, statsmodels
print('All critical imports succeeded')
```

If any import fails, see Section 5.

### Downloading the pretrained model

The pretrained SNM-1000 model and derived normative charts are available as a release asset. To download:

```bash
wget https://github.com/sina-mansour/normative_brain_charts/releases/download/v1.0/pretrained_SNM_1000_V1.0.tar.gz
tar -xzvf pretrained_SNM_1000_V1.0.tar.gz
```

The extracted directory contains the model parameters needed by the figure-generating notebooks (Section 4).

---

## 3. Runtime expectations

The pipeline is organised into nine sequential stages, with computational costs that vary by several orders of magnitude. The table below summarises approximate runtime and hardware requirements for each stage.

*[TODO: fill in the table once you have measured runtimes. Indicative entries below.]*

| Stage | Description | Typical runtime | Hardware |
|-------|-------------|-----------------|----------|
| 01 | Data import (per dataset) | 10 min – several hours | Laptop or server (some datasets require HPC for FreeSurfer processing) |
| 02 | Data aggregation | ~30 min | Laptop |
| 03 | Direct normative model fit (sanity check) | ~2 hours | Server with ≥16 GB RAM |
| 04 | Connectome eigenmode construction | ~15 min | Laptop |
| 05 | Spectral normative model fit (k = 1000) | ~6 hours | HPC node with ≥32 GB RAM |
| 05 | Spectral normative model fit (k = 10 000) | ~9 hours | HPC node with ≥64 GB RAM |
| 06 | Performance evaluations | ~1 hour | Laptop (uses cached outputs) |
| 07 | Growth gradient analyses | ~30 min | Laptop |
| 08 | Clinical evaluation (AD) | ~30 min | Laptop |
| 09 | Data sharing and chart export | ~15 min | Laptop |

For the figure-generating notebooks specifically (stages 06–09), runtimes are tractable on a laptop because they operate on the pretrained model rather than retraining it. Reproducing the training pipeline from scratch (stages 02–05) requires access to high-performance computing infrastructure.

---

## 4. Figure generation walkthrough

This section maps each manuscript figure to the notebook(s) that produce it. For each figure, the table provides:

- **Notebook**: the source notebook in the repository.
- **Rendered view**: an HTML rendering of the notebook served via GitHub Pages, with anchor links to the relevant section. This allows the figure-generating code and its outputs to be viewed in a browser without local execution.
- **Dependencies**: any upstream notebooks whose outputs are required.
- **Runtime**: approximate execution time, assuming the dependencies are already satisfied.

### Main-text figures

*[TODO: fill in cells, anchor URLs, and runtimes. Indicative structure below.]*

| Figure | Description | Notebook | Rendered view | Dependencies | Runtime |
|--------|-------------|----------|---------------|--------------|---------|
| Fig. 1 | SNM framework schematic | (illustrative; not produced by code) | — | — | — |
| Fig. 2 | SNM vs direct model benchmark | [`06_01_systematic_performance_evaluations.ipynb`](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb) | [HTML §Fig. 2](TODO) | Pretrained SNM-1000 + direct models (stage 05) | ~30 min from cached outputs |
| Fig. 3 | High-resolution lifespan charts | [`06_02_publication_figures.ipynb`](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb) | [HTML §Fig. 3](TODO) | Pretrained SNM-1000 | ~10 min |
| Fig. 4 | Lifespan thickness growth gradients | [`07_01_growth_gradients.ipynb`](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb) | [HTML §Fig. 4](TODO) | Pretrained SNM-1000 | ~15 min |
| Fig. 5 | AD cortical atrophy signature | [`08_01_clinical_application_in_AD.ipynb`](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb) | [HTML §Fig. 5](TODO) | Pretrained SNM-1000 + MACC data | ~20 min |
| Fig. 6 | AD heterogeneity landscape | [`08_01_clinical_application_in_AD.ipynb`](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb) | [HTML §Fig. 6](TODO) | Pretrained SNM-1000 + MACC data | ~15 min |

### Extended Data figures

*[TODO: fill in. Below are placeholders for the figures discussed in the main text.]*

| Figure | Description | Notebook | Rendered view | Dependencies | Runtime |
|--------|-------------|----------|---------------|--------------|---------|
| Ext. Data Fig. E.1 | Eigenmode reconstruction accuracy | TODO | TODO | TODO | TODO |
| Ext. Data Fig. E.2–E.4 | Probabilistic gradient strata construction | TODO | TODO | TODO | TODO |
| Ext. Data Fig. E.5–E.7 | Normative trajectories along TGG1–3 | TODO | TODO | TODO | TODO |
| Ext. Data Fig. E.8 | Developmental pruning vs aging decline | TODO | TODO | TODO | TODO |

### Supplementary figures

*[TODO: include key supplementary figures, particularly those referenced in the response to reviewers. At minimum: Figs. S.1 (QC), S.2 (split distributions), S.3 (cross-basis correlations), S.5–S.9 (regional performance), and S.17–S.20 (within-cohort cognitive associations and heterogeneity robustness).]*

| Figure | Description | Notebook | Rendered view | Dependencies | Runtime |
|--------|-------------|----------|---------------|--------------|---------|
| TODO | TODO | TODO | TODO | TODO | TODO |

---

## 5. If something doesn't work

### Environment installation fails

If `conda env create` fails with a `PackagesNotFoundError`, the most likely cause is a platform mismatch; `environment-full.yml` is Linux-specific and will not resolve on macOS or Windows. Use `environment.yml` instead.

If `conda env create` fails with a `LinkError` mentioning `git-annex` or compiler tools, the conda-provided git-annex is failing to install. Either install `git-annex` separately via your system package manager and remove it from the environment file, or install the `gxx_linux-64` package into your base conda environment to provide the required toolchain.

### Figures look different from the manuscript

If figures render but text appears in a sans-serif font rather than Computer Modern, LaTeX is not being used. Verify your LaTeX installation as described in Section 1.

If figures are missing or contain only axes with no data, the most likely cause is a missing intermediate file. Check the "Dependencies" column for the relevant figure in Section 4 and run the upstream notebooks first.

### Where to get help

For issues specific to this repository, open an issue at [github.com/sina-mansour/normative_brain_charts/issues](https://github.com/sina-mansour/normative_brain_charts/issues).

For questions about the underlying SNM framework, see the [spectranorm documentation](https://sina-mansour.github.io/spectranorm/).

For correspondence about the manuscript itself, contact the corresponding authors as listed in the published paper.

[![Contact](https://img.shields.io/badge/Contact-sina.mansour.lakouraj@gmail.com-blue)](https://sina-mansour.github.io/)