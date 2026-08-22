# Reproducibility Walkthrough

This document provides a guide for reproducing the results reported in **Mansour L. *et al.* (2026). Spectral Normative Modeling of Brain Structure.** It is intended as the primary entry point for anyone wishing to verify, reproduce, or build on the manuscript's findings.

The walkthrough is organised into five sections:

1. 🧰 **Prerequisites**: system tools and access requirements before installing the Python environment.
2. 📦 **Setting up the environment**: installing dependencies and verifying the installation.
3. ⏱️ **Runtime expectations**: approximate runtimes and hardware requirements for each part of the pipeline.
4. 🖼️ **Figure generation walkthrough**: a figure-by-figure map linking each manuscript result to the notebook (and where possible, the specific cell) that produces it.
5. 🛠️ **If something doesn't work**: common issues and where to find help.

Throughout this document, references to manuscript figures use the form *Fig. N* for main-text figures, *Ext. Data Fig. E.N* for extended data figures, and *Fig. S.N* for supplementary figures.

---

## 1. 🧰 Prerequisites

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

Note that **no raw imaging data is redistributed through this repository**. Researchers wishing to re-execute the complete pipeline from scratch should obtain each dataset through its conventional access route. Where applicable, the figure-generating notebooks operate on outputs of the trained model, provided as part of the [pretrained release](https://github.com/sina-mansour/normative_brain_charts/releases/latest), and can therefore be reproduced without access to the underlying imaging data. Reproducing the full training pipeline from raw data is supported but requires substantial compute resources (see Section 3).

Even for users who do not obtain the underlying data, the notebooks in this repository serve as a transparent record of the computational pipeline behind every manuscript figure. Each analysis step, from preprocessing through model fitting, evaluation, and figure generation, is documented in executable form alongside its corresponding outputs, allowing the methodology to be audited and adapted without rerunning the pipeline end-to-end.

---

## 2. 📦 Setting up the environment

Two environment files are provided at the repository root, serving different purposes.

**`environment.yml`** is all that most users will need. It lists the direct dependencies, pinned to the versions used to produce the manuscript, and is the recommended starting point on any platform:

```bash
conda env create -f environment.yml
conda activate snm_env
```

Installation takes under five minutes on a typical laptop with a reasonable internet connection.

**`environment-full.yml`** is the complete record, for anyone who wants it: every package present in the environment used to produce the manuscript, transitive dependencies included, each pinned to an exact version. Build strings have been stripped from the export, so conda is free to select an appropriate build for the current platform. A small number of the pinned versions were resolved on Linux and may not have counterparts elsewhere, which is why `environment.yml` remains the recommended route for cross-platform installation.

### Verifying the installation

Once the environment is active, the following imports should all succeed:

```python
import numpy, scipy, pandas, sklearn, matplotlib, seaborn
import nibabel, nilearn
import pymc, pytensor, arviz, xarray, h5py
import spectranorm, neuromaps, brainsmash
import Connectome_Spatial_Smoothing, cerebro
import gdist, trimesh, diptest, statsmodels
import joblib, tqdm, pyarrow, skimage
print('All critical imports succeeded')
```

### The R environment for SigClust2

Nearly all analyses run in the Python environment described above. The one exception is the SigClust2 hierarchical significance test used in the AD clusterability analyses, which is implemented only as an R package and requires a small, separate conda environment. Setup instructions, the pinned package version, and the command to run the analysis are given in 📓 [`08_02_clinical_AD_cluster_analyses.ipynb` cell 53](code/notebooks/08_clinical_evaluation/08_02_clinical_AD_cluster_analyses.ipynb#nb-cell-53).

### Downloading the pretrained model

The pretrained SNM-1000 model and derived normative charts are available as a release asset. To download:

```bash
wget https://github.com/sina-mansour/normative_brain_charts/releases/download/v1.0/pretrained_SNM_1000_V1.0.tar.gz
tar -xzvf pretrained_SNM_1000_V1.0.tar.gz
```

The extracted directory contains the model parameters needed by the figure-generating notebooks (Section 4).

---

## 3. ⏱️ Runtime expectations

The pipeline spans computational costs that vary by several orders of magnitude, from figure-generating notebooks that run in minutes to model fits that can take nearly a day. The table below gives approximate runtimes and indicates the hardware each stage realistically requires.

All timings were measured on a compute server with 64 vCPU cores and 251 GB of RAM, and are intended as general guidance rather than benchmarks. Actual runtimes will vary with hardware, parallelism, and parameter choices. Most steps used only a small fraction of the available memory; the notable exception is eigenmode extraction (stage 04), which is memory-intensive.

**Hardware categories**

| | Meaning |
|---|---|
| 🔴 | Requires a compute server, due to memory, storage, or parallelism demands |
| 🟠 | Benefits substantially from a server, but can be run on a single machine |
| 🟢 | Runs comfortably on a personal machine of the kind commonly used for neuroimaging research |

| Stage | Description | Runtime | Hardware |
|-------|-------------|---------|----------|
| 01 | Data import and preprocessing (27 notebooks) | 1–2 days | 🔴 |
| 02 | Data aggregation and spectral encoding | 2–4 hours | 🟠 |
| 03 | Direct normative model (single phenotype) | ~2 min | 🟢 |
| 04 | Connectome eigenmode construction | ~12 h | 🔴 |
| 05 | Spectral normative model fit (*k* = 1000) | ~8 h | 🟠 |
| 05 | Spectral normative model fit (*k* = 10 000) | ~12 h | 🔴 |
| 06 | Systematic performance evaluations | ~20 h | 🔴 |
| 06 | Publication figures | ~20 min | 🟢 |
| 07 | Growth gradient analyses | ~1 h | 🟠 |
| 08 | Clinical model fine-tuning | ~11 h | 🟠 |
| 08 | Individual deviation mapping and analyses | ~1 h | 🟢 |
| 08 | Supplementary video generation | ~20 min | 🟢 |
| 08 | AD clusterability analyses | ~40 min | 🟢 |
| 09 | Pretrained model and atlas charts | ~2 min | 🟢 |
| 09 | Vertex-wise charts | ~17 min | 🟢 |

**Stage 01.** The figure above assumes preprocessed cortical thickness data is already available. It excludes the considerable time required to obtain data access approvals and to run FreeSurfer preprocessing, neither of which is performed within these notebooks.

**Stage 05.** Runtime scales with the number of retained eigenmodes and with model configuration. The 1000-mode model completes in less than a single working day and remains largely feasible on a personal machine, demonstrating the overall efficiency of the framework. The 10 000-mode model is substantially more demanding and should ideally be run on a capable compute server.

**Stage 06.** The bulk of the time is spent on systematic benchmarking across mode counts and spatial granularities, parallelised across cores. The publication figures are generated from cached outputs and their respective cells run comparably faster.

**Stage 08.** Only the transfer-learning step that adapts the pretrained model to the clinical site is computationally demanding; it accounts for approximately 11 of the stage's 13 hours. Everything downstream of that fit, including individual deviation mapping, figure generation, and the post-hoc analyses, runs considerably faster. The supplementary videos are generated by three encoding cells that a reader reproducing only the figures can skip.

Notebooks that operate on the pretrained model rather than retraining it (stages 06 figures, 07, 08 post-fit analyses, and 09) are tractable without a compute server. Reproducing the training pipeline end-to-end (stages 01–06) requires high-performance computing infrastructure.

---

## 4. 🖼️ Figure generation walkthrough

This section maps each manuscript figure to the notebook, and where possible the specific cell, that produces it. For each figure, the table provides:

- 📓 **Notebook**: the source notebook in the repository.
- 🖼️ **Rendered view**: a link to the output of the cell in the HTML rendering of that notebook, served via GitHub Pages. This allows each figure and the code that produced it to be inspected in a browser without local execution.
- 🔗 **Dependencies**: upstream notebooks or cells whose outputs are required.
- ⏱️ **Runtime**: approximate execution time for the cells listed under Rendered view.

Cells are referred to by the notebook they belong to, abbreviated by its numeric prefix. A reference to `07_01` cell 112, for example, points to the 112th cell of 📓 `07_01_growth_gradients.ipynb`, counting markdown and code cells together from the top.

Several manuscript figures are composites. The notebooks produce individual panels, and the final figure is assembled from those panels, in some cases within a later notebook cell and in others manually in Inkscape. Where a figure is assembled manually, the table links to the cells that produce its constituent panels and names the panel each corresponds to.

The Dependencies column lists the immediate upstream requirements only. Most of these have dependencies of their own, so a reader reproducing the pipeline from scratch should run the notebooks in numerical order and execute the cells within each notebook sequentially from top to bottom, rather than attempting to satisfy individual dependencies in isolation.

Runtimes are banded and refer only to the compute time of the cells listed under Rendered view, measured on the hardware described in Section 3. They exclude upstream dependencies, and they also exclude the iterative development and manual styling that preceded each final figure. Figure 1 illustrates the gap: its programmatically generated panel runs in seconds, while the two schematic panels beside it were assembled by hand over several hours.

### Main-text figures

| Figure | Description | 📓 Notebook | 🖼️ Rendered view | 🔗 Dependencies | ⏱️ Runtime |
|--------|-------------|----------|---------------|--------------|---------|
| Fig. 1 | SNM framework schematic and multi-scale charting overview | Panel A from 📓 [`02_01_aggregating_imported_datasets.ipynb`](code/notebooks/02_data_aggregation/02_01_aggregating_imported_datasets.ipynb); panels B and C are schematic diagrams assembled manually, drawing on elements from 📓 `06_02_publication_figures.ipynb` | Panel A: [`02_01` cell 20](code/notebooks/02_data_aggregation/02_01_aggregating_imported_datasets.ipynb#nb-cell-20-output) | Stage 01 import notebooks | < 1 min |
| Fig. 2 | SNM normative performance benchmarked against the direct approach across *k* | 📓 [`06_02_publication_figures.ipynb`](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb), assembled manually with panels from 📓 [`06_01_systematic_performance_evaluations.ipynb`](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb) | Composite: [`06_02` cell 24](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb#nb-cell-24-output); top panel: [`06_02` cell 12](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb#nb-cell-12-output); reconstruction panels: [`06_01` cell 161](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-161-output) | [`06_01` cell 115](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-115), [`06_01` cell 117](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-117) | 2–5 min |
| Fig. 3 | SNM-derived high-resolution charts of cortical thickness across the lifespan | 📓 [`06_02_publication_figures.ipynb`](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb) | [`06_02` cell 34](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb#nb-cell-34-output) | [`06_01` cell 190](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-190), [`06_01` cell 193](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-193) | < 1 min |
| Fig. 4 | Lifespan thickness growth gradients (TGG1–3) | 📓 [`07_01_growth_gradients.ipynb`](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb) | [`07_01` cell 112](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-112-output) | [`07_01` cell 66](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-66), [`07_01` cell 90](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-90), [`07_01` cell 93](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-93), [`07_01` cell 94](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-94) | < 1 min |
| Fig. 5 | Thickness growth gradients recapitulate established cortical hierarchies | 📓 [`07_01_growth_gradients.ipynb`](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb) | [`07_01` cell 113](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-113-output) | [`07_01` cell 97](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-97), [`07_01` cell 99](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-99), [`07_01` cell 102](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-102), [`07_01` cell 110](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-110) | < 1 min |
| Fig. 6 | Cortical signature of atrophy in Alzheimer's disease and its cognitive correlates | 📓 [`08_01_clinical_application_in_AD.ipynb`](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb) and 📓 [`06_02_publication_figures.ipynb`](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb), assembled manually in Inkscape | Panel A: [`06_02` cell 14](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb#nb-cell-14-output); Panel B: [`08_01` cell 96](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-96-output); Panel C: [`08_01` cell 56](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-56-output); Panel D: [`08_01` cell 52](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-52-output); Panel E: [`08_01` cell 54](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-54-output) | `08_01` fine-tuned clinical model (stage 08) | 1–2 min |
| Fig. 7 | Heterogeneity landscape of atrophy in AD | 📓 [`08_01_clinical_application_in_AD.ipynb`](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb), assembled manually in Inkscape | Panel A: [`08_01` cell 161](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-161-output) to [cell 178](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-178-output); Panel B: [`08_01` cell 99](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-99-output) to [cell 102](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-102-output); Panel C: [`08_01` cell 109](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-109-output); Panel D: [`08_01` cell 118](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-118-output) and [cell 120](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-120-output) to [cell 131](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-131-output) | `08_01` fine-tuned clinical model and distance matrix (stage 08) | 5–10 min |

### Extended Data figures

| Figure | Description | 📓 Notebook | 🖼️ Rendered view | 🔗 Dependencies | ⏱️ Runtime |
|--------|-------------|----------|---------------|--------------|---------|
| Ext. Data Fig. E.1 | Eigenmode-based reconstruction of cortical thickness and spatial query maps | 📓 [`06_02_publication_figures.ipynb`](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb) | [`06_02` cell 18](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb#nb-cell-18-output) | [`06_01` cell 150](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-150), [`06_01` cell 154](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-154), [`06_01` cell 161](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-161) | 1–2 min |
| Ext. Data Fig. E.2 | Construction of probabilistic gradient strata along TGG1 | 📓 [`07_01_growth_gradients.ipynb`](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb) | [`07_01` cell 76](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-76-output) | [`07_01` cell 66](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-66), [`07_01` cell 72](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-72) | < 1 min |
| Ext. Data Fig. E.3 | Construction of probabilistic gradient strata along TGG2 | 📓 [`07_01_growth_gradients.ipynb`](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb) | [`07_01` cell 77](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-77-output) | [`07_01` cell 66](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-66), [`07_01` cell 72](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-72) | < 1 min |
| Ext. Data Fig. E.4 | Construction of probabilistic gradient strata along TGG3 | 📓 [`07_01_growth_gradients.ipynb`](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb) | [`07_01` cell 78](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-78-output) | [`07_01` cell 66](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-66), [`07_01` cell 72](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-72) | < 1 min |
| Ext. Data Fig. E.5 | Normative cortical thickness trajectories along TGG1 | 📓 [`07_01_growth_gradients.ipynb`](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb) | [`07_01` cell 84](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-84-output) | [`07_01` cell 66](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-66), [`07_01` cell 72](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-72) | < 1 min |
| Ext. Data Fig. E.6 | Normative cortical thickness trajectories along TGG2 | 📓 [`07_01_growth_gradients.ipynb`](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb) | [`07_01` cell 86](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-86-output) | [`07_01` cell 66](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-66), [`07_01` cell 72](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-72) | < 1 min |
| Ext. Data Fig. E.7 | Normative cortical thickness trajectories along TGG3 | 📓 [`07_01_growth_gradients.ipynb`](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb) | [`07_01` cell 88](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-88-output) | [`07_01` cell 66](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-66), [`07_01` cell 72](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-72) | < 1 min |
| Ext. Data Fig. E.8 | Replication of the TGG1–thickness correspondence across age-restricted cohorts | 📓 [`07_01_growth_gradients.ipynb`](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb) | [`07_01` cell 118](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-118-output) | [`07_01` cell 117](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-117) | < 1 min |
| Ext. Data Fig. E.9 | Relationship between developmental cortical pruning and aging-related thickness decline | 📓 [`07_01_growth_gradients.ipynb`](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb) | [`07_01` cell 121](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-121-output) | [`07_01` cell 66](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-66) | < 1 min |

### Supplementary figures

| Figure | Description | 📓 Notebook | 🖼️ Rendered view | 🔗 Dependencies | ⏱️ Runtime |
|--------|-------------|----------|---------------|--------------|---------|
| Fig. S.1 | Quality control of cortical surface reconstruction (rescaled Euler numbers) | 📓 [`02_01_aggregating_imported_datasets.ipynb`](code/notebooks/02_data_aggregation/02_01_aggregating_imported_datasets.ipynb) | [`02_01` cell 19](code/notebooks/02_data_aggregation/02_01_aggregating_imported_datasets.ipynb#nb-cell-19-output) | Stage 01 import notebooks | < 1 min |
| Fig. S.2 | Distribution of demographic information before and after splitting | 📓 [`02_02_train_test_split_labels.ipynb`](code/notebooks/02_data_aggregation/02_02_train_test_split_labels.ipynb) | [`02_02` cell 18](code/notebooks/02_data_aggregation/02_02_train_test_split_labels.ipynb#nb-cell-18-output) | [`02_01` cell 20](code/notebooks/02_data_aggregation/02_01_aggregating_imported_datasets.ipynb#nb-cell-20) | < 1 min |
| Fig. S.3 | Sparse cross-basis correlation structure | 📓 [`06_01_systematic_performance_evaluations.ipynb`](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb) | [`06_01` cell 143](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-143-output) | [`06_01` cell 138](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-138) | < 1 min |
| Fig. S.4 | Cortical projection of sampled high-resolution vertices | 📓 [`06_01_systematic_performance_evaluations.ipynb`](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb) | [`06_01` cell 73](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-73-output) | [`06_01` cell 72](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-72) | < 1 min |
| Fig. S.5 | Regional variability of reconstruction accuracy across spatial query families | 📓 [`06_02_publication_figures.ipynb`](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb) | [`06_02` cell 19](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb#nb-cell-19-output) | [`06_01` cell 158](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-158) | < 1 min |
| Fig. S.6 | Signal reconstruction residuals | 📓 [`06_02_publication_figures.ipynb`](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb) | [`06_02` cell 20](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb#nb-cell-20-output) | [`06_01` cell 161](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-161) | 1–2 min |
| Fig. S.7 | Normative performance for brain-wide spatial queries | 📓 [`06_02_publication_figures.ipynb`](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb) | [`06_02` cell 25](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb#nb-cell-25-output) | [`06_01` cell 113](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-113), [`06_01` cell 124](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-124) | 1–2 min |
| Fig. S.8 | Normative performance for regional spatial queries (Yan-200) | 📓 [`06_02_publication_figures.ipynb`](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb) | [`06_02` cell 26](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb#nb-cell-26-output) | [`06_01` cell 113](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-113), [`06_01` cell 127](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-127) | 1–2 min |
| Fig. S.9 | Normative performance for high-resolution spatial queries (8 mm FWHM) | 📓 [`06_02_publication_figures.ipynb`](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb) | [`06_02` cell 27](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb#nb-cell-27-output) | [`06_01` cell 113](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-113), [`06_01` cell 133](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-133) | < 1 min |
| Fig. S.10 | Signal reconstruction accuracies across granularity scales | 📓 [`06_02_publication_figures.ipynb`](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb) | [`06_02` cell 22](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb#nb-cell-22-output) | [`06_01` cell 150](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-150), [`06_01` cell 154](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-154) | 1–2 min |
| Fig. S.11 | Sensitivity of regional queries to granularity scale | 📓 [`06_02_publication_figures.ipynb`](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb) | [`06_02` cell 31](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb#nb-cell-31-output) | [`06_01` cell 113](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-113), [`06_01` cell 127](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-127) | 1–2 min |
| Fig. S.12 | Sensitivity of high-resolution queries to granularity scale | 📓 [`06_02_publication_figures.ipynb`](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb) | [`06_02` cell 32](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb#nb-cell-32-output) | [`06_01` cell 113](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-113), [`06_01` cell 133](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-133) | 1–2 min |
| Fig. S.13 | Signal reconstruction accuracies for average versus asymmetric signals | 📓 [`06_02_publication_figures.ipynb`](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb) | [`06_02` cell 21](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb#nb-cell-21-output) | [`06_01` cell 150](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-150), [`06_01` cell 154](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-154) | 2–5 min |
| Fig. S.14 | Sensitivity to lateralized brain-wide spatial queries | 📓 [`06_02_publication_figures.ipynb`](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb) | [`06_02` cell 28](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb#nb-cell-28-output) | [`06_01` cell 113](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-113), [`06_01` cell 124](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-124) | 2–5 min |
| Fig. S.15 | Sensitivity to lateralized regional spatial queries | 📓 [`06_02_publication_figures.ipynb`](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb) | [`06_02` cell 29](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb#nb-cell-29-output) | [`06_01` cell 113](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-113), [`06_01` cell 127](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-127) | 1–2 min |
| Fig. S.16 | Sensitivity to lateralized high-resolution spatial queries | 📓 [`06_02_publication_figures.ipynb`](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb) | [`06_02` cell 30](code/notebooks/06_performance_evaluation/06_02_publication_figures.ipynb#nb-cell-30-output) | [`06_01` cell 113](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-113), [`06_01` cell 133](code/notebooks/06_performance_evaluation/06_01_systematic_performance_evaluations.ipynb#nb-cell-133) | 1–2 min |
| Fig. S.17 | Extended spatial correspondence between TGG1–3 and 81 cortical maps | 📓 [`07_01_growth_gradients.ipynb`](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb) | [`07_01` cell 108](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-108-output) | [`07_01` cell 101](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-101), [`07_01` cell 105](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb#nb-cell-105) | 1–2 min |
| Fig. S.18 | Cognitive association tests within the healthy cohort (HC) | 📓 [`08_01_clinical_application_in_AD.ipynb`](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb) | [`08_01` cell 88](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-88-output) to [cell 91](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-91-output) | `08_01` fine-tuned clinical model (stage 08) | 1–2 min |
| Fig. S.19 | Cognitive association tests within the MCI cohort | 📓 [`08_01_clinical_application_in_AD.ipynb`](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb) | [`08_01` cell 81](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-81-output) to [cell 85](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-85-output) | `08_01` fine-tuned clinical model (stage 08) | 1–2 min |
| Fig. S.20 | Cognitive association tests within the AD cohort | 📓 [`08_01_clinical_application_in_AD.ipynb`](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb) | [`08_01` cell 75](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-75-output) to [cell 78](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-78-output) | `08_01` fine-tuned clinical model (stage 08) | 1–2 min |
| Fig. S.21 | Cognitive association tests with a coarse brain-wide normative metric | 📓 [`08_01_clinical_application_in_AD.ipynb`](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb) | [`08_01` cell 182](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-182-output) to [cell 185](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-185-output) | `08_01` fine-tuned clinical model (stage 08) | < 1 min |
| Fig. S.22 | Robustness of the heterogeneity landscape to dimensionality reduction choices | 📓 [`08_01_clinical_application_in_AD.ipynb`](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb) | [`08_01` cell 152](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-152-output) and [cell 157](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-157-output) | [`08_01` cell 109](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb#nb-cell-109) | 1–2 min |
| Fig. S.23 | Cluster tendency of AD deviation patterns (gap statistic) | 📓 [`08_02_clinical_AD_cluster_analyses.ipynb`](code/notebooks/08_clinical_evaluation/08_02_clinical_AD_cluster_analyses.ipynb) | [`08_02` cell 41](code/notebooks/08_clinical_evaluation/08_02_clinical_AD_cluster_analyses.ipynb#nb-cell-41-output) | [`08_02` cell 33](code/notebooks/08_clinical_evaluation/08_02_clinical_AD_cluster_analyses.ipynb#nb-cell-33); `08_01` deviation maps and MDS landscape | < 1 min |
| Fig. S.24 | Stability of candidate partitions under resampling | 📓 [`08_02_clinical_AD_cluster_analyses.ipynb`](code/notebooks/08_clinical_evaluation/08_02_clinical_AD_cluster_analyses.ipynb) | [`08_02` cell 47](code/notebooks/08_clinical_evaluation/08_02_clinical_AD_cluster_analyses.ipynb#nb-cell-47-output) | `08_01` deviation maps and MDS landscape | < 1 min |
| Fig. S.25 | Generalizability of candidate partitions (prediction strength) | 📓 [`08_02_clinical_AD_cluster_analyses.ipynb`](code/notebooks/08_clinical_evaluation/08_02_clinical_AD_cluster_analyses.ipynb) | [`08_02` cell 50](code/notebooks/08_clinical_evaluation/08_02_clinical_AD_cluster_analyses.ipynb#nb-cell-50-output) | `08_01` deviation maps and MDS landscape | < 1 min |
| Fig. S.26 | Statistical significance of cluster splits (SigClust2) | 📓 [`08_02_clinical_AD_cluster_analyses.ipynb`](code/notebooks/08_clinical_evaluation/08_02_clinical_AD_cluster_analyses.ipynb) | [`08_02` cell 55](code/notebooks/08_clinical_evaluation/08_02_clinical_AD_cluster_analyses.ipynb#nb-cell-55-output) | [`08_02` cell 52](code/notebooks/08_clinical_evaluation/08_02_clinical_AD_cluster_analyses.ipynb#nb-cell-52), [`08_02` cell 54](code/notebooks/08_clinical_evaluation/08_02_clinical_AD_cluster_analyses.ipynb#nb-cell-54); requires the R environment described in Section 2 | < 1 min |

---

## 5. 🛠️ If something doesn't work

### Environment installation

If `conda env create` reports a `PackagesNotFoundError`, the most likely explanation is that `environment-full.yml` was used on a platform other than Linux, and one of its exact pinned versions has no counterpart there. Use `environment.yml`, which is the cross-platform file.

If `conda env create` reports a `LinkError` mentioning `git-annex` or compiler tools, the conda-provided git-annex is failing to install. Either install `git-annex` separately via your system package manager and remove it from the environment file, or install the `gxx_linux-64` package into your base conda environment to provide the required toolchain.

### Figures look different from the manuscript

If figures render but text appears in a sans-serif font rather than Computer Modern, LaTeX is not being used. Verify your LaTeX installation as described in Section 1.

If figures are missing or contain only axes with no data, the most likely cause is a missing intermediate file. Check the Dependencies column for the relevant figure in Section 4 and run the upstream notebooks first.

### Where to get help

For issues specific to this repository, open an issue at [github.com/sina-mansour/normative_brain_charts/issues](https://github.com/sina-mansour/normative_brain_charts/issues).

For questions about the underlying SNM framework, see the [spectranorm documentation](https://sina-mansour.github.io/spectranorm/).

For correspondence about the manuscript itself, contact the corresponding authors as listed in the published paper.

[![Contact](https://img.shields.io/badge/Contact-sina.mansour.lakouraj@gmail.com-blue)](https://sina-mansour.github.io/)