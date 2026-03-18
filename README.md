# Spectral Normative Modeling of Brain Structure

![Spectral Normative Modeling](assets/banner.webp)

<p align="center">
  <a href="https://doi.org/10.1101/2025.01.16.25320639"><img src="https://img.shields.io/badge/Preprint-medRxiv-red?style=flat-square" alt="medRxiv preprint"></a>
  <a href="https://sina-mansour.github.io/spectranorm/"><img src="https://img.shields.io/badge/Package-SpectraNorm-blue?style=flat-square" alt="SpectraNorm"></a>
  <a href="https://pypi.python.org/pypi/spectranorm/"><img src="https://img.shields.io/pypi/v/spectranorm?style=flat-square" alt="PyPI"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-AGPLv3%20%2F%20Commercial-green?style=flat-square" alt="License"></a>
</p>

This repository provides the complete, reproducible analysis code accompanying the manuscript:

> **Mansour L., S., et al. (2025). Spectral Normative Modeling of Brain Structure. *medRxiv*.** DOI: [10.1101/2025.01.16.25320639](https://doi.org/10.1101/2025.01.16.25320639)

The study introduces **Spectral Normative Models (SNMs)** — a framework for normative modeling of high-dimensional cortical phenotypes using spectral (eigenmode-based) representations. The repository contains every step of the analysis pipeline, from raw data preprocessing through model fitting, performance evaluation, and clinical application, enabling full reproduction of all manuscript figures and results.

---

## 🧰 SpectraNorm Package

The analysis code in this repository is built on top of **SpectraNorm** — a dedicated Python package we developed to facilitate the implementation of conventional and spectral normative models.

[![SpectraNorm Docs](https://img.shields.io/badge/Docs-sina--mansour.github.io/spectranorm-blue?style=flat-square)](https://sina-mansour.github.io/spectranorm/)
[![PyPI](https://img.shields.io/pypi/v/spectranorm?style=flat-square)](https://pypi.python.org/pypi/spectranorm/)

```bash
pip install spectranorm --upgrade
```

For tutorials, API reference, and further documentation, visit the [SpectraNorm documentation site](https://sina-mansour.github.io/spectranorm/).

---

## 📁 Repository Structure

The analysis is organized into **9 sequential stages**, each corresponding to a directory of Jupyter notebooks:

### 1. [Data Import Scripts](code/notebooks/01_data_import)

Imports and preprocesses cortical thickness data from several large-scale neuroimaging datasets, comprising over **78,000 healthy brains** from **30 different datasets**. Includes:
- Applying exclusion criteria to ensure data quality and consistency
- Processing cortical thickness maps and transforming them to **fs_LR 32k** surface space
- Preparing data for downstream analysis

### 2. [Data Aggregation Scripts](code/notebooks/02_data_aggregation)

Aggregates and prepares the cleaned data for normative modeling. Includes:
- Combining cleaned data from multiple datasets
- Final round of quality control
- Random train/test split for model validation
- Generating spectral coefficients of cortical thickness from vertex-wise data

### 3. [Direct Normative Models](code/notebooks/03_direct_model)

Specifies a direct normative modeling architecture as a backbone for the spectral model. Includes:
- Implementing a **Hierarchical Bayesian Regression** normative model (via [PyMC](https://www.pymc.io/))
- Fitting the model for a single variable (mean cortical thickness) as a function of age, sex, and site — used as a sanity check
- This architecture is reused in subsequent notebooks for both direct and spectral implementations

### 4. [Spectral Basis Construction](code/notebooks/04_eigenmode_basis)

Generates the spectral basis functions (connectome eigenmodes) used to encode high-resolution cortical phenotypes. Includes:
- Computing connectome-based brain eigenmodes via SVD of a random walk graph Laplacian shift operator

### 5. [Spectral Normative Model Fitting](code/notebooks/05_spectral_model)

Fits the full Spectral Normative Model (SNM) using the spectral coefficients of cortical thickness. Includes:
- Fitting the **Spectral Normative Model (SNM)**
- Verifying the generation of normative ranges
- Saving the fitted model for subsequent evaluation

### 6. [Performance Evaluations](code/notebooks/06_performance_evaluation)

Evaluates SNM performance and benchmarks it against the direct normative model. Includes:
- Evaluating accuracy across varying numbers of modes and spatial granularities
- Comparing SNM and direct model performance
- Reproducing comparative figures from the manuscript

### 7. [Cortical Growth Gradients](code/notebooks/07_growth_gradients)

Uses the fitted SNM to investigate cortical growth gradients across the lifespan. Includes:
- Mapping three distinct cortical growth gradients
- Comparing growth gradients to known cortical organization hierarchies
- Characterizing lifespan trajectories of change

### 8. [Clinical (AD) Evaluations](code/notebooks/08_clinical_evaluation)

Applies the fitted SNM to a clinical Alzheimer's Disease (AD) sample. Includes:
- Loading data from the MACC Harmonization dataset
- Fine-tuning the healthy SNM via normative adaptation (site harmonization)
- Mapping individual deviations from the normative model
- Characterizing the normative signature of AD-related cortical atrophy

### 9. [Data Sharing](code/notebooks/09_data_sharing)

Shares the data produced in this manuscript. Includes:
- A pretrained SNM for cortical thickness
- Normative brain charts for multiple cortical parcellations

---

## 📄 License

This repository is **dual licensed**:

- **Non-commercial / Academic use**: [GNU AGPLv3](LICENSE)
- **Commercial use**: A separate commercial license is required

See the [LICENSE](LICENSE) file for full details.

---

## 📬 Contact

If you have questions, encounter issues, or would like to collaborate, feel free to reach out:

[![Contact](https://img.shields.io/badge/Contact-sina.mansour.lakouraj@gmail.com-blue)](https://sina-mansour.github.io/)