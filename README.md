# Spectral Normative Modeling (SNM)

This repository hosts the scripts and code used for our manuscript on **[Spectral Normative Modeling of Brain Structure](https://doi.org/10.1101/2025.01.16.25320639)**. The study focuses on normative modeling of brain development using spectral methods, with applications to both healthy and clinical populations.

![Spectral Normative Modeling](assets/banner.png)

---

## Overview

The repository is organized into **7 sets of notebooks**, each corresponding to a specific step in the analytical pipeline. Below is a summary of what each step accomplishes:

### 1. **Data Import Scripts** (`code/notebooks/01_data_import`)
   - Process and clean cortical thickness data from healthy participants in the HCP lifespan datasets.
   - Prepare data for downstream analysis.

### 2. **Data Aggregation Scripts** (`code/notebooks/02_data_aggregation`)
   - Combine cleaned data from multiple datasets.
   - Randomly split data into training and test sets for model validation.

### 3. **Normative Models** (`code/notebooks/03_normative_models`)
   - Implement a basic normative model architecture using **Hierarchical Bayesian Regression** (via the PyMC Python package).
   - Normative model for a single variable (e.g., mean cortical thickness) as a function of covariates (age, sex, site).
   - This architecture is used for both direct and spectral implementations.

### 4. **Spectral Basis Construction** (`code/notebooks/04_normative_kernels`)
   - Construct spectral kernels/basis functions to encode high-resolution cortical phenotypes.
   - Map connectome-based brain eigenmodes via singular value decomposition of a random walk graph Laplacian shift operator (see notebook `04_04_02`).

### 5. **Spectral Normative Model Fitting** (`code/05_kernel_normative_models`)
   - Fit a prototype of the **Spectral Normative Model (SNM)**.
   - Verify the generation of normative ranges (see notebook `05_01_05`).

### 6. **Performance Evaluation Scripts** (`code/06_performance_evaluation`)
   - Evaluate the accuracy of SNMs with varying numbers of modes.
   - Compare SNM performance to a direct model.
   - Generate comparative figures presented in the manuscript.

### 7. **Clinical (AD) Evaluations** (`code/07_clinical_evaluation`)
   - Clean and preprocess data from a clinical Alzheimer's Disease (AD) sample (MACC Harmonization dataset).
   - Fine-tune the healthy SNM to the new site by learning harmonization parameters.
   - Perform clinical evaluations and reproduce figures reported in the manuscript.

---

## Citation

If you use this repository in your research, please cite our manuscript:

> Mansour L., S., et al. (2025). Spectral Normative Modeling of Brain Structure. medRxiv. DOI: 10.1101/2025.01.16.25320639

---

## Contact:

If you have any questions or need assistance, please don't hesitate to reach out.

[![sina \[dot\] mansour \[dot\] lakouraj \[at\] gmail](https://img.shields.io/badge/Contact-sina%20[dot]%20mansour%20[dot]%20lakouraj%20[at]%20gmail-blue)](https://sina-mansour.github.io/)
