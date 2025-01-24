# Spectral Normative Modeling (SNM)

---

This repository hosts the associated scripts used for our manuscript on [Spectral Normative Modeling of Brain Structure](https://doi.org/10.1101/2025.01.16.25320639)
A repository hosting codes used for normative models of brain development.

---

The scripts are structured into 7 different sets of notebooks, together containing all required analytical steps utilized in our study.

The following is a description of what notebooks in each step accomplish:

1. **Data Import Scripts**: The scripts from `code/notebooks/01_data_import` were used to process and clean the cortical thickness data of healthy participants from the HCP lifespan datasets.
2. **Data Aggregation Scripts**: The scripts form `code/notebooks/02_data_aggregation` were used to combine the clean data from different datasets, in addition to scripts utilized to randomly split the data into training and test sets.
3. **Normative Models**: The scripts from `code/notebooks/03_normative_models` were used to implement a basic normative model architecture to model a single variable as a function of covariates (age, sex, site). We implemented a Hierarchical Bayesian Regression architecture utilizing the PyMC Python package for probabilistic programming. The same model architecture was used for both the direct and the spectral implementations.
4. **Spectral Basis Construction**: The scripts form `code/notebooks/04_normative_kernels` were used to construct spectral kernels/basis functions to encode the high-resolution phenotypes. Notably, we mapped connectome-based brain eigenmodes via singular value decomposition of a random walk graph Laplacian shift operator (see notebook 04_04_02).
5. **Spectral Normative Model Fitting**: The scripts from `code/05_kernel_normative_models` were used to fit a prototype of the Spectral Normative Model and verify its generation of normative ranges (see notebook 05_01_05).
6. **Performance Evaluation Scripts**: The scripts from `code/06_performance_evaluation` were used to evaluate the accuracy of SNMs with different number of modes and compare that to a direct model. This directory contains all relevant scripts to reproduce the comparative figures presented in our manuscript.
7. **Clinical (AD) Evaluations**: The scripts from `code/07_clinical_evaluation` were used to (i) clean the clinical AD sample's data (MACC Harmonization dataset), (ii) fine-tune the healthy SNM to this new site (learning harmonization parameters), and (iii) perform all clinical evaluations reported in the manuscript (including reproducible figures).

# Spectral Normative Modeling (SNM)

This repository hosts the scripts and code used for our manuscript on **[Spectral Normative Modeling of Brain Structure](https://doi.org/10.1101/2025.01.16.25320639)**. The study focuses on normative modeling of brain development using spectral methods, with applications to both healthy and clinical populations.

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

> Mansour L., et al. (2025). Spectral Normative Modeling of Brain Structure. medRxiv. DOI: 10.1101/2025.01.16.25320639

---

## Contact:

If you have any questions or need assistance, please don't hesitate to reach out.

[![sina \[dot\] mansour \[dot\] lakouraj \[at\] gmail](https://img.shields.io/badge/Contact-sina%20[dot]%20mansour%20[dot]%20lakouraj%20[at]%20gmail-blue)](https://sina-mansour.github.io/)
