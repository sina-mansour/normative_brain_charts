# Spectral Normative Modeling (SNM)

![Spectral Normative Modeling](assets/banner.webp)

This repository hosts the scripts and code used for our manuscript on **[Spectral Normative Modeling of Brain Structure](https://doi.org/10.1101/2025.01.16.25320639)**. The study focuses on normative modeling of brain development using spectral methods, with applications to both healthy and clinical populations.

**Mansour L, Sina, et al. "Spectral normative modeling of brain structure." medRxiv (2025): 2025-01.**

---

## Overview

This repository contains a collection of Jupyter notebooks that provide a transparent and reproducible workflow for the development and evaluation of Spectral Normative Models (SNMs) as described in the manuscript. The notebooks are organized to reflect the sequential steps taken in the analysis, from data preprocessing to model fitting and clinical evaluation. Notably, the code provided in this repository showcases every step of the analysis, from data preprocessing to model fitting and evaluation, as well as the generation of all figures presented in the original manuscript.

Notably, the code in this repository uses the `Spectranorm` package, which is a dedicated Python package developed to facilitate the implementation of conventional and spectral normative models. For more information on the `Spectranorm` package, please refer to the [Spectranorm Documentation](https://sina-mansour.github.io/spectranorm/).

The following is a brief overview of the analysis pipeline and the corresponding notebooks included in this repository:

### 1. [**Data Import Scripts**](code/notebooks/01_data_import)

Notebooks in this directory are responsible for importing and preprocessing data from several large-scale neuroimaging datasets, comprising over 78,000 healthy brains from 30 different datasets. The notebooks include scripts to:
   - Apply required exclusion criteria to ensure data quality and consistency.
   - Process and clean cortical thickness data and transform cortical thickness maps to the fs_LR 32k space.
   - Prepare data for downstream analysis.

### 2. [**Data Aggregation Scripts**](code/notebooks/02_data_aggregation)

Notebooks in this directory are responsible for aggregating and preparing the cleaned data for subsequent normative modeling. The notebooks include scripts to:
   - Combine cleaned data from multiple datasets.
   - Apply a final round of quality control.
   - Randomly split data into training and test sets for model validation.
   - Generate the spectral coefficients of thickness from the vertex-wise data.

### 3. [**Direct Normative Models**](code/notebooks/03_direct_model)
The [notebook](code/notebooks/03_direct_model/03_01_direct_thickness_model.ipynb) in this directory specifies a direct normative modeling architecture that serves as a backbone for the subsequent spectral normative model. The notebook includes scripts to:
   - Implement a conventional normative modeling architecture using **Hierarchical Bayesian Regression**.
   - Test the normative model for a single variable (e.g., mean cortical thickness) as a function of covariates (age, sex, site). [meant as a sanity check]
   - This architecture is used in subsequent notebooks for both direct and spectral implementations.

### 4. [**Spectral Basis Construction**](code/notebooks/04_eigenmode_basis)
The [notebook](code/notebooks/04_eigenmode_basis/04_01_cortical_connectome_eigenmodes.ipynb) in this directory generates the spectral basis functions (connectome eigenmodes) used to encode high-resolution cortical phenotypes in the spectral normative model.

### 5. [**Spectral Normative Model Fitting**](code/notebooks/05_spectral_model)
The notebook in this directory fits the Spectral Normative Model (SNM) using the spectral coefficients of cortical thickness (derived in [2.3](code/notebooks/02_data_aggregation/02_03_encode_thickness_data.ipynb)). The notebook includes scripts to:
   - Fit an example of the **Spectral Normative Model (SNM)**.
   - Verify the generation of normative ranges (see notebook `05_01_05`).
   - Save the fitted model for subsequent evaluation and analysis.

### 6. [**Performance Evaluations**](code/notebooks/06_performance_evaluation)
The notebooks in this directory evaluate the performance of the fitted SNM and compare it to the direct normative model across a range of spatial queries (i.e., varying granularities of the cortical phenotype) and numbers of modes. The notebooks include scripts to:
   - Evaluate the accuracy of SNMs with varying numbers of modes.
   - Compare SNM performance to a direct model.
   - Generate multiple figures presented in the manuscript.

### 7. [**Cortical Growth Gradients**](code/notebooks/07_growth_gradients)
The [notebook](code/notebooks/07_growth_gradients/07_01_growth_gradients.ipynb) in this directory utilizes the fitted SNM to investigate cortical growth gradients across the lifespan. The notebook includes scripts to:
   - Map three distinct cortical growth gradients across the lifespan.
   - Compare the identified growth gradients to known cortical organization hierarchies.
   - Characterize the lifespan trajectories of change for the identified growth gradients.
   - Generate multiple figures presented in the manuscript.

### 8. [**Clinical (AD) Evaluations**](code/notebooks/08_clinical_evaluation)
The [notebook](code/notebooks/08_clinical_evaluation/08_01_clinical_application_in_AD.ipynb) in this directory applies the fitted SNM to a clinical sample of Alzheimer's Disease (AD) patients. The notebook includes scripts to:
   - Load data from a clinical Alzheimer's Disease (AD) sample (MACC Harmonization dataset).
   - Fine-tune the healthy SNM to the new site by learning harmonization parameters (normative adaptation).
   - Map individual deviations from the healthy normative model in the AD sample.
   - Utilize the identified deviation maps to characterize the normative signature of AD-related cortical atrophy.
   - Generate multiple figures reported in the manuscript.

### 9. [**Data Sharing**](code/notebooks/09_data_sharing)
The notebooks in this directory are responsible for sharing the data used in the manuscript, including a pretrained version of the SNM for cortical thickness, as well as charts generated for multiple alternative cortical parcelations.

---

## References

> Mansour L., S., et al. (2025). Spectral Normative Modeling of Brain Structure. medRxiv. DOI: 10.1101/2025.01.16.25320639

---

## Contact:

If you have any questions or need assistance, please don't hesitate to reach out.

[![sina \[dot\] mansour \[dot\] lakouraj \[at\] gmail](https://img.shields.io/badge/Contact-sina%20[dot]%20mansour%20[dot]%20lakouraj%20[at]%20gmail-blue)](https://sina-mansour.github.io/)
