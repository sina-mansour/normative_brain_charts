# Data Import

The training sample comprises 78,405 scans drawn from 30 datasets across 189
imaging sites. Each notebook below handles one data source: locating and
downloading the imaging data where possible, applying the study's exclusion
criteria, processing cortical thickness maps, and transforming them into
`fs_LR 32k` surface space for downstream aggregation.

Datasets differ substantially in their access conditions. Some are openly
available and are fetched programmatically within the notebook; others require
approval from a data access committee or an executed data transfer agreement,
in which case the notebook documents the processing applied to data obtained
through those channels. No raw imaging data is redistributed through this
repository.

| Notebook | Source |
| --- | --- |
| [01_01](code/notebooks/01_data_import/01_01_hcp_young_adult.ipynb) | HCP Young Adult |
| [01_02](code/notebooks/01_data_import/01_02_hcp_development.ipynb) | HCP Development |
| [01_03](code/notebooks/01_data_import/01_03_hcp_aging.ipynb) | HCP Aging |
| [01_04](code/notebooks/01_data_import/01_04_abcd.ipynb) | ABCD |
| [01_05](code/notebooks/01_data_import/01_05_ukb.ipynb) | UK Biobank |
| [01_06](code/notebooks/01_data_import/01_06_rbc.ipynb) | Reproducible Brain Charts (multiple constituent studies) |
| [01_07](code/notebooks/01_data_import/01_07_BGSP.ipynb) | BGSP |
| [01_08](code/notebooks/01_data_import/01_08_DLBS.ipynb) | DLBS |
| [01_09](code/notebooks/01_data_import/01_09_adhd200.ipynb) | ADHD-200 |
| [01_10](code/notebooks/01_data_import/01_10_aomic.ipynb) | AOMIC |
| [01_11](code/notebooks/01_data_import/01_11_sald.ipynb) | SALD |
| [01_12](code/notebooks/01_data_import/01_12_ixi.ipynb) | IXI |
| [01_13](code/notebooks/01_data_import/01_13_asrb.ipynb) | ASRB |
| [01_14](code/notebooks/01_data_import/01_14_wayne.ipynb) | Wayne State |
| [01_15](code/notebooks/01_data_import/01_15_adni.ipynb) | ADNI |
| [01_16](code/notebooks/01_data_import/01_16_oasis.ipynb) | OASIS |
| [01_17](code/notebooks/01_data_import/01_17_aibl.ipynb) | AIBL |
| [01_18](code/notebooks/01_data_import/01_18_singer.ipynb) | SINGER |
| [01_19](code/notebooks/01_data_import/01_19_sg70.ipynb) | SG70 |
| [01_20](code/notebooks/01_data_import/01_20_tcp.ipynb) | TCP |
| [01_21](code/notebooks/01_data_import/01_21_gusto.ipynb) | GUSTO |
| [01_22](code/notebooks/01_data_import/01_22_sglifespan.ipynb) | SG Lifespan |
| [01_23](code/notebooks/01_data_import/01_23_spresto.ipynb) | S-PRESTO |
| [01_24](code/notebooks/01_data_import/01_24_abide.ipynb) | ABIDE |
| [01_25](code/notebooks/01_data_import/01_25_corr.ipynb) | CoRR |
| [01_26](code/notebooks/01_data_import/01_26_slim.ipynb) | SLIM |
| [01_27](code/notebooks/01_data_import/01_27_macc.ipynb) | MACC (clinical cohort, used in Sections 2.5–2.6) |

