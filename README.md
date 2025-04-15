# Complementary modes of resistance to EGFR TKI in lung adenocarcinoma through MAPK activation and cellular plasticity

This repository contains code related to the manuscript:

> Zatzman, M., Quintanal-Villalonga, A. _et al_, Complementary modes of resistance to EGFR TKI in lung adenocarcinoma through MAPK activation and cellular plasticity

## Overview

This repository is structued as an installable R package. The `R` folder contains functions related to data processing, visualization, and other helper functions.   The `code` folder contains notebooks related to analysis steps. Explanations of these notebooks can be found within the `code` folder `README.md` file. Vignettes detailing figure generation can be found [here](https://shahcompbio.github.io/egfr_nucseq), or at the following links:

* [Histotime Analysis](https://shahcompbio.github.io/egfr_nucseq/articles/23_histotime.html): Contains code related to histotime analysis in Figure 5
* [IMPACT Analysis](https://shahcompbio.github.io/egfr_nucseq/articles/31_impact_analysis.html): Contains code related to MSK-IMPACT bulk DNA panel sequencing
* [Lung Cancer Cell Atlas Analysis](https://shahcompbio.github.io/egfr_nucseq/articles/41_luca_analysis.html): Contains code related to the lung cancer cell atlas dataset analysis
* [TCGA/CPTAC Analysis](https://shahcompbio.github.io/egfr_nucseq/articles/42_tcga_cptac_analysis.html): Contains code related to bulk RNA LUAD datasets from TCGA and CPTAC
* [Main Plotting Code](https://shahcompbio.github.io/egfr_nucseq/articles/99_plots.html): Contains most of the remaining plotting code for the manuscript

## Data Availablility

Raw sequencing data generated in this project has been deposited to the GEO Expression Omnibus under accession [GSE292469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292469). 

Cohort, epithelial, and tumor histotime scRNA data objects can be interactively visualized and downloaded from [CELLxGENE](https://cellxgene.cziscience.com/collections/df8ef04c-7a87-49d4-a127-f8d653d89d91).

cBioPortal study to browse mutation data can be found [here (link to come)](https://www.cbioportal.org/).
