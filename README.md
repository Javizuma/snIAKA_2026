# snIAKA_2026

This repository contains code used for the analysis of single-nucleus RNA-seq data from the mouse intra-amygdala kainic acid (IAKA) model of temporal lobe epilepsy.

The repository contains scripts for preprocessing, cell-type annotation, differential expression analysis (single-cell and pseudobulk), pathway enrichment analysis, cell–cell communication analysis using CellChat, and generation of figures included in the manuscript.

For full details of the experimental design and analysis, please refer to the Methods section of the associated publication.

## Publication

Villegas Salmerón J. et al. (2026).  
**Single-nucleus transcriptomic analysis reveals subfield-specific cell-to-cell synaptic reorganisation of the mouse hippocampus in focal temporal lobe epilepsy.**  
*Journal of Translational Medicine.*

https://doi.org/10.1186/s12967-026-08978-2

Processed analysis outputs used in the manuscript, including differential expression and CellChat results, are available in the supplementary files associated with the publication.

## Workflow

The preprocessing and analysis scripts are intended to be run in the following order:

**QC → Doublet removal → Filtering → Cell-type annotation → Differential expression → CellChat**

Some downstream scripts use objects generated in previous steps.

Expected doublet rates used during preprocessing were based on the corresponding 10x Genomics documentation.


Cell-type annotation includes an external [MapMyCells](https://knowledge.brain-map.org/mapmycells/process) step. The annotation script generates an `.h5ad` file for upload to MapMyCells, and the resulting annotation file is then imported back into R before continuing with downstream analyses. Smaller MapMyCells subclasses were grouped into broader parent classes to obtain the cell subtype categories used for downstream analyses.

Pathway enrichment analysis was performed using Ingenuity Pathway Analysis (IPA; QIAGEN). As this analysis was carried out within the IPA software, no standalone R script is included in this repository. The pathway enrichment outputs used in the manuscript are provided in the supplementary files associated with the publication.

## Sample naming

Sample identifiers in the GEO/raw data correspond to the following analysis labels:

- **KA1–KA4** = **Epi1–Epi4**
- **PBS1–PBS4** = **Control1–Control4**

These refer to the same biological samples; only the naming convention differs between the raw data and some downstream analysis scripts.

## Software environment

Analyses were performed in R. The main packages and versions used are listed in package_versions.txt. R packages required for each analysis are listed at the beginning of the corresponding scripts. Packages can be installed from CRAN, Bioconductor, or their respective GitHub repositories, as appropriate. Please refer to the package documentation for installation instructions and version requirements.

Full R session information is available in sessionInfo.txt.

## Data availability

Single-nucleus RNA-seq data are available through GEO under accession **GSE343089**.

Users interested in independently interrogating the dataset can download the data from GEO and apply their own analysis workflow. The scripts in this repository document the workflow used for the analyses presented in the publication.

## Contact

For questions about the dataset, analysis, code, or science, feel free to get in touch!

**Javier Villegas Salmerón**  
javiervillegass22@rcsi.ie
