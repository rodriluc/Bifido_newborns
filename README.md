# IMPRINT & Born immune data analyses 
## Bifidobacteria-mediated immune system imprinting early in life (Henrick et al.) 

## Table of contents
* [General info](#general-info)
* [Dependencies](#dependencies)
* [Repo description](#repo-description)
	* [Figure 1](#figure-1)
	* [Figure 2](#figure-2)
	* [Figure 3](#figure-3)
	* [Figure 4](#figure-4)
	* [Figure 5](#figure-5)
	* [Figure 6](#figure-6)

## General info
This project used multiple omics data:
- Plasma protein expression (Olink - NPX values)
- Cell abundance (CyTOF - Grid cell abundance)
- Metagenome data (Shotgun metagenomic sequencing - CPM)
- Targeted single-cell RNA-seq (BD Rhapsody - MolsPerCell)
	
## Dependencies
Project is created with:
* RStudio version: 4.0.2
* Python version: 3.8

## Repo description

## Figure 2
### Metagenomics
- 2A: PCoA from bray-curtis distance matrix at family level created in R using ... Script available in `figure2/2A.R`.
- 2B: areaplot of family level relative abundances created in python using `pandas`, `numpy`, `matplotlib` and `seaborn`. Guide available in `figure2/2B.html`.
- 2C: ...plot of Bifidobacterium species relative abundance trajectory per individual. Created in R using ... Script available in `figure2/2C.R`. 

## Figure 5
### Cytokines and Metagenomics
- 5C: Radar plot of 0-1 normalized cytokine medians created in R using `tidyverse`, `ggsignif`, `ggiraphExtra`, `ggpubr` and `forcats`. Script available in `figure5/Cytokines_and_Metagenomics.R`.
- 5D: Heatmap of cytokine and bacterial family relative abundances created in R using `tidyverse`, `RColorBrewer`, `ggsignif`, `ggpubr`, `psych` and `forcats`. Guide available in `figure5/Cytokines_and_Metagenomics.R`.

## Figure 6
### Targeted BD Rhapsody data
- 6B & 6D: UMAP of each cytokine condition and volcano plot displaying differentially expresed mRNA as well as processing of output files found in ```BD_CD4Tpolarization_exp1.R```
- 6F: Dot plot of different conditions for all cytokine cultures as well as processing of output files found in ```BD_CD4Tpolarization_exp2.R```
