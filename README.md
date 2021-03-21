# Bifidobacteria-mediated immune system imprinting early in life (Henrick et al.) 
## Born immune data analyses 

## Table of contents
* [General info](#general-info)
* [Dependencies](#dependencies)
* [Repo description](#repo-description)
* [Figures](#figures)

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
- */Figure 1* 
- */Figure 2* 
- */Figure 3* 
- */Figure 4* 
- */Figure 5* 
- */Figure 6*

## Figures
### Spearman Correlation 
- ```cell_corr.R``` uses cell abundance dataframe built from Grid that is sub-setted by grouped days for correlation
- Option to re-order one matrix in accordance to the other for comparison purposes
