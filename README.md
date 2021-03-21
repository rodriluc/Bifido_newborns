# Born immune data analyses 
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

### Figure 1 
- ```cell_corr.R``` uses cell abundance dataframe built from Grid that is sub-setted by grouped days for correlation
- Option to re-order one matrix in accordance to the other for comparison purposes
