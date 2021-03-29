###--- Load packages
library(Seurat)
library(tidyverse)
library(data.table)
library(MAST)
library(circlize)
library(cowplot)
library(doMC)
library(scamp)
library(reticulate)
library(RColorBrewer)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MAST")
BiocManager::install("SingleCellExperiment")
BiocManager::install("zinbwave")
BiocManager::install("BiocParallel")
BiocManager::install("scamp")
BiocManager::install("Seurat")


##-- Set working directory
setwd('~/Documents/BabyBifido/re-analysis_BD_Rhaps/All_runsC1234/merge_files')

##-- Data Import
Abseq.cartridge_1 <- fread(input = 'Combined_P20403_1001_DBEC_MolsPerCell_withSampleTag.csv') 
Abseq.cartridge_2 <- fread(input = 'Combined_P20403_1002_DBEC_MolsPerCell_withSampleTag.csv') 
Abseq.cartridge_3 <- fread(input = 'Combined_P20403_1003-1006_DBEC_MolsPerCell_withSampleTag.csv') 
Abseq.cartridge_4 <- fread(input = 'Combined_P20403_1004-1005_DBEC_MolsPerCell_withSampleTag.csv') 

###--- Data processing
#- RNA matrices
Abseq.cartridge.RNA_1 <- Abseq.cartridge_1[, str_detect(string = colnames(Abseq.cartridge_1), pattern = 'pAbO|Cell_Index|Sample_Tag|Sample_Name', negate=TRUE), with = FALSE] # 272 genes
Abseq.cartridge.RNA_1[, Cell_ID := Abseq.cartridge_1$Cell_Index]
Abseq.cartridge.RNA_1[, Sample_Name := Abseq.cartridge_1$Sample_Name]
Abseq.cartridge.RNA_1[, Sample_Tag := Abseq.cartridge_1$Sample_Tag]
Abseq.cartridge.RNA_1[, Batch := 1]

Abseq.cartridge.RNA_2 <- Abseq.cartridge_2[, str_detect(string = colnames(Abseq.cartridge_2), pattern = 'pAbO|Cell_Index|Sample_Tag|Sample_Name', negate=TRUE), with = FALSE] # 272 genes
Abseq.cartridge.RNA_2[, Cell_ID := Abseq.cartridge_2$Cell_Index]
Abseq.cartridge.RNA_2[, Sample_Name := Abseq.cartridge_2$Sample_Name]
Abseq.cartridge.RNA_2[, Sample_Tag := Abseq.cartridge_2$Sample_Tag]
Abseq.cartridge.RNA_2[, Batch := 2]

Abseq.cartridge.RNA_3 <- Abseq.cartridge_3[, str_detect(string = colnames(Abseq.cartridge_3), pattern = 'pAbO|Cell_Index|Sample_Tag|Sample_Name', negate=TRUE), with = FALSE] # 272 genes
Abseq.cartridge.RNA_3[, Cell_ID := Abseq.cartridge_3$Cell_Index]
Abseq.cartridge.RNA_3[, Sample_Name := Abseq.cartridge_3$Sample_Name]
Abseq.cartridge.RNA_3[, Sample_Tag := Abseq.cartridge_3$Sample_Tag]
Abseq.cartridge.RNA_3[, Batch := 3]

Abseq.cartridge.RNA_4 <- Abseq.cartridge_4[, str_detect(string = colnames(Abseq.cartridge_4), pattern = 'pAbO|Cell_Index|Sample_Tag|Sample_Name', negate=TRUE), with = FALSE] # 272 genes
Abseq.cartridge.RNA_4[, Cell_ID := Abseq.cartridge_4$Cell_Index]
Abseq.cartridge.RNA_4[, Sample_Name := Abseq.cartridge_4$Sample_Name]
Abseq.cartridge.RNA_4[, Sample_Tag := Abseq.cartridge_4$Sample_Tag]
Abseq.cartridge.RNA_4[, Batch := 4]

#- Ab matrices
Abseq.cartridge.Ab_1 <- Abseq.cartridge_1[, str_detect(string = colnames(Abseq.cartridge_1), pattern = 'pAbO'), with = FALSE] # 10 markers
Abseq.cartridge.Ab_1[, Cell_ID := Abseq.cartridge_1$Cell_Index]
Abseq.cartridge.Ab_1[, Sample_Name := Abseq.cartridge_1$Sample_Name]
Abseq.cartridge.Ab_1[, Sample_Tag := Abseq.cartridge_1$Sample_Tag]
Abseq.cartridge.Ab_1[, Batch := 1]

Abseq.cartridge.Ab_2 <- Abseq.cartridge_2[, str_detect(string = colnames(Abseq.cartridge_2), pattern = 'pAbO'), with = FALSE] # 10 markers
Abseq.cartridge.Ab_2[, Cell_ID := Abseq.cartridge_2$Cell_Index]
Abseq.cartridge.Ab_2[, Sample_Name := Abseq.cartridge_2$Sample_Name]
Abseq.cartridge.Ab_2[, Sample_Tag := Abseq.cartridge_2$Sample_Tag]
Abseq.cartridge.Ab_2[, Batch := 2]

Abseq.cartridge.Ab_3 <- Abseq.cartridge_3[, str_detect(string = colnames(Abseq.cartridge_3), pattern = 'pAbO'), with = FALSE] # 10 markers
Abseq.cartridge.Ab_3[, Cell_ID := Abseq.cartridge_3$Cell_Index]
Abseq.cartridge.Ab_3[, Sample_Name := Abseq.cartridge_3$Sample_Name]
Abseq.cartridge.Ab_3[, Sample_Tag := Abseq.cartridge_3$Sample_Tag]
Abseq.cartridge.Ab_3[, Batch := 3]

Abseq.cartridge.Ab_4 <- Abseq.cartridge_4[, str_detect(string = colnames(Abseq.cartridge_4), pattern = 'pAbO'), with = FALSE] # 10 markers
Abseq.cartridge.Ab_4[, Cell_ID := Abseq.cartridge_4$Cell_Index]
Abseq.cartridge.Ab_4[, Sample_Name := Abseq.cartridge_4$Sample_Name]
Abseq.cartridge.Ab_4[, Sample_Tag := Abseq.cartridge_4$Sample_Tag]
Abseq.cartridge.Ab_4[, Batch := 4]

#- Features - RNA
identical(names(Abseq.cartridge.RNA_1), names(Abseq.cartridge.RNA_2)) # TRUE
RNA_names <- str_replace(string = names(Abseq.cartridge.RNA_2)[1:259], pattern = '\\|[^(PolyA)]*',  replacement = '_')
RNA_names <- str_replace(string = RNA_names, pattern = '_?$',  replacement = '.rna')
setnames(x = Abseq.cartridge.RNA_1, old = names(Abseq.cartridge.RNA_1), new = c(RNA_names, 'Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch'))
setnames(x = Abseq.cartridge.RNA_2, old = names(Abseq.cartridge.RNA_2), new = c(RNA_names, 'Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch'))
setnames(x = Abseq.cartridge.RNA_3, old = names(Abseq.cartridge.RNA_3), new = c(RNA_names, 'Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch'))
setnames(x = Abseq.cartridge.RNA_4, old = names(Abseq.cartridge.RNA_4), new = c(RNA_names, 'Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch'))

#- Features - Ab
identical(names(Abseq.cartridge.Ab_1), names(Abseq.cartridge.Ab_2)) # TRUE
Ab_names <- sapply(X = str_split(string = names(Abseq.cartridge.Ab_2)[1:10], pattern = '\\|'), 
                   FUN = function(x) paste(x[1], x[2], sep = '|'))
Ab_names <- str_replace(string = Ab_names, pattern = '[|]', replacement = '.')
Ab_names <- paste(Ab_names, 'ab', sep = '.')
setnames(x = Abseq.cartridge.Ab_1, old = names(Abseq.cartridge.Ab_1), new = c(Ab_names, 'Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch'))
setnames(x = Abseq.cartridge.Ab_2, old = names(Abseq.cartridge.Ab_2), new = c(Ab_names, 'Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch'))
setnames(x = Abseq.cartridge.Ab_3, old = names(Abseq.cartridge.Ab_3), new = c(Ab_names, 'Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch'))
setnames(x = Abseq.cartridge.Ab_4, old = names(Abseq.cartridge.Ab_4), new = c(Ab_names, 'Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch'))

###--- Seurat v3
#- Setup a Seurat object
Abseq.RNA <- bind_rows(Abseq.cartridge.RNA_1, Abseq.cartridge.RNA_2, Abseq.cartridge.RNA_3, Abseq.cartridge.RNA_4)
Abseq.RNA.rawData <- t(Abseq.RNA[, -c('Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch')])
colnames(Abseq.RNA.rawData) <- paste('cell', 1:dim(Abseq.RNA.rawData)[2], sep = '-')
PBMC.RNA2 <- Seurat::CreateSeuratObject(counts = Abseq.RNA.rawData, min.cells = 3, min.features = 50, assay = 'RNA', project = 'Abseq bifido 2021')
PBMC.RNA2 # 255 genes across 27,857 cells 

#- Add features - meta.data
PBMC.RNA2@meta.data$cell.id <- rownames(PBMC.RNA2@meta.data)
head(Abseq.RNA)
Abseq.RNA2 <- t(Abseq.RNA)
colnames(Abseq.RNA2) <- paste('cell', 1:dim(Abseq.RNA2)[2], sep = '-')
Abseq.RNA2 <- as.data.frame(t(Abseq.RNA2)) 

Abseq.RNA3 <- subset(Abseq.RNA2, rownames(Abseq.RNA2) %in% rownames(as.data.frame(PBMC.RNA2@meta.data)))
PBMC.RNA2@meta.data$sample.tag <- Abseq.RNA3$Sample_Tag
PBMC.RNA2@meta.data$sample.name <- Abseq.RNA3$Sample_Name
PBMC.RNA2@meta.data$batch <- Abseq.RNA3$Batch

head(x = PBMC.RNA2[[]])

PBMC.RNA2@meta.data$cytokine_culture <- plyr::mapvalues(x = PBMC.RNA2@meta.data$sample.tag, 
                                                        from = unique(PBMC.RNA2@meta.data$sample.tag), 
                                                        to = c("Multiplet","Th1","Th0","Th2", 
                                                               "Th0","Th17","Th1","Th2",               
                                                               "Th17","Th0","Undetermined","Th1",    
                                                               "Th1","Th1","Th0","Th2",  
                                                               "Th17","Th2","Th17","Th0",    
                                                               "Th2","Th1","Th2","Th17" ,  
                                                               "Th17","Th0","Th2","Th0",           
                                                               "Th1","Th2","Th0","Th17" ,          
                                                               "Th1","Th0","Th2","Th17" ))
PBMC.RNA2@meta.data$groups <- plyr::mapvalues(x = PBMC.RNA2@meta.data$sample.tag, 
                                              from = unique(PBMC.RNA2@meta.data$sample.tag), 
                                              to = c("Multiplet","BIFIDOhigh","BIFIDOlow","BIFIDOhigh", 
                                                     "BIFIDOhigh","BIFIDOhigh","CTRL","CTRL",               
                                                     "CTRL","CTRL","Undetermined","BBreveHMO",    
                                                     "BIFIDOlow","BInfantisHMO","BInfantisHMO","BIFIDOlow",  
                                                     "BIFIDOlow","BInfantisHMO","BInfantisHMO","BBreveHMO",    
                                                     "BBreveHMO","pooledHMO","pooledHMO","BBreveHMO" ,  
                                                     "pooledHMO","pooledHMO","10mM_ILA","10mM_ILA",           
                                                     "1mM_ILA","1mM_ILA","1mM_ILA","1mM_ILA" ,          
                                                     "0.1mM_ILA","0.1mM_ILA","0.1mM_ILA","0.1mM_ILA"))
PBMC.RNA2@meta.data$Usamples <- plyr::mapvalues(x = PBMC.RNA2@meta.data$sample.tag, 
                                              from = unique(PBMC.RNA2@meta.data$sample.tag), 
                                              to = c("Multiplet","FW_BIFIDOhigh1:100_Th1","FW_BIFIDOlow1:100_Th0","FW_BIFIDOhigh1:100_Th2", 
                                                     "FW_BIFIDOhigh1:100_Th0","FW_BIFIDOhigh1:100_Th17","CTRL_Th1","CTRL_Th2",               
                                                     "CTRL_Th17","CTRL_Th0","Undetermined","BBreveHMO_1:500_Th1",    
                                                     "FW_BIFIDOlow1:100_Th1","BInfantisHMO_1:500_Th1","BInfantisHMO_1:500_Th0","FW_BIFIDOlow1:100_Th2",  
                                                     "FW_BIFIDOlow1:100_Th17","BInfantisHMO_1:500_Th2","BInfantisHMO_1:500_Th17","BBreveHMO_1:500_Th0",    
                                                     "BBreveHMO_1:500_Th2","pooledHMO_1:500_Th1","pooledHMO_1:500_Th2","BBreveHMO_1:500_Th17" ,  
                                                     "pooledHMO_1:500_Th17","pooledHMO_1:500_Th0","10mM_ILA_Th2","10mM_ILA_Th0",           
                                                     "1mM_ILA_Th1","1mM_ILA_Th2","1mM_ILA_Th0","1mM_ILA_Th17" ,          
                                                     "0.1mM_ILA_Th1","0.1mM_ILA_Th0","0.1mM_ILA_Th2","0.1mM_ILA_Th17"))
ggplot(PBMC.RNA2@meta.data, aes(cytokine_culture))+geom_bar(stat="count")+geom_text(stat="count",aes(label =..count.., vjust = -0.2))

#- Add Ab expression to Seurat object
Abseq.Ab <- bind_rows(Abseq.cartridge.Ab_1, Abseq.cartridge.Ab_2,Abseq.cartridge.Ab_3,Abseq.cartridge.Ab_4)
Abseq.Ab.rawData <- t(Abseq.Ab[, -c('Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch')])
colnames(Abseq.Ab.rawData) <- paste('cell', 1:dim(Abseq.Ab.rawData)[2], sep = '-')
Abseq.Ab.rawData <- Abseq.Ab.rawData[, colnames(PBMC.RNA2)]
PBMC.RNA2[['Ab']] <- CreateAssayObject(counts = Abseq.Ab.rawData)

#- Subsetting and QC
PBMC.RNA2 <- SetIdent(object = PBMC.RNA2, value = 'groups')
#PBMC.RNA2 <- subset(x = PBMC.RNA2, idents = c('BIFIDOlow', 'BIFIDOhigh'))
PBMC.RNA2 <- subset(x = PBMC.RNA2, idents = c('CTRL','BIFIDOlow', 'BIFIDOhigh','BBreveHMO',    
                                                     "BInfantisHMO","pooledHMO","10mM_ILA","1mM_ILA","0.1mM_ILA"))

FeatureScatter(object = PBMC.RNA2, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')#, cex.use = 0.5) 
# Visualize QC metrics as a violin plot
VlnPlot(PBMC.RNA2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
PBMC.RNA2 <- subset(x = PBMC.RNA2, cells = setdiff(colnames(PBMC.RNA2), rownames(PBMC.RNA2@meta.data)[which(PBMC.RNA2@meta.data$nCount_RNA > 2000 & PBMC.RNA2@meta.data$nFeature_RNA < 100)]))
PBMC.RNA2 <- subset(x = PBMC.RNA2, cells = setdiff(colnames(PBMC.RNA2), rownames(PBMC.RNA2@meta.data)[which(PBMC.RNA2@meta.data$nCount_RNA > 4000)])) # 255 genes across 10,003 cells 

#- RNA Normalization
PBMC.RNA2 <- NormalizeData(object = PBMC.RNA2, assay = 'RNA', normalization.method = 'LogNormalize', scale.factor = 10000)

#- Ab Normalization
PBMC.RNA2 <- NormalizeData(object = PBMC.RNA2, assay = 'Ab', normalization.method = 'CLR')


#- Highly variable genes
PBMC.RNA2 <- FindVariableFeatures(object = PBMC.RNA2, verbose = FALSE, selection.method = "vst")
PBMC.RNA2@assays$RNA@var.features # 255 genes

####### Identify the 20 most highly variable genes
top20 <- head(x = VariableFeatures(object = PBMC.RNA2), 
              n =20)
# Plot variable features with labels
plot1 <- VariableFeaturePlot(object = PBMC.RNA2)
LabelPoints(plot = plot1, 
            points = top20, 
            repel = TRUE)

#- Remove uninteresting sources of variation (nUMI as proxy for CDR)
sum(GetAssayData(object = PBMC.RNA2, slot = "data", assay='RNA')['GAPDH.rna',]>0)
all_genes <- rownames(x = PBMC.RNA2)

#GAPDH
percent.GAPDH<-grep(pattern = "GAPDH.rna", x=rownames(x=PBMC.RNA2),value = TRUE)
percent.GAPDHgene <-Matrix::colSums(PBMC.RNA2[percent.GAPDH])/Matrix::colSums(PBMC.RNA2)
PBMC.RNA2<-AddMetaData(object=PBMC.RNA2,metadata = percent.GAPDHgene,col.name = "percent.GAPDHgene")
VlnPlot(PBMC.RNA2, features = 'percent.GAPDHgene')
PBMC.RNA2 <- RunPCA(PBMC.RNA2, features = VariableFeatures(PBMC.RNA2))
PBMC.RNA2 <- RunPCA(PBMC.RNA2, features = percent.GAPDHgene)
DimPlot(PBMC.RNA2)
head(x = PBMC.RNA2[[]])

PBMC.RNA2 <- ScaleData(object = PBMC.RNA2, assay = 'RNA', features= all_genes, 
                       vars.to.regress = c('nCount_RNA','batch','percent.GAPDHgene'))

sum(GetAssayData(object = PBMC.RNA2, slot = "data", assay='RNA')['GAPDH.rna',]>0)
sum(GetAssayData(object = PBMC.RNA2, slot = "data", assay='RNA')['HLA-C.rna',]>0)

ggplot(PBMC.RNA2@meta.data, aes(cytokine_culture))+geom_bar(stat="count")

data_save <- t(GetAssayData(object = PBMC.RNA2, assay = 'RNA'))
data_save_meta <- PBMC.RNA2@meta.data
#write.csv(data_save, file= "RNA_C123456.csv")
#write.csv(data_save_meta, file= "RNA_C123456_meta.csv")

#- PCA
PBMC.RNA2 <- RunPCA(object = PBMC.RNA2, assay = 'RNA')
ElbowPlot(object = PBMC.RNA2, ndims = 50) # 30 PCs
nPC <- 30

#- Clustering
PBMC.RNA2 <- FindNeighbors(object = PBMC.RNA2, reduction = 'pca', dims = 1:nPC, k.param = 30, force.recalc = T)
PBMC.RNA2 <- FindClusters(object = PBMC.RNA2, dims.use = 1:nPC, verbose = TRUE, n.start = 100)

#- Visualization (UMAP)
PBMC.RNA2 <- RunUMAP(object = PBMC.RNA2, reduction = 'pca', dims = 1:nPC, min_dist = 0.2, seed.use = 42, n_neighbors = 30, metric = 'correlation')
DimPlot(object = PBMC.RNA2, reduction = 'umap', group.by = 'cytokine_culture', pt.size = 0.5)
DimPlot(object = PBMC.RNA2, reduction = 'umap', group.by = 'groups', pt.size = 0.5) 

#DE
immune.combined <- PBMC.RNA2

Idents(object = immune.combined)
colnames(x = immune.combined[[]])
Idents(object = immune.combined) <- 'cytokine_culture'
levels(x = immune.combined)

theme_set(theme_cowplot())
cyto_cult_Th0 <- subset(immune.combined, idents = 'Th0') 
Idents(cyto_cult_Th0) <- "groups"
avg.Th0.cells <- log1p(AverageExpression(cyto_cult_Th0, verbose = FALSE)$RNA)
#avg.Th0.cells$gene <- rownames(avg.Th0.cells)
avg.Th0.cells <- dplyr::as_data_frame(avg.Th0.cells, rownames = "gene")

cyto_cult_Th1 <- subset(immune.combined, idents = 'Th1') 
Idents(cyto_cult_Th1) <- "groups"
avg.Th1.cells <- log1p(AverageExpression(cyto_cult_Th1, verbose = FALSE)$RNA)
avg.Th1.cells <-dplyr::as_data_frame(avg.Th1.cells, rownames = "gene")

cyto_cult_Th2 <- subset(immune.combined, idents = 'Th2') 
Idents(cyto_cult_Th2) <- "groups"
avg.Th2.cells <- log1p(AverageExpression(cyto_cult_Th2, verbose = FALSE)$RNA)
avg.Th2.cells <-dplyr::as_data_frame(avg.Th2.cells, rownames = "gene")

cyto_cult_Th17 <- subset(immune.combined, idents = 'Th17') 
Idents(cyto_cult_Th17) <- "groups"
avg.Th17.cells <- log1p(AverageExpression(cyto_cult_Th17, verbose = FALSE)$RNA)
avg.Th17.cells <-dplyr::as_data_frame(avg.Th17.cells, rownames = "gene")

ggplot(immune.combined@meta.data, aes(groups, fill=cytokine_culture))+geom_bar(stat="count") + geom_text(stat="count",aes(label =..count..),size = 5, position = position_stack(vjust = 0.5))#geom_text(stat="count",aes(label =..count..))
#genes.to.label = c("LGALS1.rna")

p1 <- ggplot(avg.Th0.cells, aes(`0.1mM_ILA`, CTRL)) + geom_point()+ geom_point(data=filter(avg.Th0.cells,gene=='CD247.rna'| gene=='TYMS.rna'|gene=='LGALS1.rna'|gene=='PYCR1.rna'|gene=='SPOCK2.rna'),colour='red') + ggtitle("Th0") + geom_text_repel(data=filter(avg.Th0.cells,gene=='CD247.rna'| gene=='TYMS.rna'|gene=='LGALS1.rna'|gene=='PYCR1.rna'|gene=='SPOCK2.rna'), aes(label=gene)) #+ geom_text(aes(label=gene))
p2 <- ggplot(avg.Th2.cells, aes(`0.1mM_ILA`, CTRL)) + geom_point() + geom_point(data=filter(avg.Th2.cells,gene=='CD4.rna'| gene=='PYCR1.rna'|gene=='SPOCK2.rna'|gene=='ZAP70.rna'|gene=='LGALS1.rna'), colour='red')+ ggtitle("Th2") + geom_text_repel(data=filter(avg.Th2.cells,gene=='CD4.rna'| gene=='PYCR1.rna'|gene=='SPOCK2.rna'|gene=='ZAP70.rna'|gene=='LGALS1.rna'), aes(label=gene))
p3 <- ggplot(avg.Th17.cells, aes(`0.1mM_ILA`, CTRL)) + geom_point() + geom_point(data=filter(avg.Th17.cells,gene=='LGALS1.rna'| gene=='CD4.rna'|gene=='SPOCK2.rna'), colour='red') + ggtitle("Th17") + geom_text_repel(data=filter(avg.Th17.cells,gene=='LGALS1.rna'| gene=='CD4.rna'|gene=='SPOCK2.rna'), aes(label=gene))
p4 <- ggplot(avg.Th1.cells, aes(`0.1mM_ILA`, CTRL)) + geom_point() + geom_point(data=filter(avg.Th1.cells,gene=='TYMS.rna'| gene=='PYCR1.rna'|gene=='LGALS1.rna'|gene=='LTA.rna'), colour='red') + ggtitle("Th1") + geom_text_repel(data=filter(avg.Th1.cells,gene=='TYMS.rna'| gene=='PYCR1.rna'|gene=='LGALS1.rna'|gene=='LTA.rna'), aes(label=gene))
plot_grid(p1, p2, p3, p4)

#top30
p1 <- ggplot(avg.Th0.cells, aes(`1mM_ILA`, CTRL)) + geom_point()+ geom_point(data=filter(avg.Th0.cells,gene=='TK1.rna'| gene=='PTTG2.rna'|gene=='CXCR3.rna'|gene=='UBE2C.rna'|gene=='AURKB.rna'|gene=='TOP2A.rna'|gene=='HMGB2.rna'|gene=='NCR3.rna'|gene=='MKI67.rna'|gene=='CD70.rna'|gene=='HMMR.rna'|gene=='CD9.rna'|gene=='CCNB1.rna'|gene=='TYMS.rna'|gene=='IL12RB2.rna'|gene=='LGALS1.rna'|gene=='CXCR6.rna'|gene=='PRF1.rna'|gene=='CHI3L2.rna'|gene=='TNFSF10.rna'|gene=='KIT.rna'|gene=='IL9R.rna'|gene=='TCF7.rna'|gene=='SPOCK2.rna'|gene=='LCK.rna'|gene=='NKG7.rna'|gene=='ICAM1.rna'|gene=='APOBEC3G.rna'|gene=='OAS1.rna'|gene=='FOXP3.rna'),colour='red') + ggtitle("Th0") + geom_text_repel(data=filter(avg.Th0.cells,gene=='TK1.rna'| gene=='PTTG2.rna'|gene=='CXCR3.rna'|gene=='UBE2C.rna'|gene=='AURKB.rna'|gene=='TOP2A.rna'|gene=='HMGB2.rna'|gene=='NCR3.rna'|gene=='MKI67.rna'|gene=='CD70.rna'|gene=='HMMR.rna'|gene=='CD9.rna'|gene=='CCNB1.rna'|gene=='TYMS.rna'|gene=='IL12RB2.rna'|gene=='LGALS1.rna'|gene=='CXCR6.rna'|gene=='PRF1.rna'|gene=='CHI3L2.rna'|gene=='TNFSF10.rna'|gene=='KIT.rna'|gene=='IL9R.rna'|gene=='TCF7.rna'|gene=='SPOCK2.rna'|gene=='LCK.rna'|gene=='NKG7.rna'|gene=='ICAM1.rna'|gene=='APOBEC3G.rna'|gene=='OAS1.rna'|gene=='FOXP3.rna'), aes(label=gene)) #+ geom_text(aes(label=gene))
p2 <- ggplot(avg.Th2.cells, aes(`1mM_ILA`, CTRL)) + geom_point() + geom_point(data=filter(avg.Th2.cells,gene=='TK1.rna'| gene=='PTTG2.rna'|gene=='AURKB.rna'|gene=='TOP2A.rna'|gene=='UBE2C.rna'|gene=='LGALS1.rna'|gene=='NKG7.rna'|gene=='TYMS.rna'|gene=='HMGB2.rna'|gene=='NCR3.rna'|gene=='MKI67.rna '|gene=='CCNB1.rna'|gene=='CXCR3.rna'|gene=='TNFSF10.rna'|gene=='HMMR.rna'|gene=='PRF1.rna'|gene=='IL12RB2.rna'|gene=='CD70.rna'|gene=='IL9R.rna'|gene=='ICAM1.rna'|gene=='CHI3L2.rna'|gene=='CD9.rna'|gene=='OAS1.rna'|gene=='IL13.rna'|gene=='GZMB.rna'|gene=='LIF.rna'|gene=='LAP3.rna'|gene=='CTSW.rna'|gene=='TCF7.rna'|gene=='CSF2.rna'), colour='red')+ ggtitle("Th2") + geom_text_repel(data=filter(avg.Th2.cells,gene=='TK1.rna'| gene=='PTTG2.rna'|gene=='AURKB.rna'|gene=='TOP2A.rna'|gene=='UBE2C.rna'|gene=='LGALS1.rna'|gene=='NKG7.rna'|gene=='TYMS.rna'|gene=='HMGB2.rna'|gene=='NCR3.rna'|gene=='MKI67.rna '|gene=='CCNB1.rna'|gene=='CXCR3.rna'|gene=='TNFSF10.rna'|gene=='HMMR.rna'|gene=='PRF1.rna'|gene=='IL12RB2.rna'|gene=='CD70.rna'|gene=='IL9R.rna'|gene=='ICAM1.rna'|gene=='CHI3L2.rna'|gene=='CD9.rna'|gene=='OAS1.rna'|gene=='IL13.rna'|gene=='GZMB.rna'|gene=='LIF.rna'|gene=='LAP3.rna'|gene=='CTSW.rna'|gene=='TCF7.rna'|gene=='CSF2.rna'), aes(label=gene))
p3 <- ggplot(avg.Th17.cells, aes(`1mM_ILA`, CTRL)) + geom_point() + geom_point(data=filter(avg.Th17.cells,gene=='TK1.rna'| gene=='GZMB.rna'|gene=='CD70.rna'|gene=='PTTG2.rna'|gene=='UBE2C.rna'|gene=='TYMS.rna'|gene=='AURKB.rna'|gene=='CD9.rna'|gene=='TOP2A.rna'|gene=='CXCR3.rna'|gene=='MKI67.rna'|gene=='IL12RB2.rna '|gene=='JUNB.rna'|gene=='NKG7.rna'|gene=='HMMR.rna'|gene=='HMGB2.rna'|gene=='LGALS1.rna'|gene=='PRF1.rna'|gene=='CCL5.rna'|gene=='CCNB1.rna'|gene=='CHI3L2.rna'|gene=='ZBED2.rna'|gene=='IL9R.rna'|gene=='NCR3.rna'|gene=='CCR4.rna'|gene=='IL2RA.rna'|gene=='LIF.rna'|gene=='IL13.rna'|gene=='HAVCR2.rna'|gene=='TNFRSF8.rna'), colour='red') + ggtitle("Th17") + geom_text_repel(data=filter(avg.Th17.cells,gene=='TK1.rna'| gene=='GZMB.rna'|gene=='CD70.rna'|gene=='PTTG2.rna'|gene=='UBE2C.rna'|gene=='TYMS.rna'|gene=='AURKB.rna'|gene=='CD9.rna'|gene=='TOP2A.rna'|gene=='CXCR3.rna'|gene=='MKI67.rna'|gene=='IL12RB2.rna '|gene=='JUNB.rna'|gene=='NKG7.rna'|gene=='HMMR.rna'|gene=='HMGB2.rna'|gene=='LGALS1.rna'|gene=='PRF1.rna'|gene=='CCL5.rna'|gene=='CCNB1.rna'|gene=='CHI3L2.rna'|gene=='ZBED2.rna'|gene=='IL9R.rna'|gene=='NCR3.rna'|gene=='CCR4.rna'|gene=='IL2RA.rna'|gene=='LIF.rna'|gene=='IL13.rna'|gene=='HAVCR2.rna'|gene=='TNFRSF8.rna'), aes(label=gene))
p4 <- ggplot(avg.Th1.cells, aes(`1mM_ILA`, CTRL)) + geom_point() + geom_point(data=filter(avg.Th1.cells,gene=='PTTG2.rna'| gene=='IL12RB2.rna'|gene=='TK1.rna'|gene=='AURKB.rna'|gene=='CCNB1.rna'|gene=='UBE2C.rna'|gene=='CD9.rna'|gene=='TNFRSF8.rna'|gene=='TOP2A.rna'|gene=='GZMB.rna'|gene=='PIK3IP1.rna'|gene=='CSF2.rna'|gene=='HMMR.rna'|gene=='HAVCR2.rna'|gene=='HMGB2.rna'|gene=='CCL5.rna'|gene=='OAS1.rna'|gene=='JUNB.rna'|gene=='TYMS.rna'|gene=='ZBED2.rna'|gene=='MKI67.rna'|gene=='TNFSF10.rna'|gene=='CXCR3.rna'|gene=='NKG7.rna'|gene=='ICAM1.rna'|gene=='LGALS1.rna'|gene=='IFNG.rna'|gene=='CD70.rna'|gene=='LAP3.rna'|gene=='CXCR5.rna'), colour='red') + ggtitle("Th1") + geom_text_repel(data=filter(avg.Th1.cells,gene=='PTTG2.rna'| gene=='IL12RB2.rna'|gene=='TK1.rna'|gene=='AURKB.rna'|gene=='CCNB1.rna'|gene=='UBE2C.rna'|gene=='CD9.rna'|gene=='TNFRSF8.rna'|gene=='TOP2A.rna'|gene=='GZMB.rna'|gene=='PIK3IP1.rna'|gene=='CSF2.rna'|gene=='HMMR.rna'|gene=='HAVCR2.rna'|gene=='HMGB2.rna'|gene=='CCL5.rna'|gene=='OAS1.rna'|gene=='JUNB.rna'|gene=='TYMS.rna'|gene=='ZBED2.rna'|gene=='MKI67.rna'|gene=='TNFSF10.rna'|gene=='CXCR3.rna'|gene=='NKG7.rna'|gene=='ICAM1.rna'|gene=='LGALS1.rna'|gene=='IFNG.rna'|gene=='CD70.rna'|gene=='LAP3.rna'|gene=='CXCR5.rna'), aes(label=gene))
plot_grid(p1, p2, p3, p4)

#top20
p1 <- ggplot(avg.Th0.cells, aes(`1mM_ILA`, CTRL)) + geom_point()+ geom_point(data=filter(avg.Th0.cells,gene=='TK1.rna'| gene=='PTTG2.rna'|gene=='CXCR3.rna'|gene=='UBE2C.rna'|gene=='AURKB.rna'|gene=='TOP2A.rna'|gene=='HMGB2.rna'|gene=='NCR3.rna'|gene=='MKI67.rna'|gene=='CD70.rna'|gene=='HMMR.rna'|gene=='CD9.rna'|gene=='CCNB1.rna'|gene=='TYMS.rna'|gene=='IL12RB2.rna'|gene=='LGALS1.rna'|gene=='CXCR6.rna'|gene=='PRF1.rna'|gene=='CHI3L2.rna'|gene=='TNFSF10.rna'),colour='red') + ggtitle("Th0") + geom_text_repel(data=filter(avg.Th0.cells,gene=='TK1.rna'| gene=='PTTG2.rna'|gene=='CXCR3.rna'|gene=='UBE2C.rna'|gene=='AURKB.rna'|gene=='TOP2A.rna'|gene=='HMGB2.rna'|gene=='NCR3.rna'|gene=='MKI67.rna'|gene=='CD70.rna'|gene=='HMMR.rna'|gene=='CD9.rna'|gene=='CCNB1.rna'|gene=='TYMS.rna'|gene=='IL12RB2.rna'|gene=='LGALS1.rna'|gene=='CXCR6.rna'|gene=='PRF1.rna'|gene=='CHI3L2.rna'|gene=='TNFSF10.rna'), aes(label=gene)) #+ geom_text(aes(label=gene))
p2 <- ggplot(avg.Th2.cells, aes(`1mM_ILA`, CTRL)) + geom_point() + geom_point(data=filter(avg.Th2.cells,gene=='TK1.rna'| gene=='PTTG2.rna'|gene=='AURKB.rna'|gene=='TOP2A.rna'|gene=='UBE2C.rna'|gene=='LGALS1.rna'|gene=='NKG7.rna'|gene=='TYMS.rna'|gene=='HMGB2.rna'|gene=='NCR3.rna'|gene=='MKI67.rna '|gene=='CCNB1.rna'|gene=='CXCR3.rna'|gene=='TNFSF10.rna'|gene=='HMMR.rna'|gene=='PRF1.rna'|gene=='IL12RB2.rna'|gene=='CD70.rna'|gene=='IL9R.rna'|gene=='ICAM1.rna'), colour='red')+ ggtitle("Th2") + geom_text_repel(data=filter(avg.Th2.cells,gene=='TK1.rna'| gene=='PTTG2.rna'|gene=='AURKB.rna'|gene=='TOP2A.rna'|gene=='UBE2C.rna'|gene=='LGALS1.rna'|gene=='NKG7.rna'|gene=='TYMS.rna'|gene=='HMGB2.rna'|gene=='NCR3.rna'|gene=='MKI67.rna '|gene=='CCNB1.rna'|gene=='CXCR3.rna'|gene=='TNFSF10.rna'|gene=='HMMR.rna'|gene=='PRF1.rna'|gene=='IL12RB2.rna'|gene=='CD70.rna'|gene=='IL9R.rna'|gene=='ICAM1.rna'), aes(label=gene))
p3 <- ggplot(avg.Th17.cells, aes(`1mM_ILA`, CTRL)) + geom_point() + geom_point(data=filter(avg.Th17.cells,gene=='TK1.rna'| gene=='GZMB.rna'|gene=='CD70.rna'|gene=='PTTG2.rna'|gene=='UBE2C.rna'|gene=='TYMS.rna'|gene=='AURKB.rna'|gene=='CD9.rna'|gene=='TOP2A.rna'|gene=='CXCR3.rna'|gene=='MKI67.rna'|gene=='IL12RB2.rna '|gene=='JUNB.rna'|gene=='NKG7.rna'|gene=='HMMR.rna'|gene=='HMGB2.rna'|gene=='LGALS1.rna'|gene=='PRF1.rna'|gene=='CCL5.rna'|gene=='CCNB1.rna'), colour='red') + ggtitle("Th17") + geom_text_repel(data=filter(avg.Th17.cells,gene=='TK1.rna'| gene=='GZMB.rna'|gene=='CD70.rna'|gene=='PTTG2.rna'|gene=='UBE2C.rna'|gene=='TYMS.rna'|gene=='AURKB.rna'|gene=='CD9.rna'|gene=='TOP2A.rna'|gene=='CXCR3.rna'|gene=='MKI67.rna'|gene=='IL12RB2.rna '|gene=='JUNB.rna'|gene=='NKG7.rna'|gene=='HMMR.rna'|gene=='HMGB2.rna'|gene=='LGALS1.rna'|gene=='PRF1.rna'|gene=='CCL5.rna'|gene=='CCNB1.rna'), aes(label=gene))
p4 <- ggplot(avg.Th1.cells, aes(`1mM_ILA`, CTRL)) + geom_point() + geom_point(data=filter(avg.Th1.cells,gene=='PTTG2.rna'| gene=='IL12RB2.rna'|gene=='TK1.rna'|gene=='AURKB.rna'|gene=='CCNB1.rna'|gene=='UBE2C.rna'|gene=='CD9.rna'|gene=='TNFRSF8.rna'|gene=='TOP2A.rna'|gene=='GZMB.rna'|gene=='PIK3IP1.rna'|gene=='CSF2.rna'|gene=='HMMR.rna'|gene=='HAVCR2.rna'|gene=='HMGB2.rna'|gene=='CCL5.rna'|gene=='OAS1.rna'|gene=='JUNB.rna'|gene=='TYMS.rna'|gene=='ZBED2.rna'), colour='red') + ggtitle("Th1") + geom_text_repel(data=filter(avg.Th1.cells,gene=='PTTG2.rna'| gene=='IL12RB2.rna'|gene=='TK1.rna'|gene=='AURKB.rna'|gene=='CCNB1.rna'|gene=='UBE2C.rna'|gene=='CD9.rna'|gene=='TNFRSF8.rna'|gene=='TOP2A.rna'|gene=='GZMB.rna'|gene=='PIK3IP1.rna'|gene=='CSF2.rna'|gene=='HMMR.rna'|gene=='HAVCR2.rna'|gene=='HMGB2.rna'|gene=='CCL5.rna'|gene=='OAS1.rna'|gene=='JUNB.rna'|gene=='TYMS.rna'|gene=='ZBED2.rna'), aes(label=gene))
plot_grid(p1, p2, p3, p4)

DE_Th1 <- FindMarkers(cyto_cult_Th1, ident.1 = "1mM_ILA", ident.2 = "CTRL", verbose = FALSE, min.diff.pct = 0.1)
head(DE_Th1, n = 30)
DE_Th17 <- FindMarkers(cyto_cult_Th17, ident.1 = "1mM_ILA", ident.2 = "CTRL", verbose = FALSE, min.diff.pct = 0.1) #min.diff.pct = 0.1
head(DE_Th17, n = 30)
DE_Th2 <- FindMarkers(cyto_cult_Th2, ident.1 = "1mM_ILA", ident.2 = "CTRL", verbose = FALSE, min.diff.pct = 0.1)
head(DE_Th2, n = 30)
DE_Th0 <- FindMarkers(cyto_cult_Th0, ident.1 = "1mM_ILA", ident.2 = "CTRL", verbose = FALSE, min.diff.pct = 0.1)
head(DE_Th0, n = 30)
