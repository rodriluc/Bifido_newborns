###--- Load packages
library(Seurat)
library(tidyverse)
library(data.table)
library(MAST)
library(circlize)
library(cowplot)
library(doMC)
library(reticulate)
library(RColorBrewer)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MAST")
BiocManager::install("SingleCellExperiment")
BiocManager::install("zinbwave")
BiocManager::install("BiocParallel")
BiocManager::install("scamp")

##-- Set working directory
setwd('~/Documents/BabyBifido/BD_Rhapsody_NaiveCD4T/Results')

##-- Data Import
Abseq.cartridge_1 <- fread(input = 'Combined_P17458_1001_DBEC_MolsPerCell_withSampleTag.csv') 
Abseq.cartridge_2 <- fread(input = 'Combined_P17458_1002_DBEC_MolsPerCell_withSampleTag.csv') 
Abseq.cartridge_3 <- fread(input = 'Combined_P17458_1003_DBEC_MolsPerCell_withSampleTag.csv') 
Abseq.cartridge_4 <- fread(input = 'Combined_P17458_1004_DBEC_MolsPerCell_withSampleTag.csv') 

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
identical(names(Abseq.cartridge.RNA_4), names(Abseq.cartridge.RNA_2)) # TRUE
RNA_names <- str_replace(string = names(Abseq.cartridge.RNA_2)[1:259], pattern = '\\|[^(PolyA)]*',  replacement = '_')
RNA_names <- str_replace(string = RNA_names, pattern = '_?$',  replacement = '.rna')
setnames(x = Abseq.cartridge.RNA_1, old = names(Abseq.cartridge.RNA_1), new = c(RNA_names, 'Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch'))
setnames(x = Abseq.cartridge.RNA_2, old = names(Abseq.cartridge.RNA_2), new = c(RNA_names, 'Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch'))
setnames(x = Abseq.cartridge.RNA_3, old = names(Abseq.cartridge.RNA_3), new = c(RNA_names, 'Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch'))
setnames(x = Abseq.cartridge.RNA_4, old = names(Abseq.cartridge.RNA_4), new = c(RNA_names, 'Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch'))
#- Features - Ab
identical(names(Abseq.cartridge.Ab_3), names(Abseq.cartridge.Ab_2)) # TRUE
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
PBMC.RNA2 <- CreateSeuratObject(counts = Abseq.RNA.rawData, min.cells = 3, min.features = 50, assay = 'RNA', project = 'Abseq bifido 2020')
PBMC.RNA2 # 255 genes across 55,421 cells 

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
                                                       to = c("Th1", "Th17", "Multiplet", "Th0", "Th0", "IFNb", "iTreg", "Th2", 
                                                              "Undetermined", "Th0", "Th2", "Th17", "Th1", 
                                                              "Th1", "IFNb", "iTreg", "Th2", "IFNb", 
                                                              "Th17",  "Th17", "Th1", "iTreg", "iTreg", 
                                                              "Th0", "Th2", "Th0", "IFNb", "Th17", 
                                                              "Th1", "Th2","iTreg", "IFNb"))
PBMC.RNA2@meta.data$bifido <- plyr::mapvalues(x = PBMC.RNA2@meta.data$sample.tag, 
                                             from = unique(PBMC.RNA2@meta.data$sample.tag), 
                                             to = c("CTRL", "CTRL", "Multiplet", "CTRL", "BIFIDOpos", "CTRL", "CTRL", "CTRL", 
                                                    "Undetermined", "BIFIDOpos", "BIFIDOpos", "BIFIDOpos", "BIFIDOpos", 
                                                    "BIFIDOpos", "BIFIDOpos", "BIFIDOpos", "BIFIDOpos", "BIFIDOpos", 
                                                    "BIFIDOneg",  "BIFIDOpos", "BIFIDOneg", "BIFIDOpos", "BIFIDOneg", 
                                                    "BIFIDOneg", "BIFIDOneg", "BIFIDOneg", "BIFIDOneg", "BIFIDOneg", 
                                                    "BIFIDOneg", "BIFIDOneg","BIFIDOneg", "BIFIDOneg"))
PBMC.RNA2@meta.data$concentration <- plyr::mapvalues(x = PBMC.RNA2@meta.data$sample.tag, 
                                                    from = unique(PBMC.RNA2@meta.data$sample.tag), 
                                                    to = c("CTRL", "CTRL", "Multiplet", "CTRL", "1:100", "CTRL", "CTRL", "CTRL", 
                                                           "Undetermined", "1:333", "1:333", "1:100", "1:333", 
                                                           "1:100", "1:100", "1:100", "1:100", "1:333", 
                                                           "1:100",  "1:333", "1:100", "1:333", "1:100", 
                                                           "1:100", "1:100", "1:333", "1:333", "1:333", 
                                                           "1:333", "1:333","1:333", "1:100"))
PBMC.RNA2@meta.data$Usamples <- plyr::mapvalues(x = PBMC.RNA2@meta.data$sample.tag, 
                                            from = unique(PBMC.RNA2@meta.data$sample.tag), 
                                            to = c("Th1_CTRL", "Th17_CTRL", "Multiplet", "Th0_CTRL", "Th0_BIFIDOpos1:100", "IFNb_CTRL", "iTreg_CTRL", "Th2_CTRL", 
                                                   "Undetermined", "Th0_BIFIDOpos1:333", "Th2_BIFIDOpos1:333", "Th17_BIFIDOpos1:100", "Th1_BIFIDOpos1:333", 
                                                   "Th1_BIFIDOpos1:100", "IFNb_BIFIDOpos1:100", "iTreg_BIFIDOpos1:100", "Th2_BIFIDOpos1:100", "IFNb_BIFIDOpos1:333", 
                                                   "Th17_BIFIDOneg1:100",  "Th17_BIFIDOpos1:333", "Th1_BIFIDOneg_1:100", "iTreg_BIFIDOpos1:333", "iTreg_BIFIDOneg1:100", 
                                                   "Th0_BIFIDOneg1:100", "Th2_BIFIDOneg1:100", "Th0_BIFIDOneg1:333", "IFNb_BIFIDOneg1:333", "Th17_BIFIDOneg1:333", 
                                                   "Th1_BIFIDOneg1:333", "Th2_BIFIDOneg1:333","iTreg_BIFIDOneg1:333", "IFNb_BIFIDOneg1:100"))
PBMC.RNA2@meta.data$bifido_conc <- plyr::mapvalues(x = PBMC.RNA2@meta.data$sample.tag, 
                                               from = unique(PBMC.RNA2@meta.data$sample.tag), 
                                               to = c("CTRL", "CTRL", "Multiplet", "CTRL", "BIFIDOpos1:100", "CTRL", "CTRL", "CTRL", 
                                                      "Undetermined", "BIFIDOpos1:333", "BIFIDOpos1:333", "BIFIDOpos1:100", "BIFIDOpos1:333", 
                                                      "BIFIDOpos1:100", "BIFIDOpos1:100", "BIFIDOpos1:100", "BIFIDOpos1:100", "BIFIDOpos1:333", 
                                                      "BIFIDOneg1:100",  "BIFIDOpos1:333", "BIFIDOneg1:100", "BIFIDOpos1:333", "BIFIDOneg1:100", 
                                                      "BIFIDOneg1:100", "BIFIDOneg1:100", "BIFIDOneg1:333", "BIFIDOneg1:333", "BIFIDOneg1:333", 
                                                      "BIFIDOneg1:333", "BIFIDOneg1:333","BIFIDOneg1:333", "BIFIDOneg1:100"))

ggplot(PBMC.RNA2@meta.data, aes(cytokine_culture))+geom_bar(stat="count")+geom_text(stat="count",aes(label =..count.., vjust = -0.2))

#- Add Ab expression to Seurat object
Abseq.Ab <- bind_rows(Abseq.cartridge.Ab_1, Abseq.cartridge.Ab_2, Abseq.cartridge.Ab_3, Abseq.cartridge.Ab_4)
Abseq.Ab.rawData <- t(Abseq.Ab[, -c('Sample_Name', 'Sample_Tag', 'Cell_ID', 'Batch')])
colnames(Abseq.Ab.rawData) <- paste('cell', 1:dim(Abseq.Ab.rawData)[2], sep = '-')
Abseq.Ab.rawData <- Abseq.Ab.rawData[, colnames(PBMC.RNA2)]
PBMC.RNA2[['Ab']] <- CreateAssayObject(counts = Abseq.Ab.rawData)

#- Subsetting and QC
PBMC.RNA2 <- SetIdent(object = PBMC.RNA2, value = 'bifido')
PBMC.RNA2 <- subset(x = PBMC.RNA2, idents = c('BIFIDOneg', 'BIFIDOpos'))
FeatureScatter(object = PBMC.RNA2, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')#, cex.use = 0.5) 
# Visualize QC metrics as a violin plot
VlnPlot(PBMC.RNA2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
PBMC.RNA2 <- subset(x = PBMC.RNA2, cells = setdiff(colnames(PBMC.RNA2), rownames(PBMC.RNA2@meta.data)[which(PBMC.RNA2@meta.data$nCount_RNA > 2000 & PBMC.RNA2@meta.data$nFeature_RNA < 100)]))
PBMC.RNA2 <- subset(x = PBMC.RNA2, cells = setdiff(colnames(PBMC.RNA2), rownames(PBMC.RNA2@meta.data)[which(PBMC.RNA2@meta.data$nCount_RNA > 4000)])) # 256 genes across 42,750 cells 

#- RNA Normalization
PBMC.RNA2 <- NormalizeData(object = PBMC.RNA2, assay = 'RNA', normalization.method = 'LogNormalize', scale.factor = 10000)

#- Ab Normalization
PBMC.RNA2 <- NormalizeData(object = PBMC.RNA2, assay = 'Ab', normalization.method = 'CLR')

#- Highly variable genes
PBMC.RNA2 <- FindVariableFeatures(object = PBMC.RNA2, verbose = FALSE, selection.method = "vst")
PBMC.RNA2@assays$RNA@var.features 

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

#- PCA
PBMC.RNA2 <- RunPCA(object = PBMC.RNA2, assay = 'RNA')
ElbowPlot(object = PBMC.RNA2, ndims = 50) 
nPC <- 30

#- Clustering
PBMC.RNA2 <- FindNeighbors(object = PBMC.RNA2, reduction = 'pca', dims = 1:nPC, k.param = 30, force.recalc = T)
PBMC.RNA2 <- FindClusters(object = PBMC.RNA2, dims.use = 1:nPC, verbose = TRUE, n.start = 100)

#- Visualization (UMAP)
PBMC.RNA2 <- RunUMAP(object = PBMC.RNA2, reduction = 'pca', dims = 1:nPC, min_dist = 0.2, seed.use = 42, n_neighbors = 30, metric = 'correlation')
DimPlot(object = PBMC.RNA2, reduction = 'umap', group.by = 'cytokine_culture', pt.size = 0.5)
DimPlot(object = PBMC.RNA2, reduction = 'umap', group.by = 'bifido', pt.size = 0.5) 
DimPlot(object = PBMC.RNA2, reduction = 'umap', group.by = 'concentration', pt.size = 0.5) 
DimPlot(object = PBMC.RNA2, reduction = 'umap', group.by = 'bifido_conc', pt.size = 0.5) 

#DE
immune.combined <- PBMC.RNA2

Idents(object = immune.combined)
colnames(x = immune.combined[[]])
Idents(object = immune.combined) <- 'cytokine_culture'
levels(x = immune.combined)

theme_set(theme_cowplot())
cyto_cult_Th0 <- subset(immune.combined, idents = 'Th0') 
Idents(cyto_cult_Th0) <- "bifido"
avg.Th0.cells <- log1p(AverageExpression(cyto_cult_Th0, verbose = FALSE)$RNA)
avg.Th0.cells$gene <- rownames(avg.Th0.cells)

cyto_cult_Th1 <- subset(immune.combined, idents = 'Th1') 
Idents(cyto_cult_Th1) <- "bifido"
avg.Th1.cells <- log1p(AverageExpression(cyto_cult_Th1, verbose = FALSE)$RNA)
avg.Th1.cells$gene <- rownames(avg.Th1.cells)

cyto_cult_Th2 <- subset(immune.combined, idents = 'Th2') 
Idents(cyto_cult_Th2) <- "bifido"
avg.Th2.cells <- log1p(AverageExpression(cyto_cult_Th2, verbose = FALSE)$RNA)
avg.Th2.cells$gene <- rownames(avg.Th2.cells)

cyto_cult_Th17 <- subset(immune.combined, idents = 'Th17') 
Idents(cyto_cult_Th17) <- "bifido"
avg.Th17.cells <- log1p(AverageExpression(cyto_cult_Th17, verbose = FALSE)$RNA)
avg.Th17.cells$gene <- rownames(avg.Th17.cells)

cyto_cult_IFNb <- subset(immune.combined, idents = 'IFNb') 
Idents(cyto_cult_IFNb) <- "bifido"
avg.IFNb.cells <- log1p(AverageExpression(cyto_cult_IFNb, verbose = FALSE)$RNA)
avg.IFNb.cells$gene <- rownames(avg.IFNb.cells)

cyto_cult_iTreg <- subset(immune.combined, idents = 'iTreg') 
Idents(cyto_cult_iTreg) <- "bifido"
avg.iTreg.cells <- log1p(AverageExpression(cyto_cult_iTreg, verbose = FALSE)$RNA)
avg.iTreg.cells$gene <- rownames(avg.iTreg.cells)

#genes.to.label = c("CXCL13.rna", "CD6.rna", "CD2.rna")
p1 <- ggplot(avg.Th0.cells, aes(BIFIDOpos, BIFIDOneg)) + geom_point() + ggtitle("Th0") #+ geom_text(aes(label=gene))
p2 <- ggplot(avg.Th2.cells, aes(BIFIDOpos, BIFIDOneg)) + geom_point() + ggtitle("Th2")
p3 <- ggplot(avg.Th17.cells, aes(BIFIDOpos, BIFIDOneg)) + geom_point() + ggtitle("Th17")
p4 <- ggplot(avg.Th1.cells, aes(BIFIDOpos, BIFIDOneg)) + geom_point() + ggtitle("Th1")
p5 <- ggplot(avg.IFNb.cells, aes(BIFIDOpos, BIFIDOneg)) + geom_point() + ggtitle("IFNb")
p6 <- ggplot(avg.iTreg.cells, aes(BIFIDOpos, BIFIDOneg)) + geom_point() + ggtitle("iTreg")
#p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
plot_grid(p1, p2, p3, p4, p5, p6)

DE_iTreg <- FindMarkers(cyto_cult_iTreg, ident.1 = "BIFIDOpos", ident.2 = "BIFIDOneg", verbose = FALSE, min.diff.pct = 0.1) #min.pct = 0, logfc.threshold = 0
head(DE_iTreg, n = 15)
FeaturePlot(cyto_cult_iTreg, features = c("GAPDH.rna", "GIMAP7.rna", "CD3G.rna"), split.by = "bifido", max.cutoff = 3, 
            cols = c("grey", "red"))
DE_IFNb <- FindMarkers(cyto_cult_IFNb, ident.1 = "BIFIDOpos", ident.2 = "BIFIDOneg", verbose = FALSE, min.diff.pct = 0.1)
head(DE_IFNb, n = 15)
FeaturePlot(cyto_cult_iTreg, features = c("GAPDH.rna", "IL23R.rna", "HLA-C.rna"), split.by = "concentration", max.cutoff = 3, 
            cols = c("grey", "red"))
DE_Th1 <- FindMarkers(cyto_cult_Th1, ident.1 = "BIFIDOpos", ident.2 = "BIFIDOneg", verbose = FALSE, min.diff.pct = 0.1)
head(DE_Th1, n = 15)
FeaturePlot(cyto_cult_Th1, features = c("GAPDH.rna", "IL12RB2.rna", "HLA-C.rna"), split.by = "bifido", max.cutoff = 3, 
            cols = c("grey", "red"))
DE_Th17 <- FindMarkers(cyto_cult_Th17, ident.1 = "BIFIDOpos", ident.2 = "BIFIDOneg", verbose = FALSE, min.diff.pct = 0.1) #min.diff.pct = 0.1
head(DE_Th17, n = 15)
FeaturePlot(cyto_cult_Th17, features = c("GAPDH.rna", "IL23R.rna", "GIMAP7.rna"), split.by = "bifido", max.cutoff = 3, 
            cols = c("grey", "red"))
DE_Th2 <- FindMarkers(cyto_cult_Th2, ident.1 = "BIFIDOpos", ident.2 = "BIFIDOneg", verbose = FALSE, min.diff.pct = 0.1)
head(DE_Th2, n = 15)
FeaturePlot(cyto_cult_Th2, features = c("GAPDH.rna", "GIMAP7.rna", "HLA-C.rna"), split.by = "bifido", max.cutoff = 3, 
            cols = c("grey", "red"))
DE_Th0 <- FindMarkers(cyto_cult_Th0, ident.1 = "BIFIDOpos", ident.2 = "BIFIDOneg", verbose = FALSE, min.diff.pct = 0.1)
head(DE_Th0, n = 15)
head(FindMarkers(cyto_cult_Th0, ident.1 = "BIFIDOpos", ident.2 = "BIFIDOneg", test.use = "DESeq2", max.cells.per.ident = 50))
FeaturePlot(cyto_cult_Th0, features = c("GAPDH.rna", "GIMAP7.rna", "HLA-C.rna"), split.by = "bifido", max.cutoff = 3, 
            cols = c("grey", "red"))

## Volcano pots DE
devtools::install_github('kevinblighe/EnhancedVolcano')
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

EnhancedVolcano(DE_Th0,
                lab = rownames(DE_Th0),
                x = 'avg_logFC',
                y = 'p_val',
                title = 'Th0: bifido+/-',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
)

library(ggplot2)
library(patchwork)
library(cowplot)
immune.combined <- PBMC.RNA2

Idents(object = immune.combined)
colnames(x = immune.combined[[]])
Idents(object = immune.combined) <- 'cytokine_culture'
levels(x = immune.combined)

cyto_cult_Th0 <- subset(immune.combined, idents = 'Th0') 
Idents(cyto_cult_Th0) <- "bifido"

markers.RNA <- FindMarkers(immune.combined, ident.1='Th0', ident.2 = NULL,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25, # log-FC
                              test.use = 'MAST',
                              verbose = FALSE,
                              assay = 'RNA',
                              latent.vars = 'nCount_RNA') # nUMI as proxy for CDR
top10RNA <- markers.RNA %>% 
  group_by('cytokine_culture') %>%
  top_n(n = 10)

features <- c('STAT1.rna', 'CD48.rna', 'ITGB2.rna', 'NCR3.rna', 'STAT4.rna', 'GAPDH.rna', 'IRF8.rna', 'TNFSF10.rna', 'BCL2.rna', 'CBLB.rna')
features <- c('IL23R.rna', 'IL12RB2.rna')
plots <- VlnPlot(immune.combined, features = features , split.by = "concentration", group.by = "cytokine_culture", 
                 pt.size = 0, combine = FALSE) 
wrap_plots(plots = plots, ncol = 1)
RidgePlot(immune.combined, features = features ,  group.by = "cytokine_culture") 

###--- UMAP --------- STOP
data_ggplot <- data.table(PBMC.RNA2@meta.data, Embeddings(object = PBMC.RNA2, reduction = 'umap'))
plot_UMAP <- ggplot(data_ggplot %>% arrange(sample(x = cell.id, replace = FALSE)), aes(x = UMAP_1, y = UMAP_2)) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

#- Batch
plot_UMAP + geom_point(aes(color = as.character(batch)), size = 0.1) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE)
#There is no batch (cartridges) effect

#- Sequencing Depth
plot_UMAP + geom_point(aes(color = nCount_RNA), size = 0.2) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + 
  scale_color_continuous(guide = FALSE, low = 'royalblue', high = 'red')
#There is no sequencing effect (except for MPS)

#- cytokines
plot_UMAP + geom_point(aes(color = concentration), size = 0.1) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE) +
  facet_wrap(~ cytokine_culture)

#- bifido
plot_UMAP + geom_point(aes(color = cytokine_culture), size = 0.1) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE) +
  facet_wrap(~ bifido)

#- concentrations
plot_UMAP + geom_point(aes(color = concentration), size = 0.1) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE) +
  facet_wrap(~ concentration)

#- samples
plot_UMAP + geom_point(aes(color = Usamples), size = 0.1) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE) +
  facet_wrap(~ Usamples)

#- bifido_conc
plot_UMAP + geom_point(aes(colour = bifido_conc), size = 0.1) + #scale_color_manual(values = cytokine_culture)+
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE) +
  facet_wrap(~ bifido_conc)

##-- Marker genes
DimPlot(object = PBMC.RNA2, reduction = 'umap', group.by = 'RNA_snn_res.0.8', pt.size = 0.5) # 19 clusters
plot_UMAP + geom_point(aes(colour = RNA_snn_res.0.8), size = 0.1) + #scale_color_manual(values = cytokine_culture)+
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE) +
  facet_wrap(~ cytokine_culture)

#- Find markers for every cluster compared to all remaining cells, report only the positive ones
markers.RNA <- FindAllMarkers(object = PBMC.RNA2, 
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25, # log-FC
                              test.use = 'MAST',
                              verbose = FALSE,
                              assay = 'RNA',
                              latent.vars = 'nCount_RNA') # nUMI as proxy for CDR
top10RNA <- markers.RNA %>% 
  group_by(cluster) %>%
  top_n(n = 10)

cluster0.markers <- c('IL12RB2.rna','IL18RAP.rna','IL18R1.rna','TK1.rna','TYMS.rna',
                      'S1PR1.rna','STAT1.rna','PTTG2.rna','IL32.rna','NINJ2.rna')
cluster1.markers <- c('SELPLG.rna','VNN2.rna','TCF7.rna','TRAC.rna','TRIB2.rna',
                      'SPOCK2.rna','TRBC2.rna','STAT6.rna','TIAF1.rna','TXK.rna')
cluster2.markers <- c('SEMA7A.rna','SELL.rna','ZBED2.rna','TXK.rna','STAT6.rna',
                      'SPOCK2.rna','STAT4.rna','STAT3.rna','TRBC2.rna','VNN2.rna')
cluster3.markers <- c('TNFSF10.rna','TYMS.rna','OAS1.rna','NKG7.rna','STAT4.rna',
                      'TNFRSF1B.rna','STAT1.rna','TK1.rna','RUNX3.rna','PYCR1.rna')
cluster4.markers <- c('UBE2C.rna','LGALS3.rna','TOP2A.rna','PTTG2.rna','ZBED2.rna',
                      'TYMS.rna','MKI67.rna','TK1.rna','SEMA7A.rna','SELL.rna')
cluster5.markers <- c('HMMR.rna','TOP2A.rna','UBE2C.rna','MKI67.rna','HMGB2.rna',
                      'PTTG2.rna','TYMS.rna','LEF1.rna','TRIB2.rna','TCF7.rna')
cluster6.markers <- c('TNFRSF8.rna','ZBED2.rna','RORC.rna','TNF.rna','SEMA7A.rna',
                      'SLAMF1.rna','TK1.rna','TYMS.rna','STAT3.rna','TNFRSF25.rna')
cluster7.markers <- c('PRDM1.rna','SELPLG.rna','TNFRSF1B.rna','TNFSF10.rna','TRAT1.rna',
                      'SLAMF1.rna','TNF.rna','PRF1.rna','SELL.rna','TRIB2.rna')
cluster8.markers <- c('TIGIT.rna','TNFRSF1B.rna','PRDM1.rna','TNFRSF8.rna','SELL.rna',
                      'PASK.rna','TNF.rna','PRF1.rna','SLAMF1.rna','PYCR1.rna')
cluster9.markers <- c('PTGDR2.rna','TRIB2.rna','NCR3.rna','SELPLG.rna','VNN2.rna',
                      'TRAC.rna','S1PR1.rna','PIK3IP1.rna','NKG7.rna','PRF1.rna')
cluster10.markers <- c('NAMPT.rna','MYC.rna','TNF.rna','SEMA7A.rna','PYCR1.rna',
                       'SLAMF1.rna','ZBED2.rna','STAT3.rna','STAT5A.rna','NKG7.rna')
cluster11.markers <- c('LGALS3.rna','PIK3IP1.rna','TXK.rna','NKG7.rna','RUNX3.rna',
                       'SEMA7A.rna','LTB.rna','TRAT1.rna','VNN2.rna','ZAP70.rna')
cluster12.markers <- c('NKG7.rna','PRF1.rna','TARP-refseq.rna','LAG3.rna','TBX21.rna',
                       'PRDM1.rna','LTA.rna','TNF.rna','SLAMF1.rna','LAIR2.rna')
cluster13.markers <- c('MKI67.rna','TOP2A.rna','UBE2C.rna','ZAP70.rna','LTB.rna',
                       'TXK.rna','TK1.rna','IFNGR1.rna','VNN2.rna','JUNB.rna')

protein_levels <- c('CD103.ITGAE.ab', 'CD123.IL3RA.ab', 'CD161:DX12.KLRB1.ab', #'CD34:581.CD34.ab',
                   'CD38.CD38.ab', 'CD39.ENTPD1.ab', 'CD45RA:HI100.PTPRC.ab', #'CD99.CD99.ab',
                   'HLA-DR.CD74.ab', 'TCR-gamma-delta:B1.TRD-TRG.ab')

Th1.markers <- c('CXCR3.rna', 'FBXO22.rna', 'HAVCR2.rna')
Th17.markers <- c('IL17.rna', 'RORA.rna', 'RORC.rna')
Th2.markers <- c('GATA3.rna', 'IL25.rna', 'PTGDR2.rna', 'STAT6.rna')
exhausted.markers <- c('CD274.rna', 'FOSB.rna', 'HAVCR2.rna', 'LAG3.rna')
CD4naive.markers <- 'LAT.rna'
CD4.markers <- c('C10orf54.rna', 'CCR10.rna', 'CD4.rna', 'IL9.rna', 'LAT.rna', 'PMCH.rna')
gdT.markers <- c('CD300A.rna', 'TARP.rna', 'TRDC.rna', 'VNN2.rna')

VlnPlot(object = PBMC.RNA2, features = cluster13.markers)#, pt=0)
FeaturePlot(object = PBMC.RNA2, features = gdT.markers, cols = c("grey", "blue"), 
            reduction = "umap")#, min.cutoff = "q05", max.cutoff = "q95")
RidgePlot(object=PBMC.RNA2, features = protein_levels, ncol = 3)

DimPlot(object = PBMC.RNA2, reduction = 'umap', group.by = 'RNA_snn_res.0.8', pt.size = 0.5) # 19 clusters
plot_UMAP + geom_point(aes(colour = RNA_snn_res.0.8), size = 0.1) + #scale_color_manual(values = cytokine_culture)+
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE) +
  facet_wrap(~ cytokine_culture)



###- Ab expression
#- Setup seurat object
PBMC.RNA_Ab <- PBMC.RNA2
DefaultAssay(object = PBMC.RNA_Ab) <- 'Ab'
PBMC.AB <- PBMC.RNA_Ab

all_Ab <- rownames(x = PBMC.AB)

#- Normalization
PBMC.AB <- NormalizeData(object = PBMC.AB, assay = "Ab", normalization.method = 'CLR')

#- Highly variABle genes
PBMC.AB <- FindVariableFeatures(object = PBMC.AB, verbose = FALSE)
PBMC.AB@assays$Ab@var.features # 8 features

#- Scale data
PBMC.AB <- ScaleData(object = PBMC.AB, assay = 'Ab', 
                     features= all_Ab, vars.to.regress = 'batch')

#- PCA
PBMC.AB <- RunPCA(object = PBMC.AB, assay = 'Ab', npcs = 30)
ElbowPlot(object = PBMC.AB, ndims = 7) # 20 PCs
nPC <- 6

#- Clustering
PBMC.AB <- FindNeighbors(object = PBMC.AB, reduction = 'pca', dims = 1:nPC)
PBMC.AB <- FindClusters(object = PBMC.AB, dims.use = 1:nPC, verbose = TRUE)

#PrintFindClustersParams(object = PBMC.AB)

#- Visualization (UMAP)
PBMC.AB <- RunUMAP(object = PBMC.AB, reduction = 'pca', dims = 1:nPC, min_dist = 0.2, seed.use = 42, n_neighbors = 30, metric = 'correlation')
DimPlot(object = PBMC.AB, reduction = 'umap', group.by = 'cytokine_culture', pt.size = 0.5)
DimPlot(object = PBMC.AB, reduction = 'umap', group.by = 'Ab_snn_res.0.8', pt.size = 0.5) # 17 clusters

###--- UMAP
data_ggplot <- data.table(PBMC.AB@meta.data, Embeddings(object = PBMC.AB, reduction = 'umap'))
plot_UMAP <- ggplot(data_ggplot %>% arrange(sample(x = cell.id, replace = FALSE)), aes(x = UMAP_1, y = UMAP_2)) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())

#- Batch
plot_UMAP + geom_point(aes(color = as.character(batch)), size = 0.2) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE)
#There is no batch (cartridges) effect. 

#- Sequencing depth
plot_UMAP + geom_point(aes(color = nCount_Ab), size = 0.2) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + 
  scale_color_continuous(guide = FALSE, low = 'royalblue', high = 'red')
#There is no sequencing effect (except for MPS).  

#- cytokines
plot_UMAP + geom_point(aes(color = cytokine_culture), size = 0.1) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE) +
  facet_wrap(~ cytokine_culture)

#- bifido
plot_UMAP + geom_point(aes(color = bifido), size = 0.1) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE) +
  facet_wrap(~ bifido)

#- concentrations
plot_UMAP + geom_point(aes(color = concentration), size = 0.1) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE) +
  facet_wrap(~ concentration)

##-- Marker genes
#- Find markers for every cluster compared to all remaining cells, report only the positive ones
markers.AB <- FindAllMarkers(object = PBMC.AB, 
                             only.pos = TRUE,
                             min.pct = 0.25,
                             logfc.threshold = 0.25, # log-FC
                             test.use = 'MAST',
                             verbose = FALSE,
                             assay = 'Ab') 
markers.AB %>% 
  group_by(cluster) %>%
  top_n(n = 10)

gdT.markers <- c('CD300A.rna', 'TARP_refseq.rna', 'TRDC.rna', 'VNN2.rna', 'TCR-gamma-delta:B1.TRD-TRG.ab')
Ab.markers <-c('CD39.ENTPD1.ab', 'TCR-gamma-delta:B1.TRD-TRG.ab', 'CD103.ITGAE.ab', 'CD161:DX12.KLRB1.ab', 'CD38.CD38.ab', 'CD123.IL3RA.ab', 'CD45RA:HI100.PTPRC.ab', 'HLA-DR.CD74.ab')
Activated.markers <- c('HLA-DR.CD74.ab', 'CCR8.rna', 'CD69.rna', 'GHR.rna',  
                       'LRRC32.rna', 'PYCR1.rna', 'SEMA7A.rna', 'TNFRSF18.rna',
                       'TNFRSF25.rna', 'TNFRSF8.rna')
Th1.markers <- c('CXCR3.rna', 'FBXO22.rna', 'HAVCR2.rna')
Th17.markers <- c('IL17.rna', 'RORA.rna', 'RORC.rna')
Th2.markers <- c('GATA3.rna', 'IL25.rna', 'PTGDR2.rna', 'STAT6.rna')
exhausted.markers <- c('CD274.rna', 'FOSB.rna', 'HAVCR2.rna', 'LAG3.rna')
CD4naive.markers <- 'LAT.rna'
CD4.markers <- c('C10orf54.rna', 'CCR10.rna', 'CD4.rna', 'IL9.rna', 'LAT.rna', 'PMCH.rna')
naive.markers <- c('SELL.rna', 'CCR7.rna', 'LRRN3.rna')
helperT.markers <- c('CCL20.rna', 'IL13.rna', 'IL17F.rna', 'IL2.rna', 'IL21.rna', 'IL23R.rna',
                     'IL3.rna', 'IL4.rna', 'IL4R.rna', 'IL5.rna', 'IL6.rna', 'LIF.rna',
                     'SELL.rna', 'STAT3.rna', 'TNF.rna', 'ZBED2.rna')

RidgePlot(object = PBMC.AB, features = gdT.markers, assay = 'Ab')
FeaturePlot(object = PBMC.AB, #min.cutoff = "q05", max.cutoff = "q95", 
            features = gdT.markers)
VlnPlot(object = PBMC.AB, features = gdT.markers, pt=0)

########################################
PBMC.AB@meta.data$new.cluster.ids <- plyr::mapvalues(x = PBMC.AB@meta.data$Ab_snn_res.0.8, 
                                                   from = unique(PBMC.AB@meta.data$Ab_snn_res.0.8), 
                                                   to = c('1','3','2','gdT cells','gdT cells','7','9','12','5','8','17','4','10','11', '6','15',
                                                          '18','14','19','0','13','16','23','22'))

DimPlot(object = PBMC.AB, reduction = 'umap', group.by = 'new.cluster.ids', pt.size = 0.5) 

PBMC.AB@meta.data$cc_bifido <- plyr::mapvalues(x = PBMC.AB@meta.data$sample.tag, 
                                                from = unique(PBMC.AB@meta.data$sample.tag), 
                                                to = c('Th0_BIFIDOpos',   'Th2_BIFIDOpos',   'Th0_BIFIDOpos',   'Th17_BIFIDOpos', 
                                                       'Th1_BIFIDOpos',   'Th1_BIFIDOpos',   'IFNb_BIFIDOpos',  'iTreg_BIFIDOpos',
                                                       'Th2_BIFIDOpos',   'Th17_BIFIDOneg',  'Th1_BIFIDOneg',  'IFNb_BIFIDOpos', 
                                                       'Th17_BIFIDOpos',  'iTreg_BIFIDOpos', 'iTreg_BIFIDOneg', 'Th0_BIFIDOneg', 
                                                       'Th2_BIFIDOneg',   'IFNb_BIFIDOneg',  'Th17_BIFIDOneg',  'Th1_BIFIDOneg', 
                                                       'Th0_BIFIDOneg',   'Th2_BIFIDOneg',   'iTreg_BIFIDOneg', 'IFNb_BIFIDOneg'))


data_ggplot <- data.table(PBMC.AB@meta.data, Embeddings(object = PBMC.AB, reduction = 'umap'))
plot_UMAP <- ggplot(data_ggplot %>% arrange(sample(x = cell.id, replace = FALSE)), aes(x = UMAP_1, y = UMAP_2)) + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())
plot_UMAP + geom_point(aes(color = cc_bifido), size = 0.1) +
  labs(x = 'UMAP 1', y = 'UMAP 2')  + scale_color_discrete(guide = FALSE) +
  facet_wrap(~ cc_bifido)

# stacked barplot
devtools::install_github("dtm2451/dittoSeq")
library(dittoSeq)

dittoBarPlot(
  object = PBMC.AB,
  var = "Ab_snn_res.0.8",
  group.by = "cc_bifido",
  scale = 'percent',
  main = 'Relative abundance of clusters (12 conditions)')
  #main = 'Relative abundance of clusters (cytokine cultures)',
  #xlab = 'Clusters')

########################################
# Remove gdT cells
Idents(PBMC.AB) <- "new.cluster.ids"
PBMC.AB_gdTremoved <- subset(x = PBMC.AB, idents = 'gdT cells', invert = TRUE) 
Idents(PBMC.AB_gdTremoved) <- "concentration"
PBMC.AB_gdTremoved_conc <- subset(x = PBMC.AB_gdTremoved, idents = '1:100') 

head(x = PBMC.RNA_gdTremoved_conc100[[]])
head(x = PBMC.RNA2[[]])
#write.csv(PBMC.AB@meta.data, file='PBMC.AB-CLUSTER.csv')
#write.csv(PBMC.RNA2@meta.data, file='PBMC.RNA-AbCLUSTER.csv')
PBMC.RNA2@meta.data$Ab.cluster.ids <- PBMC.AB@meta.data$new.cluster.ids
Idents(PBMC.RNA2) <- "Ab.cluster.ids"
PBMC.RNA_gdTremoved <- subset(x = PBMC.RNA2, idents = 'gdT cells', invert = TRUE) 

Idents(PBMC.RNA_gdTremoved) <- "concentration"
PBMC.RNA_gdTremoved_conc100 <- subset(x = PBMC.RNA_gdTremoved, idents = '1:100') 

#Idents(PBMC.RNA_gdTremoved) <- "concentration"
#PBMC.RNA_gdTremoved_conc100 <- subset(x = PBMC.RNA_gdTremoved, idents = '1:100') 

library(data.table)
data_save <- t(GetAssayData(object = PBMC.RNA_gdTremoved_conc100, assay = 'RNA'))
data_save_meta <- PBMC.RNA_gdTremoved_conc100@meta.data
#write.csv(data_save, file= "RNA_conc1:100.csv")
#write.csv(data_save_meta, file= "RNA_meta1:100.csv")

#- Visualization (UMAP)
DimPlot(object = PBMC.AB_gdTremoved, reduction = 'umap', group.by = 'Ab_snn_res.0.8', pt.size = 0.5) 

# DE Analysis
immune.combinedAb <- PBMC.RNA_gdTremoved_conc100 #conc

Idents(object = immune.combinedAb)
colnames(x = immune.combinedAb[[]])
Idents(object = immune.combinedAb) <- 'cytokine_culture'
levels(x = immune.combinedAb)

theme_set(theme_cowplot())
cyto_cult_Th0 <- subset(immune.combinedAb, idents = 'Th0') 
Idents(cyto_cult_Th0) <- "bifido"
avg.Th0.cells <- log1p(AverageExpression(cyto_cult_Th0, verbose = FALSE)$RNA)
avg.Th0.cells$gene <- rownames(avg.Th0.cells)

cyto_cult_Th1 <- subset(immune.combinedAb, idents = 'Th1') 
Idents(cyto_cult_Th1) <- "bifido"
avg.Th1.cells <- log1p(AverageExpression(cyto_cult_Th1, verbose = FALSE)$RNA)
avg.Th1.cells$gene <- rownames(avg.Th1.cells)

cyto_cult_Th2 <- subset(immune.combinedAb, idents = 'Th2') 
Idents(cyto_cult_Th2) <- "bifido"
avg.Th2.cells <- log1p(AverageExpression(cyto_cult_Th2, verbose = FALSE)$RNA)
avg.Th2.cells$gene <- rownames(avg.Th2.cells)

cyto_cult_Th17 <- subset(immune.combinedAb, idents = 'Th17') 
Idents(cyto_cult_Th17) <- "bifido"
avg.Th17.cells <- log1p(AverageExpression(cyto_cult_Th17, verbose = FALSE)$RNA)
avg.Th17.cells$gene <- rownames(avg.Th17.cells)

cyto_cult_IFNb <- subset(immune.combinedAb, idents = 'IFNb') 
Idents(cyto_cult_IFNb) <- "bifido"
avg.IFNb.cells <- log1p(AverageExpression(cyto_cult_IFNb, verbose = FALSE)$RNA)
avg.IFNb.cells$gene <- rownames(avg.IFNb.cells)

cyto_cult_iTreg <- subset(immune.combinedAb, idents = 'iTreg') 
Idents(cyto_cult_iTreg) <- "bifido"
avg.iTreg.cells <- log1p(AverageExpression(cyto_cult_iTreg, verbose = FALSE)$RNA)
avg.iTreg.cells$gene <- rownames(avg.iTreg.cells)

genes.to.label = c("IL23R.rna", "IL12RB2.rna")
p1 <- ggplot(avg.Th0.cells, aes(BIFIDOpos, BIFIDOneg)) + geom_point() + ggtitle("Th0") #+ geom_text(aes(label=gene))
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.Th2.cells, aes(BIFIDOpos, BIFIDOneg)) + geom_point() + ggtitle("Th2")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
p3 <- ggplot(avg.Th17.cells, aes(BIFIDOpos, BIFIDOneg)) + geom_point() + ggtitle("Th17")
p3 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE)
p4 <- ggplot(avg.Th1.cells, aes(BIFIDOpos, BIFIDOneg)) + geom_point() + ggtitle("Th1")
p4 <- LabelPoints(plot = p4, points = genes.to.label, repel = TRUE)
p5 <- ggplot(avg.IFNb.cells, aes(BIFIDOpos, BIFIDOneg)) + geom_point() + ggtitle("IFNb")
p5 <- LabelPoints(plot = p5, points = genes.to.label, repel = TRUE)
p6 <- ggplot(avg.iTreg.cells, aes(BIFIDOpos, BIFIDOneg)) + geom_point() + ggtitle("iTreg")
p6 <- LabelPoints(plot = p6, points = genes.to.label, repel = TRUE)
plot_grid(p1, p2, p3, p4, p5, p6)

DE_iTreg <- FindMarkers(cyto_cult_iTreg, ident.1 = "BIFIDOpos", ident.2 = "BIFIDOneg", verbose = FALSE, min.pct = 0, logfc.threshold = 0) #min.pct = 0, logfc.threshold = 0
head(DE_iTreg, n = 15)
FeaturePlot(cyto_cult_iTreg, features = c("GAPDH.rna", "GIMAP7.rna", "CD3G.rna"), split.by = "bifido", max.cutoff = 3, 
            cols = c("grey", "red"))
DE_IFNb <- FindMarkers(cyto_cult_IFNb, ident.1 = "BIFIDOpos", ident.2 = "BIFIDOneg", verbose = FALSE, min.pct = 0, logfc.threshold = 0)
head(DE_IFNb, n = 15)
FeaturePlot(cyto_cult_iTreg, features = c("GAPDH.rna", "IL23R.rna", "HLA-C.rna"), split.by = "concentration", max.cutoff = 3, 
            cols = c("grey", "red"))
DE_Th1 <- FindMarkers(cyto_cult_Th1, ident.1 = "BIFIDOpos", ident.2 = "BIFIDOneg", verbose = FALSE, min.pct = 0, logfc.threshold = 0)
head(DE_Th1, n = 15)
FeaturePlot(cyto_cult_Th1, features = c("GAPDH.rna", "IL12RB2.rna", "HLA-C.rna"), split.by = "bifido", max.cutoff = 3, 
            cols = c("grey", "red"))
DE_Th17 <- FindMarkers(cyto_cult_Th17, ident.1 = "BIFIDOpos", ident.2 = "BIFIDOneg", verbose = FALSE, min.pct = 0, logfc.threshold = 0) #min.diff.pct = 0.1
head(DE_Th17, n = 15)
FeaturePlot(cyto_cult_Th17, features = c("GAPDH.rna", "IL23R.rna", "GIMAP7.rna"), split.by = "bifido", max.cutoff = 3, 
            cols = c("grey", "red"))
DE_Th2 <- FindMarkers(cyto_cult_Th2, ident.1 = "BIFIDOpos", ident.2 = "BIFIDOneg", verbose = FALSE, min.pct = 0, logfc.threshold = 0)
head(DE_Th2, n = 15)
FeaturePlot(cyto_cult_Th2, features = c("GAPDH.rna", "GIMAP7.rna", "HLA-C.rna"), split.by = "bifido", max.cutoff = 3, 
            cols = c("grey", "red"))
DE_Th0 <- FindMarkers(cyto_cult_Th0, ident.1 = "BIFIDOpos", ident.2 = "BIFIDOneg", verbose = FALSE, min.pct = 0, logfc.threshold = 0)
head(DE_Th0, n = 15)
head(FindMarkers(cyto_cult_Th0, ident.1 = "BIFIDOpos", ident.2 = "BIFIDOneg", test.use = "DESeq2", max.cells.per.ident = 50))
FeaturePlot(cyto_cult_Th0, features = c("GAPDH.rna", "GIMAP7.rna", "HLA-C.rna"), split.by = "bifido", max.cutoff = 3, 
            cols = c("grey", "red"))

## Volcano pots DE
devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)

p1 <- EnhancedVolcano(DE_iTreg,
                lab = rownames(DE_iTreg),
                x = 'avg_logFC',
                y = 'p_val',
                title = 'iTreg 1:100: bifido+/-',
                pCutoff = 10e-16,
                FCcutoff = 0,
                drawConnectors = TRUE
                #selectLab = c('IL17A.rna','IL23R.rna', 'IL4.rna', 'IL13.rna', 'CCL20.rna', 'IL12A.rna', 'IL21.rna', 'IL31.rna', 'LEF1.rna')
)
p1 +
  ggplot2::coord_cartesian(xlim=c(-2.5, 2.5)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-2.5,2.5, 1))

library(ggplot2)
library(patchwork)
library(cowplot)
immune.combined <- PBMC.RNA2

Idents(object = immune.combined)
colnames(x = immune.combined[[]])
Idents(object = immune.combined) <- 'cytokine_culture'
levels(x = immune.combined)

cyto_cult_Th0 <- subset(immune.combined, idents = 'Th0') 
Idents(cyto_cult_Th0) <- "bifido"

markers.RNA <- FindMarkers(immune.combined, ident.1='Th0', ident.2 = NULL,
                           only.pos = TRUE,
                           min.pct = 0.25,
                           logfc.threshold = 0.25, # log-FC
                           test.use = 'MAST',
                           verbose = FALSE,
                           assay = 'RNA',
                           latent.vars = 'nCount_RNA') # nUMI as proxy for CDR
top10RNA <- markers.RNA %>% 
  group_by('cytokine_culture') %>%
  top_n(n = 10)

features <- c('STAT1.rna', 'CD48.rna', 'ITGB2.rna', 'NCR3.rna', 'STAT4.rna', 'GAPDH.rna', 'IRF8.rna', 'TNFSF10.rna', 'BCL2.rna', 'CBLB.rna')
features <- c('IL23R.rna', 'IL12RB2.rna')
plots <- VlnPlot(immune.combinedAb, features = 'GAPDH.rna' , split.by = "bifido_conc", group.by = "cytokine_culture", 
                 pt.size = 1, combine = FALSE) 
wrap_plots(plots = plots, ncol = 1)
RidgePlot(immune.combined, features = features ,  group.by = "cytokine_culture") 
