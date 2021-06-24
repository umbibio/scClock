library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(RColorBrewer)
library(sctransform)
library(glassoFast)
library(igraph)
library(ggraph)
library(graphlayouts)



## Read the WT bdiv and create an initial Seurat object
file.counts <- read.csv("../Input/scRNAseqBdivCS/bd0hr_count.csv")
genes <- file.counts$X
expr <- file.counts[,-1]
rownames(expr) <- genes

S.O <- CreateSeuratObject(counts = expr, min.cells = 10, min.features = 100)


#VlnPlot(S.O, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#FeatureScatter(S.O, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

cutoffs <- quantile(S.O$nCount_RNA, probs = c(0.01, 0.9))
print(cutoffs)
S.O <- subset(S.O, subset = nFeature_RNA > cutoffs[1] & nFeature_RNA < cutoffs[2] )

## Downsample the cells
S.O.sub <- subset(x = S.O, downsample = 3000)

pheno <- data.frame(Sample = names(S.O.sub$orig.ident))
spp <- 'BDiv0hr'
pheno$spp <- spp
pheno$NAME <- paste(pheno$spp, pheno$Sample, sep = '_')
rownames(pheno) <- names(S.O.sub$orig.ident)
S.O.sub <- AddMetaData(S.O.sub, metadata = pheno) 

## Find variable Features and fit a KNN graph on PCA coordinates
S.O.sub <- NormalizeData(S.O.sub) 
## Set to 500 genes
S.O.sub <- FindVariableFeatures(S.O.sub, selection.method = "vst", nfeatures = 800)
S.O.sub <- ScaleData(S.O.sub, verbose = FALSE)
S.O.sub <- RunPCA(S.O.sub, npcs = 30, verbose = FALSE)
S.O.sub <- RunUMAP(S.O.sub, reduction = "pca", dims = 1:30)
S.O.sub <- FindNeighbors(S.O.sub, reduction = "pca", dims = 1:10)
S.O.sub <- FindClusters(S.O.sub, resolution = 0.2)



### Do a network smoothing

## Smoothing over KNN
AdjacencyMat <- as.matrix(S.O.sub@graphs$RNA_nn)

con <- colSums(AdjacencyMat)
con <- unlist(lapply(con, function(x) ifelse(x == 0, 0, 1/x)))

## Scaling adjacancy
for(i in 1:ncol(AdjacencyMat)){
  AdjacencyMat[,i] <- con[i] * AdjacencyMat[,i]
}

## Get the expression data
Idents(S.O.sub) <- 'Sample'
exprs <- AverageExpression(S.O.sub, features = VariableFeatures(object = S.O.sub))

## Smoothing for a few iterations
max.smoothing <- 2
alpha <- 0.4
scDat <- exprs$RNA
Ft <- scDat
for(i in 1:max.smoothing){
  Ft <- alpha * Ft %*% AdjacencyMat + (1 - alpha) * scDat
}


sds <- apply(Ft, 1, sd)
rm.ind <- which(sds == 0)

if(length(rm.ind) > 0){
  Ft <- Ft[-rm.ind,]
}

genes <- rownames(Ft)
cells <- colnames(Ft)

Ft.scale <- scale(t(Ft))
S <- cov(Ft.scale)

saveRDS(S, '../Input/scClock/S_800.RData')

