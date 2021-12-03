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


source('util_funcs.R')

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
S.O.sub <- prep_S.O(S.O, res = 0.1, var.features = T, down.sample = T, smooth.data = T, network = T)



pheno <- data.frame(Sample = names(S.O.sub$orig.ident))
spp <- 'BDiv0hr'
pheno$spp <- spp
pheno$NAME <- paste(pheno$spp, pheno$Sample, sep = '_')
rownames(pheno) <- names(S.O.sub$orig.ident)
S.O.sub <- AddMetaData(S.O.sub, metadata = pheno) 

Ft <- S.O.sub@assays$smooth@data
genes <- rownames(Ft)
cells <- colnames(Ft)

Ft.scale <- scale(t(Ft))
S <- cov(Ft.scale)
glassoFast()
saveRDS(S, '../Input/scClock/glasso/S_smooth_3200.RData')
saveRDS(S.O.sub, '../Input/scClock/glasso/S.O.sub_smooth_3200.RData')


## Find variable Features and fit a KNN graph on PCA coordinates
S.O.sub <- subset(x = S.O, downsample = 3200)
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



####

### b. divergens
input.dir.bdiv <- "../Input/scRNAseqBdiv/"
bdiv.count.file <- "bdiv.expr.csv"
## Reading b. divergense data
bdiv.count <- read.csv(paste(input.dir.bdiv, bdiv.count.file, sep = ''))
genes <- bdiv.count$X
bd.expr <- bdiv.count[,-1]
rownames(bd.expr) <- genes
# Set initial Seurat clusters
S.O.bd <- CreateSeuratObject(counts = bd.expr, min.cells = 10, min.features = 100)
S.O.bd <- subset(S.O.bd, subset = nFeature_RNA > 200 & nFeature_RNA < 1200 )
S.O.bd <- prep_S.O(S.O.bd, res = 0.1, var.features = T, down.sample = T, smooth.data = T, network = T)
Ft.bd <- S.O.bd@assays$smooth@data
Ft.bd.scale <- scale(t(Ft.bd))
S.bd <- cov(Ft.bd.scale)

saveRDS(S.bd, '../Input/scClock/glasso/S_bd_smooth_2000.RData')

### T. gondii
input.dir.tg <- '../Input/scRNAseqToxoAtlas/kz/'
rh384.expr.file <- 'rh384_expression_filtered.csv'
## Reading T. gondii data
rh384.expr <- read.csv(paste(input.dir.tg, rh384.expr.file, sep = ''))

processCounts <- function(expr){
  cols <- expr[,1]
  expr <- t(expr[,-1])
  colnames(expr) <- cols
  return(expr)
}

rh384.expr <- processCounts(rh384.expr)
# Set initial Seurat clusters
S.O.tg <- CreateSeuratObject(counts = rh384.expr, min.cells = 10, min.features = 100)
#VlnPlot(S.O.tg, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#FeatureScatter(S.O.tg, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
S.O.tg <- subset(S.O.tg, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 )
S.O.tg <- prep_S.O(S.O.tg, res = 0.1, var.features = T, down.sample = T, smooth.data = T, network = T)
Ft.tg <- S.O.tg@assays$smooth@data
Ft.tg.scale <- scale(t(Ft.tg))
S.tg <- cov(Ft.tg.scale)
saveRDS(S.tg, '../Input/scClock/glasso/S_tg_smooth_2000.RData')


