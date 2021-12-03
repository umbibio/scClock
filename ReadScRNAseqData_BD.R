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
library(parallel)
#library(sctransform)


source('./util_funcs.R')


### b. divergense
input.dir.bdiv <- "../Input/scRNAseqBdiv/"
bdiv.count.file <- "bdiv.expr.csv"


## Reading b. divergense data
bdiv.count <- read.csv(paste(input.dir.bdiv, bdiv.count.file, sep = ''))


genes <- bdiv.count$X
bd.expr <- bdiv.count[,-1]
rownames(bd.expr) <- genes


bdiv.pheno <- data.frame(X = colnames(bd.expr))
bdiv.pheno$spp <- 'BDiv' 


# Set initial Seurat clusters
S.O.bd <- CreateSeuratObject(counts = bd.expr, min.cells = 10, min.features = 100)
#S.O.bd <- CreateSeuratObject(counts = bd.expr)

#VlnPlot(S.O.bd, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#FeatureScatter(S.O.bd, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

S.O.bd <- subset(S.O.bd, subset = nFeature_RNA > 200 & nFeature_RNA < 1200 )
S.O.bd <- prep_S.O(S.O.bd)

bd.clust.info <- data.frame(X=as.character(names(S.O.bd$orig.ident)),
                            cluster=as.character(S.O.bd$seurat_clusters))
bdiv.pheno <- inner_join(bdiv.pheno, bd.clust.info, by = 'X')
bdiv.pheno$cells <- paste('BDiv', bdiv.pheno$cluster, sep = '')

bdiv.pheno$NAME <- paste(bdiv.pheno$cells, 
                         bdiv.pheno$X, sep = '_')

colnames(bdiv.pheno) <- c('Sample', 'spp', 'cluster', 'cells', 'NAME')

bd.ind <- match(colnames(bd.expr), bdiv.pheno$X)
colnames(bd.expr) <- bdiv.pheno$NAME[bd.ind]

## down-sample the data to make it more manageable
set.seed(100)
S.O.bd.filt <- subset(x = S.O.bd, downsample = 800)

#saveRDS(S.O.bd.filt, '../Input/scClock/S.O.bd.filt.smooth.RData')
saveRDS(S.O.bd.filt, '../Input/scClock/S.O.bd.filt.RData')
# ## Differential gene expression
BD.markers <- FindAllMarkers(object = S.O.bd.filt, only.pos = TRUE, min.pct = 0)
 
BD.markers$GeneID <- gsub('-', '_', BD.markers$gene)
BD.markers.top <- BD.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
FeaturePlot(object = S.O.bd.filt, 
            features = BD.markers.top$gene, 
            cols = c("grey", "blue"), reduction = "pca")
 
BD.markers.sig <- BD.markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)

saveRDS(BD.markers.sig, '../Input/scClock/BD.markers.sig.RData')
ss <- BD.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())
print(ss)

