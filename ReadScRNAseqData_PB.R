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


### P Berghei
input.dir.10x <- "../Input/MalariaCellAtlas/Expression_Matrices/10X/pb10xIDC/"
plasmodium.10x.count.file <- 'pb10xIDC_counts.csv'
plasmodium.10x.pheno.file <- 'pb10xIDC_pheno.csv'


## Reading Plasmodium data 10x pberghei
plasmodium.10x.count <- read.csv(paste(input.dir.10x, plasmodium.10x.count.file, sep = ''))
plasmodium.10x.pheno <- read.csv(paste(input.dir.10x, plasmodium.10x.pheno.file, sep = ''))

genes.10x <- plasmodium.10x.count$X
pb.10x.expr <- plasmodium.10x.count[,-1]
rownames(pb.10x.expr) <- genes.10x

# Figure 3B
plasmodium.10x.pheno <- plasmodium.10x.pheno %>% 
  mutate(cells = case_when(absclust == 0 ~ "TrpE",
                           absclust == 1 ~ "TrpM",
                           absclust == 2 ~ "RngL",
                           absclust == 3 ~ "TrpL",
                           absclust == 4 ~ "SchE",
                           absclust == 5 ~ "SchL",
                           absclust == 6 ~ "RngE",
                           absclust == 7 ~ "SchM"))



plasmodium.10x.pheno$spp <- 'PBer'
plasmodium.10x.pheno$NAME <- paste(plasmodium.10x.pheno$cells, plasmodium.10x.pheno$X, sep = '_')

plasmodium.pheno <- plasmodium.10x.pheno %>% transmute(Sample = X, spp = spp, cells = cells, Name = NAME)

# Set initial Seurat clusters
S.O.pb <- CreateSeuratObject(counts = pb.10x.expr, min.cells = 10, min.features = 100)


VlnPlot(S.O.pb, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.pb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

S.O.pb <- subset(S.O.pb, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 )
S.O.pb <- prep_S.O(S.O.pb)
plasmodium.pheno <- plasmodium.pheno %>% dplyr::filter(Sample %in% colnames(S.O.pb))
rownames(plasmodium.pheno) <- plasmodium.pheno$Sample
S.O.pb <- AddMetaData(S.O.pb, plasmodium.pheno)

## down-sample the data to make it more manageable
set.seed(100)
S.O.pb.filt <- subset(x = S.O.pb, downsample = 800)

saveRDS(S.O.pb.filt, '../Input/scClock/S.O.pb.filt.RData')

# ## Differential gene expression
Idents(S.O.pb.filt) <- 'cells'
Pb.markers <- FindAllMarkers(object = S.O.pb.filt, only.pos = TRUE, min.pct = 0)

Pb.markers$GeneID <- gsub('-', '_', Pb.markers$gene)
Pb.markers.top <- Pb.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
FeaturePlot(object = S.O.pb.filt, 
            features = Pb.markers.top$gene, 
            cols = c("grey", "blue"), reduction = "pca")

Pb.markers.sig <- Pb.markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)

saveRDS(Pb.markers.sig, '../Input/scClock/PB.markers.sig.RData')
ss <- Pb.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())
print(ss)

