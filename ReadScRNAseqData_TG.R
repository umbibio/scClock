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
library(Matrix)
#library(sctransform)


source('./util_funcs.R')

### Tg

matrix_dir = "../Input/scRNA_toxo/Tg1_GT1/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)

colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

tg.count <- as.data.frame(as.matrix(mat))



genes <- rownames(tg.count)
tg.expr <- tg.count


tg.pheno <- data.frame(X = colnames(tg.expr))
tg.pheno$spp <- 'Tg' 


# Set initial Seurat clusters
S.O.tg<- CreateSeuratObject(counts = tg.expr, min.cells = 3, min.features = 50)

VlnPlot(S.O.tg, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(S.O.tg, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

S.O.tg <- subset(S.O.tg, subset = nFeature_RNA > 20 & nFeature_RNA < 200 )
S.O.tg <- prep_S.O(S.O.tg)

tg.clust.info <- data.frame(X=as.character(names(S.O.tg$orig.ident)),
                            cluster=as.character(S.O.tg$seurat_clusters))
tg.pheno <- inner_join(tg.pheno, tg.clust.info, by = 'X')
tg.pheno$cells <- paste('Tg', tg.pheno$cluster, sep = '')

tg.pheno$NAME <- paste(tg.pheno$cells, 
                         tg.pheno$X, sep = '_')

colnames(tg.pheno) <- c('Sample', 'spp', 'cluster', 'cells', 'NAME')

tg.ind <- match(colnames(tg.expr), rownames(tg.expr))
colnames(tg.expr) <- tg.pheno$NAME[tg.ind]

## down-sample the data to make it more manageable
set.seed(100)
S.O.tg.filt <- subset(x = S.O.tg, downsample = 1000)
S.O.tg.filt <- S.O.tg

saveRDS(S.O.tg.filt, '../Input/scClock/S.O.tg.filt.RData')
# ## Differential gene expression
tg.markers <- FindAllMarkers(object = S.O.tg.filt, only.pos = TRUE, min.pct = 0)

tg.markers$GeneID <- gsub('-', '_', tg.markers$gene)
tg.markers.top <- tg.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
FeaturePlot(object = S.O.tg.filt, 
            features = tg.markers.top$gene, 
            cols = c("grey", "blue"), reduction = "pca")

tg.markers.sig <- tg.markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)

saveRDS(tg.markers.sig, '../Input/scClock/tg.markers.sig.RData')
ss <- tg.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())
print(ss)

pc.tg <- getPCA(S.O.tg.filt)
#sds.data <- getPrinCurve(pc.tg)
#sds.data <- getSlingShot(S.O.tg.filt, 'pca')
#pc.sds.tg <- left_join(pc.tg, sds.data, by = "Sample")

p1 <- ggplot(pc.tg, aes(x=PC_1,y=PC_2)) + 
  geom_point(aes(
    fill = cluster
  ), shape=21, size = 1.5)+ 
  theme_bw(base_size = 14) + 
  theme(legend.position=c(1,1),legend.justification=c(1,1), 
        legend.title = element_blank(),
        legend.background = element_rect(fill=alpha('white', 0)),
        legend.direction="vertical") + 
  ylab('PC2') + xlab('PC1') + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=0.5, linetype="solid")) + 
  theme(
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  guides(color = FALSE)

ss <- tg.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())

p2 <- ggplot(data=ss, aes(x=cluster, y=num.DEG)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=num.DEG), vjust=1.6, color="white", size=3.5)+
  theme_minimal()

plot_grid(p1, p2, labels=c("", ""), ncol = 2, nrow = 1)
