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
library(cowplot)


source("./util_funcs.R")

## Reading the data
## These files were prepared using the python script `get_h5ad_data_kz.ipynb`. 


### b. divergense
input.dir.bdiv <- "../Input/scRNAseqBdivCS/"
bdiv0.count.file <- "bd0hr_count.csv"
bdiv12.count.file <- 'bd12hr_count.csv'


## Reading b. divergense data
bdiv0.count <- read.csv(paste(input.dir.bdiv, bdiv0.count.file, sep = ''))
bdiv12.count <- read.csv(paste(input.dir.bdiv, bdiv12.count.file, sep = ''))

genes0 <- bdiv0.count$X
bd0.expr <- bdiv0.count[,-1]
rownames(bd0.expr) <- genes0


bdiv0.pheno <- data.frame(X = colnames(bd0.expr))
bdiv0.pheno$spp <- 'BDiv0' 

genes12 <- bdiv12.count$X
bd12.expr <- bdiv12.count[,-1]
rownames(bd12.expr) <- genes12


bdiv12.pheno <- data.frame(X = colnames(bd12.expr))
bdiv12.pheno$spp <- 'BDiv12' 



# Set initial Seurat clusters
S.O.bd0 <- CreateSeuratObject(counts = bd0.expr, min.cells = 10, min.features = 100)
#S.O.bd0 <- CreateSeuratObject(counts = bd0.expr)

#VlnPlot(S.O.bd0, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#FeatureScatter(S.O.bd0, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

S.O.bd0 <- subset(S.O.bd0, subset = nFeature_RNA > 200 & nFeature_RNA < 1200 )
S.O.bd0 <- prep_S.O(S.O.bd0)


bd0.clust.info <- data.frame(X=as.character(names(S.O.bd0$orig.ident)),
                            cluster=as.character(S.O.bd0$seurat_clusters))
bdiv0.pheno <- inner_join(bdiv0.pheno, bd0.clust.info, by = 'X')
bdiv0.pheno$cells <- paste('BDiv0', bdiv0.pheno$cluster, sep = '')

bdiv0.pheno$NAME <- paste(bdiv0.pheno$cells, 
                         bdiv0.pheno$X, sep = '_')

colnames(bdiv0.pheno) <- c('Sample', 'spp', 'cluster', 'cells', 'NAME')

bd0.expr.filt <- as.matrix(S.O.bd0[["RNA"]]@data)
bd0.ind <- match(colnames(bd0.expr.filt), bdiv0.pheno$Sample)
colnames(bd0.expr.filt) <- bdiv0.pheno$NAME[bd0.ind]


# Set initial Seurat clusters
S.O.bd12 <- CreateSeuratObject(counts = bd12.expr, min.cells = 10, min.features = 100)
#S.O.bd0 <- CreateSeuratObject(counts = bd0.expr)

#VlnPlot(S.O.bd12, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#FeatureScatter(S.O.bd12, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

S.O.bd12 <- subset(S.O.bd12, subset = nFeature_RNA > 200 & nFeature_RNA < 1200 )
S.O.bd12 <- prep_S.O(S.O.bd12)

bd12.clust.info <- data.frame(X=as.character(names(S.O.bd12$orig.ident)),
                             cluster=as.character(S.O.bd12$seurat_clusters))
bdiv12.pheno <- inner_join(bdiv12.pheno, bd12.clust.info, by = 'X')
bdiv12.pheno$cells <- paste('BDiv12', bdiv12.pheno$cluster, sep = '')

bdiv12.pheno$NAME <- paste(bdiv12.pheno$cells, 
                          bdiv12.pheno$X, sep = '_')

colnames(bdiv12.pheno) <- c('Sample', 'spp', 'cluster', 'cells', 'NAME')

bd12.expr.filt <- as.matrix(S.O.bd12[["RNA"]]@data)
bd12.ind <- match(colnames(bd12.expr.filt), bdiv12.pheno$Sample)
colnames(bd12.expr.filt) <- bdiv12.pheno$NAME[bd12.ind]



all.inds <- c(bd0.ind, bd12.ind)

all.samples <- data.frame(Sample = rep(NA, length(all.inds)), 
                          spp = rep(NA, length(all.inds)), stringsAsFactors = F)

all.samples$Sample <- c(as.character(bdiv0.pheno$Sample[bd0.ind]),
                        as.character(bdiv12.pheno$Sample[bd12.ind]))

all.samples$spp <- c(as.character(bdiv0.pheno$spp[bd0.ind]),
                     as.character(bdiv12.pheno$spp[bd12.ind]))

all.samples$cells <- c(as.character(bdiv0.pheno$cells[bd0.ind]),
                       as.character(bdiv12.pheno$cells[bd12.ind]))



#### Use common genes
comm.genes <- rownames(bd12.expr.filt)[rownames(bd12.expr.filt) %in% rownames(bd0.expr.filt)]
bd12.expr.filt.sub <- bd12.expr.filt[(rownames(bd12.expr.filt) %in% comm.genes), ]
bd0.expr.filt.sub <- bd0.expr.filt[(rownames(bd0.expr.filt) %in% comm.genes), ]



## Combine the datasets
all.expr <- inner_join(data.frame(GeneID=rownames(bd0.expr.filt.sub), bd0.expr.filt.sub), 
                                  data.frame(GeneID = rownames(bd12.expr.filt.sub), bd12.expr.filt.sub), 
                                  by = 'GeneID')


rownames(all.expr) <- all.expr$GeneID
all.expr <- all.expr[,-1]

all.samples$NAME <- paste(all.samples$cells, all.samples$Sample, sep = '_')
inds <- all.samples$NAME %in% colnames(all.expr)
meta.data <- data.frame(spp = all.samples$spp[inds], cells = all.samples$cells[inds])
rownames(meta.data) <- all.samples$NAME[inds]

#Create Seurat Object
all.spp <- CreateSeuratObject(all.expr, meta.data = meta.data, min.cells = 3, min.features = 200)
#VlnPlot(all.spp, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#FeatureScatter(all.spp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


## Split by spp
all.spp.list <- SplitObject(all.spp, split.by = "spp")

for (i in 1:length(all.spp.list)) {
  all.spp.list[[i]] <- prep_S.O(all.spp.list[[i]])
  all.spp.list[[i]] <- subset(all.spp.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 1200 )
}


reference.list <- all.spp.list[c("BDiv0", "BDiv12")]
all.samples.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

all.samples.integrated <- IntegrateData(anchorset = all.samples.anchors, dims = 1:30)


# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(all.samples.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
all.samples.integrated <- ScaleData(all.samples.integrated, verbose = FALSE)
all.samples.integrated <- RunPCA(all.samples.integrated, npcs = 30, verbose = FALSE)
all.samples.integrated <- RunUMAP(all.samples.integrated, reduction = "pca", dims = 1:30)

p1 <- DimPlot(all.samples.integrated, reduction = "pca", group.by = "cells", 
              #split.by = 'spp',
              pt.size = 1,
              shape.by='spp',
              label = TRUE, label.size = 6) + NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p1)

ggsave(filename="../Output/scClockFigs/Integrated_cs_scRNAseq_split.pdf", 
       plot=p1,
       width =12, height = 8, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


p2 <- DimPlot(all.samples.integrated, reduction = "umap", 
              pt.size = 1,
              shape.by='spp',
              label = TRUE, label.size = 6) + NoLegend() + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p2)





#### Individual data plots
p1 <- DimPlot(S.O.bd0.filt, reduction = "pca", group.by = "seurat_clusters", 
              pt.size = 1,
              label = TRUE, label.size = 10) + NoLegend() + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p1)

ggsave(filename="../Output/scClockFigs/Individual_scRNAseq_BDiv0.pdf", 
       plot=p1,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


p2 <- DimPlot(S.O.bd12.filt, reduction = "pca", group.by = "seurat_clusters", 
              pt.size = 1,
              label = TRUE, label.size = 10) + NoLegend() + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p2)

ggsave(filename="../Output/scClockFigs/Individual_scRNAseq_BDiv12.pdf", 
       plot=p2,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


####

library(plotly)


bd_0_12.pc <- getPCA(all.samples.integrated)
bd_0_12.pc$spp <- all.samples.integrated$spp
bd_0_12.pc$orig.ident <- all.samples.integrated$orig.ident
bd_0_12.pc$Sample <- gsub('.*_', '', bd_0_12.pc$Sample)

fig <- plot_ly(bd_0_12.pc, x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~spp, colors = colorRamp(c("red", "blue")))
#colors = c('#BF382A', '#0C4B8E'))
fig <- fig %>% add_markers(size=2)
fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC_1'),
                                   yaxis = list(title = 'PC_2'),
                                   zaxis = list(title = 'PC_3')))

fig


### Cross Comparisons of markers

# ## Differential gene expression
## increase the reolution
S.O.bd0.filt <- prep_S.O(S.O.bd0.filt,res = 0.2)
S.O.bd12.filt <- prep_S.O(S.O.bd12.filt,res = 0.2)


p1 <- DimPlot(S.O.bd0.filt, reduction = "pca", group.by = "seurat_clusters", 
              pt.size = 1,
              label = TRUE, label.size = 10) + NoLegend() + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p1)

ggsave(filename="../Output/scClockFigs/Individual_scRNAseq_BDiv0_high_res.pdf", 
       plot=p1,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

p2 <- DimPlot(S.O.bd12.filt, reduction = "pca", group.by = "seurat_clusters", 
              pt.size = 1,
              label = TRUE, label.size = 10) + NoLegend() + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p2)

ggsave(filename="../Output/scClockFigs/Individual_scRNAseq_BDiv12_high_res.pdf", 
       plot=p2,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

BD0.markers <- FindAllMarkers(object = S.O.bd0.filt, only.pos = TRUE, min.pct = 0)
BD12.markers <- FindAllMarkers(object = S.O.bd12.filt, only.pos = TRUE, min.pct = 0)

BD0.markers$GeneID <- gsub('-', '_', BD0.markers$gene)
BD0.markers.top <- BD0.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
FeaturePlot(object = S.O.bd0.filt, 
            features = BD0.markers.top$gene, 
            cols = c("grey", "blue"), reduction = "pca")


BD12.markers$GeneID <- gsub('-', '_', BD12.markers$gene)
BD12.markers.top <- BD12.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
FeaturePlot(object = S.O.bd12.filt, 
            features = BD12.markers.top$gene, 
            cols = c("grey", "blue"), reduction = "pca")

BD0.markers.sig <- BD0.markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)
BD12.markers.sig <- BD12.markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)

ss <- BD12.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())
print(ss)




## Calculate cross-cluster similarity based on presence of up-regulated genes
BD0.markers.up.genes <- BD0.markers.sig %>% 
  select(GeneID, cluster)

BD12.markers.up.genes <- BD12.markers.sig %>% 
  select(GeneID, cluster)

BD0.markers.up.genes.lists <- BD0.markers.up.genes %>% 
  group_by(cluster) %>% summarise(genes = list(GeneID))

BD12.markers.up.genes.lists <- BD12.markers.up.genes %>% 
  group_by(cluster) %>% summarise(genes = list(GeneID))

BD0.markers.up.genes.lists$dummy <- 1
BD12.markers.up.genes.lists$dummy <- 1



# RH vs pb 10x
BD0.BD12.markers.up.genes.lists <- full_join(BD0.markers.up.genes.lists, 
                                             BD12.markers.up.genes.lists, by = 'dummy')

BD0.BD12.markers.up.genes.lists <- BD0.BD12.markers.up.genes.lists %>% rowwise() %>%
  mutate(common.genes = list(intersect(unlist(genes.x), unlist(genes.y))), 
         num.common = length(intersect(unlist(genes.x), unlist(genes.y))))

cluster.sim.BD0.BD12 <- BD0.BD12.markers.up.genes.lists  %>% 
  select(cluster.x, cluster.y, common.genes, num.common) %>% unnest(common.genes)

colnames(cluster.sim.BD0.BD12) = c('BD0', 'BD12', 'GeneID', 'num.comm.genes')
prod.desc <- read.csv('../Input/BdivCellCycle/ProductDescription/BDvi_Prod_desc.csv')
prod.desc <- prod.desc %>% transmute(GeneID = Gene.ID, Product.Description = Product.Description)

cluster.sim.BD0.BD12 <- left_join(cluster.sim.BD0.BD12, prod.desc, by = 'GeneID')
write.xlsx(cluster.sim.BD0.BD12, '../Output/scClockOut/cluster_sim_2way_PB0_PB12.xlsx')


cluster.sim.BD0.BD12 <- cluster.sim.BD0.BD12 %>% 
  select('BD0', 'BD12','num.comm.genes') %>% distinct()

cluster.sim.BD0.BD12 <- cluster.sim.BD0.BD12 %>% 
  pivot_wider(names_from = 'BD12', values_from = 'num.comm.genes')

cluster.sim.BD0.BD12[is.na(cluster.sim.BD0.BD12)] <- 0

cluster.sim.BD0.BD12 <- cluster.sim.BD0.BD12 %>% 
  pivot_longer(-c('BD0'), names_to = 'BD12', values_to = 'num.comm.genes')

hm.palette <- colorRampPalette(brewer.pal(9, 'Blues'), space='Lab')
p <- ggplot(cluster.sim.BD0.BD12 , aes(x = BD0, y = BD12, fill = num.comm.genes)) + 
  geom_tile() + theme_bw() + 
  geom_text(aes(label = num.comm.genes, size=6)) +
  scale_y_discrete(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  ylab("BD12") + xlab("BD0") + 
  scale_fill_gradientn(colours = hm.palette(5)) +
  theme(
    #axis.text.x = element_blank(),
    #axis.ticks = element_blank(),
    legend.position = "none") +
  theme(strip.text.x = element_text(size=14, face="bold"))+
  # strip.background = element_rect(colour="black", fill=NA))+
  #theme_bw()+
  #facet_grid(.~RH, scales = "free_x", space='free_x', 
  #           labeller=label_wrap_gen(multi_line = TRUE)) + 
  #theme(panel.spacing = unit(0.3, "lines")) + 
  #theme(strip.text.x=element_text(angle=0, hjust=0.5, vjust=0.5))+
  
  theme(axis.text.x = element_text(face="bold", size=10, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=10, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )+
  theme(legend.position = "None")
plot(p)  


ggsave(filename="../Output/scClockFigs/cluster_sim_2way_BD0_BD12.pdf", 
       plot=p,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


### Check out some markers

BD0.top2 <- BD0.markers.sig %>% group_by(cluster)  %>% top_n(n=1, avg_log2FC)

#### Expression Heatmap Plots:
my.genes <- BD0.top2$gene


p1 <- FeaturePlot(object = S.O.bd0.filt, 
            label = T, pt.size = 0.6, label.size = 5, 
            features = my.genes,
            cols = c("grey", "blue"), reduction = "pca")

plot(p1)
ggsave(filename="../Output/scClockFigs/top_markers_BD0.pdf", 
       plot=p1,
       width = 8, height = 12, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

p2<- FeaturePlot(object = S.O.bd12.filt, 
                  label = T, pt.size = 0.6, label.size = 5, 
                  features = my.genes,
                  cols = c("grey", "blue"), reduction = "pca")

plot(p2)
ggsave(filename="../Output/scClockFigs/top_markers_BD12.pdf", 
       plot=p2,
       width = 8, height = 12, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


## Check out markers of cluster 3 & 4

expr.bd0.df <- getNormExpr(S.O.bd0.filt)

gene.id <- 'Bdiv_023510'
gene.heat <- expr.bd0.df  %>% dplyr::filter(GeneID == gene.id) 
bd_0_12.pc.heat <- left_join(bd_0_12.pc, gene.heat, by = 'Sample')
fig <- plot_ly(bd_0_12.pc.heat, x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~expr, colors = colorRamp(c("gray", "blue")))
#colors = c('#BF382A', '#0C4B8E'))
fig <- fig %>% add_markers(size=2)
fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC_1'),
                                   yaxis = list(title = 'PC_2'),
                                   zaxis = list(title = 'PC_3')))

fig
