library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(tidytext)
library(RColorBrewer)
library(cowplot)
library(gam)
library(princurve)
library(parallel)
library(sme)
library(plotly)

source('./loadRDS.R')
source('./util_funcs.R')
source('./DEAsyncTimeCourse.R')


egress_markers <- read.xlsx('../Input/BdivCellCycle/gene_function/sync_egress_kz.xlsx')


prod.desc[grep('calcium',prod.desc$Product.Description),]
CDPK4.id <- "Bdiv_024410"
prod.desc[grep('GMP',prod.desc$Product.Description),]
PKG.id <- "Bdiv_020500"

sc.tc.mus <- smoothSplineSmeFits(sc.tc.fits, unique(sc.tc.df.adj$variable), extend = F)
colnames(sc.tc.mus) <- c('GeneID', 't', 'y')
sc.tc.mus <- sc.tc.mus %>% dplyr::filter(GeneID %in% c(BD.markers.sig$GeneID, PKG.id) )%>%
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  as.data.frame()

## Generate the clusters
num.clust <- 4L


sc.hc_dtw <- dtwClustCurves(sc.tc.mus, nclust = num.clust)


sc.tc.mus.scale <- sc.tc.mus
sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)] <- scale(sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)],
                                                   center = T,scale = T)
sc.tc.mus.scale <- sc.tc.mus.scale %>%  as.data.frame() %>% 
  pivot_longer(!starts_with('t'), names_to = 'GeneID', values_to = 'y')



sc.hc_dtw.df <- data.frame(GeneID = unique(sc.tc.mus.scale$GeneID), 
                           order = as.numeric(sc.hc_dtw$order),
                           cluster = cutree(sc.hc_dtw,k = num.clust))

sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.hc_dtw.df, by = 'GeneID')

sc.tc.mus.scale$GeneID <- factor(as.character(sc.tc.mus.scale$GeneID),
                                 levels = unique(sc.tc.mus.scale$GeneID[sc.tc.mus.scale$order]))



hc_eucledian.df <- withinCalssReOrder(sc.tc.mus.scale) 
sc.tc.mus.scale <- left_join(sc.tc.mus.scale, hc_eucledian.df, by = 'GeneID')
sc.tc.mus.scale <- sc.tc.mus.scale %>% mutate(cluster = as.factor(cluster),
                                              GeneID.reord = reorder_within(GeneID, hc_eucledian.order, cluster)) 

## map the clusteres
sc.overlap <- matchClustersToPhase(sc.hc_dtw.df, BD.markers.sig)

sc.phase.match <- sc.overlap %>% group_by(cluster) %>% summarize(phase = markers[which.max(percent)])
sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.phase.match, by = 'cluster')



p4 <- ggplot(sc.tc.mus.scale, aes(x = t, y = GeneID, fill = y)) + 
  geom_tile() + 
  #scale_x_discrete(expand=c(0,0)) +
  ylab("Genes") + xlab("time/cells") + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(
    axis.text.x = element_blank(),
    #axis.ticks = element_blank(),
    #axis.text.y  = element_blank(),
    legend.position = "none") +
  facet_grid(phase~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE)) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=8, face="bold")
  ) + scale_y_discrete(breaks = c(CDPK4.id, PKG.id), labels = c('CDPK4', 'PKG'))

plot(p4)  


ggsave(filename="../Output/scClockFigs/curve_cluster_heatmap_sc_CDPK4_PKG.pdf",
       plot=p4,
       width = 5, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 600
)

## CDPK4 trend

pdf(file = "../Output/scClockFigs/CDPK4_expression_curve.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 5) # The height of the plot in inches

par(mfrow = c(1,2))
v <- CDPK4.id
ind <- which(unique(sc.tc.df.adj$variable) == v)
plot.sme(sc.tc.fits[[ind]], paste('sc', v))
ind <- which(unique(sync.tc.df$variable) == v)
plot.sme(sync.tc.fits[[ind]], paste('sync', v), conf = F)

dev.off()

pdf(file = "../Output/scClockFigs/PKG_expression_curve.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 5) # The height of the plot in inches


par(mfrow = c(1,2))
v <- PKG.id
ind <- which(unique(sc.tc.df.adj$variable) == v)
plot.sme(sc.tc.fits[[ind]], paste('sc', v))
ind <- which(unique(sync.tc.df$variable) == v)
plot.sme(sync.tc.fits[[ind]], paste('sync', v), conf = F)

dev.off()


#### Expression
Idents(S.O.bd.filt) <- 'seurat_clusters'
DefaultAssay(S.O.bd.filt) <- "RNA"


p2 <- FeaturePlot(object = S.O.bd.filt, 
                  label = F, pt.size = 0.6, label.size = 3, 
                  features = gsub('_', '-', CDPK4.id),
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p2)

p3 <- FeaturePlot(object = S.O.bd.filt, 
                  label = F, pt.size = 0.6, label.size = 3, 
                  features = gsub('_', '-', PKG.id),
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p3)



ggsave(filename="../Output/scClockFigs/expression_CDPK4.pdf",
       plot=p2,
       width = 5, height = 5,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

ggsave(filename="../Output/scClockFigs/expression_PKG.pdf",
       plot=p3,
       width = 5, height = 5,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



####

p6 <- ggplot(sc.tc.mus.scale, aes(x = t, y = GeneID, fill = y)) + 
  geom_tile() + 
  #scale_x_discrete(expand=c(0,0)) +
  ylab("Genes") + xlab("time/cells") + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(
    axis.text.x = element_blank(),
    #axis.ticks = element_blank(),
    #axis.text.y  = element_blank(),
    legend.position = "none") +
  facet_grid(phase~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE)) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  scale_y_discrete(breaks = egress_markers$GeneID, labels = egress_markers$Pf_NAMES, guide=guide_axis(n.dodge=3)) + 
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=5, face="bold", color = 'black')
  ) 

plot(p6)  


ggsave(filename="../Output/scClockFigs/curve_cluster_heatmap_sc_egress.pdf",
       plot=p6,
       width = 6, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 600
)


###
S.O.bd.filt@meta.data$Sample <- rownames(S.O.bd.filt@meta.data)
PKG.expressing.ind <- which(S.O.bd.filt[["RNA"]]@data[grep(gsub('_', '-', PKG.id), rownames(S.O.bd.filt[["RNA"]]@data)),] > log2(3))
PKG.expressing <- colnames(S.O.bd.filt[["RNA"]]@data)[PKG.expressing.ind]
S.O.bd.filt@meta.data$PKG <- ifelse(S.O.bd.filt@meta.data$Sample %in% PKG.expressing, 1, 0)
S.O.bd.filt@meta.data$PKG_clust <- paste(S.O.bd.filt@meta.data$PKG, S.O.bd.filt@meta.data$seurat_clusters, sep = '_') 
Idents(S.O.bd.filt) <- 'PKG_clust'
PKG.markers <- FindAllMarkers(object = S.O.bd.filt, only.pos = TRUE, min.pct = 0)

PKG.markers$GeneID <- gsub('-', '_', PKG.markers$gene)
PKG.markers.top <- PKG.markers %>% group_by(cluster) %>% top_n(4, avg_log2FC)
PKG.markers.sig <- PKG.markers %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05 & cluster == 1)
PKG.markers.sig <- left_join(PKG.markers.sig, prod.desc, by = 'GeneID')
write.xlsx(PKG.markers.sig, '../Output/scClockOut/PKG_expressing_cells_markers.xlsx')
FeaturePlot(object = S.O.bd.filt, 
            features = PKG.markers.top$gene, 
            cols = c("grey", "blue"), reduction = "pca")


WhichCells(S.O.bd.filt, slot = 'data', expression = gsub('_', '-', PKG.id) > 0 )

subset(S.O.bd.filt, subset = gsub('_', '-', PKG.id) > 1)

xx<-colnames(S.O.bd.filt)
pc.sds.adj <- readRDS('../Input/scClock/pc.sds.adj.RData')

pc.sds.adj$seurat_clusters <- S.O.bd.filt@meta.data$seurat_clusters[match(pc.sds.adj$Sample, S.O.bd.filt@meta.data$Sample)]
pc.sds.adj <- pc.sds.adj %>% mutate(phase = ifelse(seurat_clusters == 1, 'S/M', 
                                           ifelse(seurat_clusters == 0, 'G1a',
                                                  ifelse(seurat_clusters == 3, 'G1b','G1c'))))
p <- ggplot(pc.sds.adj, aes(x=PC_1,y=PC_2)) +
  geom_point(aes(
    fill = phase
  ), shape=21, size = 1.5)+
  geom_path(aes(x=sc1[cell.ord],y=sc2[cell.ord])) +
  geom_point(aes(x=sc1[order(adj.time)][1],y=sc2[order(adj.time)][1]), col = 'black', shape=21, size = 5, stroke = 1.2)+
  geom_point(aes(x=sc1[order(adj.time)][1],y=sc2[order(adj.time)][1]), col = 'blue', shape=8, size = 4, stroke = 1.1)+
  theme_bw(base_size = 14) +
  # theme(legend.direction="horizontal",
  #       legend.position = c(0.6,0.99)) +
  theme(legend.position = "right") +
  
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

plot(p)


ggsave(filename="../Output/scClockFigs/PCA_time.pdf",
       plot=p,
       width = 6, height = 5,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

##### Matched marker analysis


S.O.bd.filt@meta.data$Sample <- rownames(S.O.bd.filt@meta.data)
PKG.expressing.ind <- which(S.O.bd.filt[["RNA"]]@data[grep(gsub('_', '-', PKG.id), rownames(S.O.bd.filt[["RNA"]]@data)),] > 0)
PKG.expressing <- colnames(S.O.bd.filt[["RNA"]]@data)[PKG.expressing.ind]
S.O.bd.filt@meta.data$PKG <- ifelse(S.O.bd.filt@meta.data$Sample %in% PKG.expressing, 1, 0)
S.O.bd.filt@meta.data$PKG_clust <- paste(S.O.bd.filt@meta.data$PKG, S.O.bd.filt@meta.data$seurat_clusters, sep = '_') 
Idents(S.O.bd.filt) <- 'PKG_clust'


objs <- unique(S.O.bd.filt@meta.data$PKG)
ref.obj <- objs[1]
query.obj <- objs[2]
clusters <- unique(S.O.bd.filt@meta.data$seurat_clusters)
ref.clusts <- data.frame(ref = paste(ref.obj, clusters, sep = '_'))
ref.clusts$cluster <- gsub('.*_', '', ref.clusts$ref)
ref.clusts$dummy <- 1
query.clusts <- data.frame(query = paste(rep(query.obj, each = length(clusters)), clusters, sep = '_'))
query.clusts$dummy <- 1
query.clusts$cluster <- gsub('.*_', '', query.clusts$query)

contrasts <- full_join(ref.clusts, query.clusts, by = 'dummy') %>% 
  dplyr::filter(cluster.x == cluster.y) %>%
  transmute(ref = ref, query = query)

contrasts$cluster <- gsub('.*_', '', contrasts$ref)

#ident.1 case, ident.2 is control
Idents(S.O.bd.filt) <- 'PKG_clust'
matched.DEGs <- mclapply(split(contrasts, seq(nrow(contrasts))), function(x){
  tmp <- FindMarkers(S.O.bd.filt, ident.1 = x$query, ident.2 = x$ref, verbose = T)
  tmp$ref <- x$ref
  tmp$query <- x$query
  tmp$cluster <- x$cluster
  tmp$gene <- rownames(tmp)
  tmp$GeneID <- gsub('-', '_', tmp$gene)
  return(tmp)
})


matched.DEGs <- bind_rows(matched.DEGs)
PKG.markers <- matched.DEGs %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)
PKG.markers$GeneID <- gsub('-', '_', PKG.markers$gene)
PKG.markers <- left_join(PKG.markers, prod.desc, by = 'GeneID')
write.xlsx(PKG.markers, '../Output/scClockOut/PKG_expressing_cells_markers_matched_clusters.xlsx')
FeaturePlot(object = S.O.bd.filt, 
            features = unique(PKG.markers$gene), 
            cols = c("grey", "blue"), reduction = "pca")

###################3
## Do a network smoothing
## Smoothing over KNN
S.O.bd.filt <- FindNeighbors(S.O.bd.filt, reduction = "pca", dims = 1:10)

AdjacencyMat <- as.matrix(S.O.bd.filt@graphs$RNA_nn)

con <- colSums(AdjacencyMat)
con <- unlist(lapply(con, function(x) ifelse(x == 0, 0, 1/x)))

## Scaling adjacancy
for(i in 1:ncol(AdjacencyMat)){
  AdjacencyMat[,i] <- con[i] * AdjacencyMat[,i]
}

## Get the expression data
expr.norm <- as.matrix(S.O.bd.filt[["RNA"]]@data)
expr.norm.filt <- expr.norm[which(rownames(expr.norm) %in% gsub('_', '-', c(CDPK4.id, PKG.id))), ]

## Smoothing for a few iterations
max.smoothing <- 20
alpha <- 0.7
scDat <- expr.norm.filt
Ft <- scDat
for(i in 1:max.smoothing){
  Ft <- alpha * Ft %*% AdjacencyMat + (1 - alpha) * scDat
}

S.O.smooth <- S.O.bd.filt
smooth_assay <- CreateAssayObject(counts = Ft)
S.O.smooth[["smooth"]] <- smooth_assay


Assays(S.O.smooth)
Idents(S.O.smooth) <- 'seurat_clusters'

DefaultAssay(S.O.smooth) <- "RNA"

p3 <- FeaturePlot(object = S.O.smooth, 
                  label = F, pt.size = 0.6, label.size = 3, 
                  features = gsub('_', '-', CDPK4.id),
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p3)



DefaultAssay(S.O.smooth) <- "smooth"

p4 <- FeaturePlot(object = S.O.smooth, 
                  label = F, pt.size = 0.6, label.size = 3, 
                  features = gsub('_', '-', CDPK4.id),
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p4)





#####
S.O.ave <- S.O.bd.filt
ave_assay <- CreateAssayObject(counts = scDat)
S.O.ave[["ave"]] <- ave_assay


Assays(S.O.ave)
Idents(S.O.ave) <- 'seurat_clusters'

DefaultAssay(S.O.ave) <- "RNA"

p3 <- FeaturePlot(object = S.O.ave, 
                  label = F, pt.size = 0.6, label.size = 3, 
                  features = gsub('_', '-', CDPK4.id),
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p3)

DefaultAssay(S.O.ave) <- "ave"

p4 <- FeaturePlot(object = S.O.ave, 
                  label = F, pt.size = 0.6, label.size = 3, 
                  features = gsub('_', '-', CDPK4.id),
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p4)

