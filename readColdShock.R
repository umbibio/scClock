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

set.seed(100)


processCount <- function(input.dir, filename, tt, rr, down.sample = T){
  file.counts <- read.csv(paste(input.dir, filename, sep = ''))
  genes <- file.counts$X
  expr <- file.counts[,-1]
  rownames(expr) <- genes
  
  pheno <- data.frame(X = colnames(expr))
  spp <- paste('BDiv', tt, rr, sep = '')
  pheno$spp <- spp
  pheno$time <- tt
  pheno$recover <- rr
  
  S.O <- CreateSeuratObject(counts = expr, min.cells = 10, min.features = 100)
  S.O <- subset(S.O, subset = nFeature_RNA > 200 & nFeature_RNA < 1200 )
  S.O <- prep_S.O(S.O, res = 0.2)
  
  if(down.sample){
    S.O <- subset(x = S.O, downsample = 600)
  }
  
  
  clust.info <- data.frame(X=as.character(names(S.O$orig.ident)),
                           cluster=as.character(S.O$seurat_clusters))
  pheno <- inner_join(pheno, clust.info, by = 'X')
  pheno$cells <- paste(pheno$spp, pheno$cluster, sep = '')
  
  pheno$NAME <- paste(pheno$spp, pheno$X, sep = '_')
  
  colnames(pheno) <- c('Sample', 'spp', 'time', 'recovere', 'cluster', 'cells', 'NAME')
  
  expr.filt <- as.matrix(S.O[["RNA"]]@data)
  ind <- match(colnames(expr.filt), pheno$Sample)
  colnames(expr.filt) <- pheno$NAME[ind]
  
  
  L <- list(pheno = pheno, S.O = S.O, ind = ind, expr.filt = expr.filt)
  
  return(L)
}

mergeS.O <- function(L){
  num.objs <- length(L)
  
  phenos <- lapply(L, `[[`, 1)
  inds <- lapply(L, `[[`, 3)
  exprs.filt <- lapply(L, `[[`, 4)
  
  all.inds <- c(unlist(inds))
  
  phenos.filt <- lapply(1:num.objs, function(i) phenos[[i]][inds[[i]],])
  
  all.samples <- bind_rows(phenos.filt)
  
  
  #### Use common genes
  genes <- lapply(exprs.filt, function(x) rownames(x))
  comm.genes <- Reduce(intersect, genes)
  
  exprs.filt.sub <- lapply(exprs.filt, function(x) data.frame(GeneID = rownames(x)[(rownames(x) %in% comm.genes)], 
                                                              x[(rownames(x) %in% comm.genes), ]))
  all.expr <- exprs.filt.sub %>% purrr::reduce(inner_join, by = "GeneID")
  
  
  
  rownames(all.expr) <- all.expr$GeneID
  all.expr <- all.expr[,-1]
  
  inds <- all.samples$NAME %in% colnames(all.expr)
  meta.data <- all.samples[inds, ]
  rownames(meta.data) <- meta.data$NAME
  
  #Create Seurat Object
  all.spp <- CreateSeuratObject(all.expr, meta.data = meta.data)
  return(all.spp)  
}

processeMerged <- function(S.O.list, file.info, ref.ind){
  all.spp <- mergeS.O(S.O.list[ref.ind])
  ## Split by spp
  all.spp.list <- SplitObject(all.spp, split.by = "spp")
  
  for (i in 1:length(all.spp.list)) {
    all.spp.list[[i]] <- prep_S.O(all.spp.list[[i]], res = 0.2)
    all.spp.list[[i]] <- subset(all.spp.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 1200 )
  }
  
  
  reference.list <- all.spp.list[file.info$spp[ref.ind]]
  all.samples.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
  
  all.samples.integrated <- IntegrateData(anchorset = all.samples.anchors, dims = 1:30)
  
  
  # switch to integrated assay. The variable features of this assay are automatically
  # set during IntegrateData
  DefaultAssay(all.samples.integrated) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  all.samples.integrated <- ScaleData(all.samples.integrated, verbose = FALSE)
  all.samples.integrated <- RunPCA(all.samples.integrated, npcs = 30, verbose = FALSE)
  all.samples.integrated <- RunUMAP(all.samples.integrated, reduction = "pca", dims = 1:30)
  
  L <- list(all.samples.integrated = all.samples.integrated, all.spp.list = all.spp.list)
}

getMarkers <- function(all.spp.list){
  all.markers.list <- mclapply(all.spp.list, function(x) FindAllMarkers(object = x, only.pos = TRUE, min.pct = 0))
  all.markers.list <- lapply(all.markers.list, function(x) {
    x$GeneID = gsub('-', '_', x$gene)
    return(x)
  })
  
  all.markers.list.top <- lapply(all.markers.list, function(x) {
    top.marker = x %>% group_by(cluster) %>% top_n(2, avg_log2FC)
    return(top.marker)
  })
  
  all.markers.list.sig <- lapply(all.markers.list, function(x) {
    sig.marker = x %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)
    return(sig.marker)
  })
  
  L <- list(all.markers.list = all.markers.list, 
            all.markers.list.top = all.markers.list.top, all.markers.list.sig = all.markers.list.sig)
  return(L)
}

crossCompareMarkers <- function(marker1, marker2, cond1, cond2){
  marker1 <- marker1 %>% 
    select(GeneID, cluster)
  
  marker2 <- marker2 %>% 
    select(GeneID, cluster)
  
  marker1.lists <- marker1 %>% 
    group_by(cluster) %>% summarise(genes = list(GeneID))
  
  marker2.lists <- marker2 %>% 
    group_by(cluster) %>% summarise(genes = list(GeneID))
  
  
  
  marker1.lists$dummy <- 1
  marker2.lists$dummy <- 1
  
  
  
  # RH vs pb 10x
  marker1.marker2.lists <- full_join(marker1.lists, marker2.lists, by = 'dummy')
  
  marker1.marker2.lists <- marker1.marker2.lists %>% rowwise() %>%
    mutate(common.genes = list(intersect(unlist(genes.x), unlist(genes.y))), 
           num.common = length(intersect(unlist(genes.x), unlist(genes.y))))
  
  cluster.sim.marker1.marker2 <- marker1.marker2.lists  %>% 
    select(cluster.x, cluster.y, common.genes, num.common) %>% unnest(common.genes)
  
  colnames(cluster.sim.marker1.marker2) = c(cond1, cond2, 'GeneID', 'num.comm.genes')
  
  cluster.sim.marker1.marker2 <- cluster.sim.marker1.marker2 %>% 
    select(cond1, cond2,'num.comm.genes') %>% distinct()
  
  cluster.sim.marker1.marker2 <- cluster.sim.marker1.marker2 %>% 
    pivot_wider(names_from = !!cond2, values_from = 'num.comm.genes')
  
  cluster.sim.marker1.marker2[is.na(cluster.sim.marker1.marker2)] <- 0
  
  cluster.sim.marker1.marker2 <- cluster.sim.marker1.marker2 %>% 
    pivot_longer(-cond1, names_to = paste0(cond2), values_to = 'num.comm.genes')
  
  return(cluster.sim.marker1.marker2)
}

prod.desc <- read.csv('../Input/BdivCellCycle/ProductDescription/BDvi_Prod_desc.csv')
prod.desc <- prod.desc %>% transmute(GeneID = Gene.ID, Product.Description = Product.Description)

saveRDS(prod.desc, '../Input/scClock/prod.desc.RData')

### b. divergense
input.dir.bdiv <- "../Input/scRNAseqBdivCS/"
count.files <- list.files(input.dir.bdiv)
num.total.files <- length(count.files)

file.info <- data.frame(time = rep(NA, num.total.files), 
                        recover = rep(NA, num.total.files), 
                        filename = rep(NA, num.total.files))

S.O.list <- list()

for(i in 1:num.total.files){
  file.info$filename[i] <- count.files[i]
  file.info$time[i] <- gsub('OUT', '', gsub('bd', '', strsplit(file.info$filename[i], split = '_')[[1]][1]))
  file.info$recover[i] = ifelse(grepl('OUT', file.info$filename[i]), 'Y', 'N')
  cat(paste('processing file', file.info$filename[i]))
  cat('\n')
  L <- processCount(input.dir.bdiv, file.info$filename[i], file.info$time[i], file.info$recover[i])
  S.O.list <- c(S.O.list, list(L))
}

file.info$spp <- paste('BDiv', file.info$time, file.info$recover, sep = '')

print(file.info)

## Compare 0h vs. 36h and 36h out
L <- processeMerged(S.O.list, file.info, c(1, 3, 4))

all.samples.integrated_0hr_36hrN_36hrY <- L$all.samples.integrated
all.spp.list_0hr_36hrN_36hrY <- L$all.spp.list

saveRDS(all.samples.integrated_0hr_36hrN_36hrY, '../Input/scClock/all.samples.integrated_0hr_36hrN_36hrY.RData')
saveRDS(all.spp.list_0hr_36hrN_36hrY, '../Input/scClock/all.spp.list_0hr_36hrN_36hrY.RData')



p1 <- DimPlot(all.samples.integrated_0hr_36hrN_36hrY, reduction = "pca", group.by = "cells", 
              split.by = 'spp',
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

ggsave(filename="../Output/scClockFigs/Integrated_0hr_36hrN_36hrY_split_pca.pdf", 
       plot=p1,
       width =12, height = 8, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


p2 <- DimPlot(all.samples.integrated_0hr_36hrN_36hrY, reduction = "umap", group.by = "cells", 
              split.by = 'spp',
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

plot(p2)

ggsave(filename="../Output/scClockFigs/Integrated_0hr_36hrN_36hrY_split_umap.pdf", 
       plot=p1,
       width =12, height = 8, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)




print(file.info)

## Compare 0h vs. 36h and 36h out
L <- processeMerged(S.O.list, file.info, c(1, 6, 7))

all.samples.integrated_0hr_7dN_7dY <- L$all.samples.integrated
all.spp.list_0hr_7dN_7dY <- L$all.spp.list

saveRDS(all.samples.integrated_0hr_7dN_7dY, '../Input/scClock/all.samples.integrated_0hr_7dN_7dY.RData')
saveRDS(all.spp.list_0hr_7dN_7dY, '../Input/scClock/all.spp.list_0hr_7dN_7dY.RData')



p1 <- DimPlot(all.samples.integrated_0hr_7dN_7dY, reduction = "pca", group.by = "cells", 
              split.by = 'spp',
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

ggsave(filename="../Output/scClockFigs/Integrated_0hr_7dN_7dY_split_pca.pdf", 
       plot=p1,
       width =12, height = 8, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


p2 <- DimPlot(all.samples.integrated_0hr_7dN_7dY, reduction = "umap", group.by = "cells", 
              split.by = 'spp',
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

plot(p2)

ggsave(filename="../Output/scClockFigs/Integrated_0hr_7dN_7dY_split_umap.pdf", 
       plot=p1,
       width =12, height = 8, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)




####

library(plotly)


pc <- getPCA(all.samples.integrated)
pc$spp <- all.samples.integrated$spp
pc$orig.ident <- all.samples.integrated$orig.ident
pc$Sample <- gsub('.*_', '', pc$Sample)

fig <- plot_ly(pc, x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~spp, colors = colorRamp(c("red", "green", "blue")))
#colors = c('#BF382A', '#0C4B8E'))
fig <- fig %>% add_markers(size=2)
fig <- fig %>% layout(scene = list(xaxis = list(title = 'PC_1'),
                                   yaxis = list(title = 'PC_2'),
                                   zaxis = list(title = 'PC_3')))

fig


### Cross Comparisons of markers

# ## Differential gene expression


p1 <- DimPlot( all.spp.list_0hr_36hrN_36hrY[[1]], reduction = "pca", group.by = "seurat_clusters", 
              pt.size = 1,
              label = TRUE, label.size = 10) + NoLegend() + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p1)

ggsave(filename="../Output/scClockFigs/Individual_scRNAseq_BDiv0_high_res1.pdf", 
       plot=p1,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

p2 <- DimPlot( all.spp.list_0hr_36hrN_36hrY[[2]], reduction = "pca", group.by = "seurat_clusters", 
              pt.size = 1,
              label = TRUE, label.size = 10) + NoLegend() + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p2)

ggsave(filename="../Output/scClockFigs/Individual_scRNAseq_BDiv36h_high_res.pdf", 
       plot=p2,
       width = 6, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)

p3 <- DimPlot( all.spp.list_0hr_36hrN_36hrY[[3]], reduction = "pca", group.by = "seurat_clusters", 
               pt.size = 1,
               label = TRUE, label.size = 10) + NoLegend() + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p3)


### Marker analysis
L <- getMarkers(all.spp.list_0hr_36hrN_36hrY)
markers_0hr_36hrN_36hrY <- L$all.markers.list.sig

saveRDS(markers_0hr_36hrN_36hrY, '../Input/scClock/markers_0hr_36hrN_36hrY.RData')

L <- getMarkers(all.spp.list_0hr_7dN_7dY)
markers_0hr_7dN_7dY <- L$all.markers.list.sig
saveRDS(markers_0hr_7dN_7dY, '../Input/scClock/markers_0hr_7dN_7dY.RData')




# FeaturePlot(object = all.spp.list[[3]], 
#             features = all.markers.list.top[[1]]$gene, 
#             cols = c("grey", "blue"), reduction = "pca")





## Calculate cross-cluster similarity based on presence of up-regulated genes
cluster.sim.BD0.BD36 <- crossCompareMarkers(markers_0hr_36hrN_36hrY[[1]], markers_0hr_36hrN_36hrY[[2]], 'BD0', 'BD36')





#cluster.sim.BD0.BD36 <- left_join(cluster.sim.BD0.BD36, prod.desc, by = 'GeneID')
#write.xlsx(cluster.sim.BD0.BD36, '../Output/scClockOut/cluster_sim_2way_PB0_PB36.xlsx')



# cluster.sim.BD0.BD36 <- cluster.sim.BD0.BD36 %>% 
#   select('BD0', 'BD36','num.comm.genes') %>% distinct()
# 
# cluster.sim.BD0.BD36 <- cluster.sim.BD0.BD36 %>% 
#   pivot_wider(names_from = 'BD36', values_from = 'num.comm.genes')
# 
# cluster.sim.BD0.BD36[is.na(cluster.sim.BD0.BD36)] <- 0
# 
# cluster.sim.BD0.BD36 <- cluster.sim.BD0.BD36 %>% 
#   pivot_longer(-c('BD0'), names_to = 'BD36', values_to = 'num.comm.genes')

hm.palette <- colorRampPalette(brewer.pal(9, 'Blues'), space='Lab')
p <- ggplot(cluster.sim.BD0.BD36 , aes(x = BD0, y = BD36, fill = num.comm.genes)) + 
  geom_tile() + theme_bw() + 
  geom_text(aes(label = num.comm.genes, size=6)) +
  scale_y_discrete(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  ylab("BD36") + xlab("BD0") + 
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
