library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(slingshot)
library(gam)
library(princurve)
library(parallel)
library(dtwclust)
library(doParallel)
library(tidyverse)
library(tidytext)
library(fda)
library(sme)
library(MyEllipsefit)


#library(sctransform)

source('./util_funcs.R')

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

## Fit a pseudo-time curve and align using sync data
S.O.integrated <- readRDS('../Input/scClock/S_O_scRNA_scATAC_integrated.RData')
S.O.integrated@meta.data$phase <- S.O.integrated@meta.data$predicted.id
Idents(S.O.integrated) <- 'orig.ident'

atac_sub <- subset(S.O.integrated, ident = 'scATAC')
rna_sub <- subset(S.O.integrated, ident = 'scRNA')


prod.desc <- read.xlsx('../Input/compScBdTgPb/genes/ProductDescription_GT1.xlsx')
prod.desc <- left_join(prod.desc, ME49_TGGT1, by = c('GeneID' = 'TGGT1'))


## Pseudo-time analysis with SLingshot

fitTime <- function(S.O, method = 'rna'){
  pc.tg <- getPCA(S.O)
  sds.data <- getPrinCurve(pc.tg)
  pc.sds.tg <- left_join(pc.tg, sds.data, by = "Sample")
  pc.sds.tg$phase <- S.O@meta.data$phase[match(pc.sds.tg$Sample, rownames(S.O@meta.data))]
  pc.sds.tg$phase <- factor(pc.sds.tg$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
  
  Y <- log2(S.O@assays$RNA@counts + 1)
  var.genes <- names(sort(apply(Y, 1, var),decreasing = TRUE))#[1:1000] 
  Y <- Y[var.genes, ]
  
  pt <- sds.data$pt
  
  ## Map the pseudo-time to 0-12:20 hours 
  #t <- (12 + 1/3) * ((as.numeric(pt) - min(as.numeric(pt)))/(max(as.numeric(pt)) - min(as.numeric(pt))))
  t <- (6 + 1/6) * ((as.numeric(pt) - min(as.numeric(pt)))/(max(as.numeric(pt)) - min(as.numeric(pt))))
  
  sds.data$t <- t
  
  ## time-index cells in 20 min intervals and identify cells in each partition
  ## They will be considered as replicates
  #time.breaks <- seq(1/3, 12 + 1/3, by = 1/3) 
  time.breaks <- seq(1/6, 6 + 1/6, by = 1/6) 
  time.idx <- rep(0, nrow(sds.data))
  
  ind <- which(sds.data$t <= time.breaks[1])
  time.idx[ind] <- 0
  
  for(i in 2:length(time.breaks)){
    ind <- which(sds.data$t > time.breaks[(i-1)] & sds.data$t <= time.breaks[i])
    time.idx[ind] <- i - 1
  }
  
  sds.data$time.idx <- time.idx
  
  ## Update the time to 20 min increments
  #sds.data$t <- (time.idx) * (1/3)
  sds.data$t <- (time.idx) * (1/6)
  
  sds.data <- sds.data %>%  
    group_by(time.idx) %>% mutate(rep = seq(1:n()))
  
  rownames(sds.data) <- sds.data$Sample
  
  
  
  ## Run a GAM regression of expression on the pseudo-time
  ## Use parallel computation to speed things up. 16 cores
  gam.pval <- mclapply(1:nrow(Y), function(z){
    d <- data.frame(z=as.numeric(Y[z,]), t=as.numeric(pt))
    tmp <- gam(z ~ lo(t), data=d)
    #p <- summary(tmp)[4][[1]][1,5]
    p <- summary(tmp)$anova$`Pr(F)`[2]
    p
  }, mc.cores = num.cores)
  
  gam.pval <- unlist(gam.pval)
  names(gam.pval) <- rownames(Y)
  ## Remove the NA's and get the best fits
  if(any(is.na(gam.pval))){
    gam.pval <- gam.pval[-which(is.na(gam.pval))]
  }
  
  gam.pval.adj <- p.adjust(gam.pval, method = 'fdr', n = length(gam.pval))
  gam.pval.sig <- gam.pval[gam.pval.adj < 0.01] 
  print(length(gam.pval.sig)) ## number of correlating genes
  
  ## Sort the cells on the pt
  cell.ord <- sds.data$cell.ord
  
  topgenes <- names(sort(gam.pval.sig, decreasing = FALSE))  
  cell.cycle.genes.expr <- as.matrix(S.O@assays$RNA@data[topgenes, cell.ord])
  #cell.cycle.genes.expr <- as.matrix(S.O.bd.filt@assays$smooth@data[topgenes, cell.ord]) ## smoothed version
  
  
  cell.cycle.genes.df <- data.frame(GeneID = rownames(cell.cycle.genes.expr),
                                    cell.cycle.genes.expr) %>% 
    pivot_longer(-c(GeneID), names_to = 'Sample', values_to = 'log2.expr')
  
  if(method = 'atac'){
    cell.cycle.genes.df$Sample <- gsub('\\.', '-', cell.cycle.genes.df$Sample)
  }
  
  cell.cycle.genes.df$GeneID <- gsub('-', '_', cell.cycle.genes.df$GeneID)
  #cell.cycle.genes.df <- left_join(cell.cycle.genes.df, sds.data, by = 'Sample')
  cell.cycle.genes.df$cluster <- S.O@meta.data$seurat_clusters[match(cell.cycle.genes.df$Sample, 
                                                                        rownames(S.O@meta.data))]
  
  cell.cycle.genes.df$phase <- S.O@meta.data$phase[match(cell.cycle.genes.df$Sample, 
                                                            rownames(S.O@meta.data))]
  cell.cycle.genes.df$phase <- factor(cell.cycle.genes.df$phase, 
                                      levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
  
  
  
  tmp <- left_join(cell.cycle.genes.df, sds.data, by = 'Sample')
  start.times <- tmp %>% group_by(GeneID) %>% arrange(phase, pt) %>% summarise(st = t[1])
  lag.time <- which(sort(unique(tmp$t)) == unique(start.times$st)) + 2 ## all are the same
  
  
  adjusted.time <- (sds.data$time.idx * 1/6) -  sort(unique(sds.data$time.idx) * 1/6)[lag.time]
  neg.ind <- ifelse(adjusted.time < 0, T, F)
  adjusted.time[neg.ind] <- adjusted.time[neg.ind] + (6 + 1/6)
  sds.data$adj.time <-  adjusted.time
  
  
  ## create fine-resolution 20 min cluster of cells
  clusters <- paste('C', 1:length(unique(sds.data$adj.time)), sep = '')
  sds.data$cluster <- clusters[as.integer((sds.data$adj.time) * 6 + 1)]
  
  
  ## Generate shifted curves
  time.breaks <- seq(1/6, 6 + 1/6, by = 1/6) 
  time.idx <- rep(0, nrow(sds.data))
  
  ind <- which(sds.data$adj.time <= time.breaks[1])
  time.idx[ind] <- 0
  
  for(i in 2:length(time.breaks)){
    ind <- which(sds.data$adj.time > time.breaks[(i-1)] & sds.data$adj.time <= time.breaks[i])
    time.idx[ind] <- i - 1
  }
  
  sds.data$adj.time.idx <- time.idx
  
  
  sds.data <- sds.data %>%  ungroup() %>%
    group_by(adj.time.idx) %>% mutate(rep = seq(1:n()))
  
  sds.data <- as.data.frame(sds.data)
  rownames(sds.data) <- sds.data$Sample
  
  ## Add the new clusters as meta-data
  S.O <- AddMetaData(S.O, sds.data)
  
  pc <- S.O[['pca']]@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc))
  pc$cluster <- S.O$cluster
  
  
  
  pc.sds.adj <- left_join(pc, sds.data, by = "Sample")
  
  lvs <- paste('C', unique(sort(as.numeric(gsub('C', '', pc.sds.adj$cluster.y)))), sep = '')
  pc.sds.adj$cluster.y <- factor(pc.sds.adj$cluster.y, levels = lvs)
  
  
  cell.cycle.genes.df.adj <- left_join(cell.cycle.genes.df, sds.data[,c('Sample', 'adj.time', 
                                                                        'adj.time.idx', 'rep', 'cluster')], by = 'Sample')
  sc.tc.df.adj <- cell.cycle.genes.df.adj %>% 
    transmute(y = log2.expr, tme = adj.time, ind = rep, variable = GeneID)
  
  L <- list(cell.cycle.genes.df = cell.cycle.genes.df, 
            cell.cycle.genes.df.adj = cell.cycle.genes.df.adj,
            sc.tc.df.adj = sc.tc.df.adj,
            sds.data = sds.data,
            S.O = S.O,
            gam.genes = names(gam.pval.sig))
  return(L)
}



L.atac <- fitTime(atac_sub, 'atac')
L.rna <- fitTime(rna_sub)





## Individual genes
my.gene <- "TGME49_300100" ## RON2
sme.fit.RON2.rna <- sme(L.rna$sc.tc.df.adj[L.rna$sc.tc.df.adj$variable==my.gene,c("y","tme","ind")],lambda.mu = 6, lambda.v = 6)
sme.fit.RON2.atac <- sme(L.atac$sc.tc.df.adj[L.atac$sc.tc.df.adj$variable==my.gene,c("y","tme","ind")],lambda.mu = 6, lambda.v = 6)

pdf(file = "../Output/scClockFigs/RON2_sme_fits_scATAC_scRNA.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches

par(mfrow = c(2,1))
plot.sme(sme.fit.RON2.rna, "RON2:scRNA")
plot.sme(sme.fit.RON2.atac, "RON2:scATAC")
dev.off()

sme(L.rna$sc.tc.df.adj[L.rna$sc.tc.df.adj$variable==my.gene,c("y","tme","ind")],
    lambda.mu = 6, lambda.v = 6)
Idents(S.O.tg) <- 'phase'
p <- DimPlot(S.O.tg, reduction = "pca", 
             #group.by = "cells", 
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

plot(p)




saveRDS(cell.cycle.genes.df, '../Input/compScBdTgPb/RData/tg_cell_cycle_genes_df_new.RData')
saveRDS(sc.tc.df.adj, '../Input/compScBdTgPb/RData/tg_sc_tc_df_adj_new.RData')


sc.tc.fits <- mclapply(unique(sc.tc.df.adj$variable),
                       function(v)
                         sme(sc.tc.df.adj[sc.tc.df.adj$variable==v,c("y","tme","ind")],
                             lambda.mu = 6, lambda.v = 6), mc.cores = num.cores)

#saveRDS(object = sc.tc.fits, file = "../Input/compScBdTgPb/RData/tg_sme_fits_sc_tc_20min.RData")
saveRDS(object = sc.tc.fits, file = "../Input/compScBdTgPb/RData/tg_sme_fits_sc_tc_10min_new.RData")



## Individual genes
v <- 'TGGT1_250800' ## AP2XII-8
sme.fit.ap2 <- sme(sc.tc.df.adj[sc.tc.df.adj$variable==v,c("y","tme","ind")],lambda.mu = 1, lambda.v = 1)
pdf(file = "../Output/compScBdTgPb/figs/AP2XII-8_sme_fits_sc_new.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
plot.sme(sme.fit.ap2, v)

dev.off()



## Plot expression too
tg_expr <- data.frame(as.matrix(S.O.tg@assays$RNA@data))
tg_expr$GeneID <- gsub('-', '_', rownames(as.matrix(S.O.tg@assays$RNA@data)))
tg_expr <- tg_expr %>% pivot_longer(-GeneID, names_to = 'Sample', values_to = 'expr')
pc <- S.O.tg@reductions$pca@cell.embeddings
pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
  transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)

umap <- S.O.tg@reductions$umap@cell.embeddings
umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
  transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)

meta.data <- data.frame(Sample = rownames(S.O.tg@meta.data), pt = S.O.tg@meta.data$pt, 
                        t = S.O.tg@meta.data$t, Phase = S.O.tg@meta.data$phase, 
                        phase = S.O.tg@meta.data$phase, 
                        time.idx = S.O.tg@meta.data$time.idx, 
                        sc1 = S.O.tg@meta.data$sc1, sc2 = S.O.tg@meta.data$sc2, rep = S.O.tg@meta.data$rep, 
                        cell.ord = S.O.tg@meta.data$cell.ord)

meta.data <- left_join(meta.data,
                       pc, by = 'Sample')
meta.data <- left_join(meta.data, umap, by = 'Sample')

#meta.data$phase <- factor(meta.data$phase, levels = c('G1', 'S', 'M', 'C'))
meta.data$Phase <- factor(meta.data$Phase, levels = c('G1.a','G1.b', 'S', 'M', 'C'))

gg <- 'TGGT1_250800' ## AP2XII-8

tg_expr.filt <- tg_expr %>% dplyr::filter(GeneID == gg)
tg_data <- left_join(meta.data, tg_expr.filt, by = 'Sample')
tg_expr.filt$expr <- exp(tg_expr.filt$expr)
mid_point <- mean(tg_expr.filt$expr)
col_range <- colorRampPalette(c('white', 'blue'))
tg_data$phase <- factor(tg_data$phase, levels = c('G1.a','G1.b', 'S', 'M', 'C'))
p <- ggplot(tg_data, aes(x=PC_1,y=PC_2)) +
  geom_point(aes(
    fill = expr, color = phase,
  ), shape=21, size = 1.5)+ 
  geom_path(aes(x=sc1[cell.ord],y=sc2[cell.ord])) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient2(low = "gray66", high = "blue", midpoint = 0.001) + 
  #scale_fill_brewer(palette = "BuPu") +
  ylab('PC2') + xlab('PC1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/AP2XII-8_single_cel_expr_new.pdf",
       plot=p,
       width = 5, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 600
)

Idents(S.O.tg) <- 'phase'
VlnPlot(S.O.tg, features = 'TGGT1-250800')
## Plot a few curves to check the alignments

vs = unique(sc.tc.df$variable)[1:16]


pdf(file = "../Output/compScBdTgPb/figs/tg_sme_fits_sc.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

par(mfrow = c(4,4))
for(v in vs){
  ind <- which(unique(sc.tc.df$variable) == v)
  plot.sme(sc.tc.fits[[ind]], v)
}

dev.off()




sc.tc.mus <- smoothSplineSmeFits(sc.tc.fits, unique(sc.tc.df.adj$variable), extend = F)
colnames(sc.tc.mus) <- c('GeneID', 't', 'y')
sc.tc.mus <- sc.tc.mus %>% dplyr::filter(GeneID %in% markers.sig$GeneID) %>%
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  as.data.frame()

## Generate the clusters
num.clust <- 8L

sc.hc_dtw <- dtwClustCurves(sc.tc.mus[,2:ncol(sc.tc.mus)], nclust = num.clust)
plot(sc.hc_dtw, type = 'sc')
plot(sc.hc_dtw, type = 'centroids')
plot(sc.hc_dtw, type = "series", clus = 2L)
plot(sc.hc_dtw, type = "centroids", clus = 2L)

saveRDS(sc.hc_dtw, '../Input/compScBdTgPb/RData/tg_sc.hc_dtw_new.RData')

## Scale the data for heatmaps
sc.tc.mus.scale <- sc.tc.mus
sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)] <- scale(sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)],
                                                   center = T,scale = T)
sc.tc.mus.scale <- sc.tc.mus.scale %>%  as.data.frame() %>% 
  pivot_longer(starts_with('TG'), names_to = 'GeneID', values_to = 'y')

ind <- which(sc.tc.mus$variable == gg)
plot(sc.tc.mus$tme[ind], sc.tc.mus$y[ind], type = 'l', col = 'red')
## Add curve cluster info

sc.hc_dtw.df <- data.frame(GeneID = unique(sc.tc.mus.scale$GeneID), 
                           order = as.numeric(sc.hc_dtw$order),
                           cluster = cutree(sc.hc_dtw,k = num.clust))

sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.hc_dtw.df, by = 'GeneID')

#sc.tc.mus.scale$GeneID <- factor(as.character(sc.tc.mus.scale$GeneID),
#                                 levels = unique(sc.tc.mus.scale$GeneID[sc.tc.mus.scale$order]))



## Reorder the genes within each cluster.
#hc_eucledian.df <- withinCalssReOrder(sc.tc.mus.scale) 
#sc.tc.mus.scale <- left_join(sc.tc.mus.scale, hc_eucledian.df, by = 'GeneID')
#sc.tc.mus.scale <- sc.tc.mus.scale %>% mutate(cluster = as.factor(cluster),
#                                              GeneID.reord = reorder_within(GeneID, hc_eucledian.order, cluster)) 
getCurvePeakLoc <- function(t, y){
  
  ## Fitting the estimated kernel with smooth splines
  spline.fit <- smooth.spline(x = t, y = y)
  
  ## Compute the derivatives of the fitted splines
  s.0 <- predict(spline.fit, spline.fit$x, deriv=0)
  s.1 <- predict(spline.fit, spline.fit$x, deriv=1)
  s.derv <- data.frame(s0=s.0$y, s1=s.1$y)
  
  ## Get the location of the extrema
  locs <- rle(den.sign <- sign(s.derv$s1))
  
  ## Maxima
  inc.ind <- which(locs$values == 1)
  if(length(inc.ind) > 1){
    maxima.ind = {}
    for(i in inc.ind){
      maxima.ind = c(maxima.ind, sum(locs$lengths[1:i]))
    }
    ## Interpolate a point between the location where derivative changes sign
    maxima = (spline.fit$x[maxima.ind] + spline.fit$x[(maxima.ind + 1)]) / 2
    maxima = maxima[!is.na(maxima)]
    ## Get the maximum values
    maxval = predict(spline.fit, maxima)
    
    ## Get the outliers
    maxima.outliers = which(maxval$y >= quantile(maxval$y, prob = 0.9))
    
    ## Peaks for entities of interest
    entity.x = maxval$x[maxima.outliers]
    entity.y = maxval$y[maxima.outliers]
  }else{
    entity.x <- spline.fit$x[which.max(spline.fit$y)]
    entity.y <- spline.fit$y[which.max(spline.fit$y)]
  }
  
  return(entity.x)
}

## map the clusters
sc.peak.order <- sc.tc.mus.scale %>% group_by(GeneID) %>% summarise(peak.ord = getCurvePeakLoc(t, y))
sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.peak.order, by = 'GeneID')

#markers.sig$cluster <- gsub('G1.b', 'G1', gsub('G1.a', 'G1', as.character(TG.markers.sig$cluster)))
markers.sig$phase <- factor(markers.sig$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
sc.tc.mus.scale <- inner_join(sc.tc.mus.scale, markers.sig, by = 'GeneID')

p2 <- ggplot(sc.tc.mus.scale, aes(x = t, y = reorder_within(GeneID, -peak.ord, phase), fill = y)) + 
  geom_tile() + 
  #scale_x_discrete(expand=c(0,0)) +
  ylab("Genes") + xlab("time/cells") + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y  = element_blank(),
    legend.position = "none") +
  #facet_grid(phase~., scales = "free",  space='free',
  #          labeller=label_wrap_gen(multi_line = TRUE)) +
  #theme(panel.spacing = unit(0.1, "lines")) + 
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p2)  


ggsave(filename="../Output/compScBdTgPb/figs/tg_curve_cluster_heatmap_scg_g1a_g1b.png",
       plot=p2,
       width = 4, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 600
)

sc.overlap$cluster <- factor(sc.overlap$cluster, levels = c(5, 4, 3, 2, 1))
sc.overlap$markers <- factor(sc.overlap$markers, levels = c('G1', 'S', 'M', 'C'))

p3 <- ggplot(sc.overlap, aes(x = cluster, y = markers, fill = percent)) + 
  geom_tile() + 
  #scale_x_discrete(expand=c(0,0)) +
  ylab("markers") + xlab("clusters") + theme_bw() + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  # theme(
  #   axis.text.x = element_blank(),
  #   axis.ticks = element_blank(),
  #   axis.text.y  = element_blank(),
  #   legend.position = "none") +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p3)  


saveRDS(sc.tc.mus, '../Input/compScBdTgPb/RData/tg_sc.tc.mus.RData')
saveRDS(sc.tc.mus.scale, '../Input/compScBdTgPb/RData/tg_sc.tc.mus.scale.RData')


### Calculate transition times

tmp <- sc.tc.mus.scale %>% ungroup() %>% group_by(phase) %>% summarise(mean.peak = median(peak.ord))
