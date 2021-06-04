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
  
  #pheno <- data.frame(X = colnames(expr))
  pheno <- data.frame(Sample = colnames(expr))
  spp <- paste('BDiv', tt, rr, sep = '')
  pheno$spp <- spp
  pheno$time <- tt
  pheno$reactivate <- rr
  
  S.O <- CreateSeuratObject(counts = expr, min.cells = 10, min.features = 100)
  
  #VlnPlot(S.O, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
  #FeatureScatter(S.O, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  cutoffs <- quantile(S.O$nCount_RNA, probs = c(0.05, 0.92))
  print(cutoffs)
  S.O <- subset(S.O, subset = nFeature_RNA > cutoffs[1] & nFeature_RNA < cutoffs[2] )
  
  if(down.sample){
    S.O <- subset(x = S.O, downsample = 2000)
  }
  
  
  #S.O <- prep_S.O(S.O, res = 0.2)
 
  
  #clust.info <- data.frame(X=as.character(names(S.O$orig.ident)),
  #                         cluster=as.character(S.O$seurat_clusters))
  #pheno <- inner_join(pheno, clust.info, by = 'X')
  #pheno$cells <- paste(pheno$spp, pheno$cluster, sep = '')
  
  #pheno$NAME <- paste(pheno$spp, pheno$X, sep = '_')
  pheno$NAME <- paste(pheno$spp, pheno$Sample, sep = '_')
  
  #colnames(pheno) <- c('Sample', 'spp', 'time', 'reactivate', 'cluster', 'cells', 'NAME')
  colnames(pheno) <- c('Sample', 'spp', 'time', 'reactivate',  'NAME')
  
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
    all.spp.list[[i]] <- NormalizeData(all.spp.list[[i]] , normalization.method = "LogNormalize", scale.factor = 10000)
    all.spp.list[[i]] <- FindVariableFeatures(all.spp.list[[i]] , selection.method = "vst", nfeatures = 2000)
  }
  
  
  reference.list <- all.spp.list[file.info$spp[ref.ind]]
  features <- SelectIntegrationFeatures(object.list = reference.list, )
  
  all.samples.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, anchor.features = features)
  
  all.samples.integrated <- IntegrateData(anchorset = all.samples.anchors, dims = 1:30)
  
  
  # switch to integrated assay. The variable features of this assay are automatically
  # set during IntegrateData
  DefaultAssay(all.samples.integrated) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  all.samples.integrated <- ScaleData(all.samples.integrated, verbose = FALSE)
  all.samples.integrated <- RunPCA(all.samples.integrated, npcs = 30, verbose = FALSE)
  all.samples.integrated <- RunUMAP(all.samples.integrated, reduction = "pca", dims = 1:30)
  all.samples.integrated <- FindNeighbors(all.samples.integrated, reduction = "pca", dims = 1:30)
  all.samples.integrated <- FindClusters(all.samples.integrated, resolution = 0.2)
  
  all.samples.integrated$phase.cond <- paste(all.samples.integrated@meta.data$spp, 
                                             Idents(all.samples.integrated), sep = "_")
  all.samples.integrated$celltype <- Idents(all.samples.integrated)
  Idents(all.samples.integrated) <- "phase.cond"
  
  return(all.samples.integrated)
}

getCellCyclePhaseMarkers <- function(all.spp.list){
  all.markers.list <- mclapply(all.spp.list, function(x) FindAllMarkers(object = x, only.pos = TRUE, min.pct = 0))
  all.markers.list <- lapply(all.markers.list, function(x) {
    x$GeneID = gsub('-', '_', x$gene)
    x$glob.clust <- gsub('.*_', '', x$cluster)
    return(x)
  })
  
  all.markers.list.top <- lapply(all.markers.list, function(x) {
    top.marker = x %>% group_by(cluster) %>% top_n(2, avg_log2FC)
    return(top.marker)
  })
  
  all.markers.list.sig <- lapply(all.markers.list, function(x) {
    sig.marker = x %>% dplyr::filter(avg_log2FC > log2(1) & p_val_adj < 0.01)
    return(sig.marker)
  })
  
  L <- list(all.markers.list = all.markers.list, 
            all.markers.list.top = all.markers.list.top, all.markers.list.sig = all.markers.list.sig)
  return(L)
}

crossCompareMarkers <- function(marker1, marker2, cond1, cond2){
  marker1 <- marker1 %>% 
    select(GeneID, glob.clust)
  
  marker2 <- marker2 %>% 
    select(GeneID, glob.clust)
  
  marker1.lists <- marker1 %>% 
    group_by(glob.clust) %>% summarise(genes = list(GeneID))
  
  marker2.lists <- marker2 %>% 
    group_by(glob.clust) %>% summarise(genes = list(GeneID))
  
  
  
  marker1.lists$dummy <- 1
  marker2.lists$dummy <- 1
  
  
  marker1.marker2.lists <- full_join(marker1.lists, marker2.lists, by = 'dummy')
  
  marker1.marker2.lists <- marker1.marker2.lists %>% rowwise() %>%
    mutate(common.genes = list(intersect(unlist(genes.x), unlist(genes.y))), 
           num.common = length(intersect(unlist(genes.x), unlist(genes.y))))
  
  cluster.sim.marker1.marker2 <- marker1.marker2.lists  %>% 
    select(glob.clust.x, glob.clust.y, common.genes, num.common) %>% unnest(common.genes)
  
  colnames(cluster.sim.marker1.marker2) = c(cond1, cond2, 'GeneID', 'num.comm.genes')
  
  cluster.sim.marker1.marker2 <- cluster.sim.marker1.marker2 %>% 
    select(all_of(cond1), all_of(cond2),'num.comm.genes') %>% distinct()
  
  cluster.sim.marker1.marker2 <- cluster.sim.marker1.marker2 %>% 
    pivot_wider(names_from = !!cond2, values_from = 'num.comm.genes')
  
  cluster.sim.marker1.marker2[is.na(cluster.sim.marker1.marker2)] <- 0
  
  cluster.sim.marker1.marker2 <- cluster.sim.marker1.marker2 %>% 
    pivot_longer(-cond1, names_to = paste0(cond2), values_to = 'num.comm.genes')
  
  return(cluster.sim.marker1.marker2)
}


markerContrasts <- function(marker1, marker2, cond1, cond2){
  marker1 <- marker1 %>% 
    select(GeneID, glob.clust)
  
  marker2 <- marker2 %>% 
    select(GeneID, glob.clust)
  
  marker1.lists <- marker1 %>% 
    group_by(glob.clust) %>% summarise(genes = list(GeneID))
  
  marker2.lists <- marker2 %>% 
    group_by(glob.clust) %>% summarise(genes = list(GeneID))
  
  
  marker1.lists$dummy <- 1
  marker2.lists$dummy <- 1
  
  marker1.marker2.lists <- full_join(marker1.lists, marker2.lists, by = 'dummy')
  marker1.marker2.lists <- marker1.marker2.lists %>% dplyr::filter(glob.clust.x == glob.clust.y)
  
  marker1.marker2.lists <- marker1.marker2.lists %>% rowwise() %>%
    mutate(common.genes = list(intersect(unlist(genes.x), unlist(genes.y))), 
           num.common = length(intersect(unlist(genes.x), unlist(genes.y))),
           cond1.only.genes = list(setdiff(unlist(genes.x), unlist(genes.y))),
           num.cond1.only.genes = length(setdiff(unlist(genes.x), unlist(genes.y))),
           cond2.only.genes = list(setdiff(unlist(genes.y), unlist(genes.x))),
           num.cond2.only.genes = length(setdiff(unlist(genes.y), unlist(genes.x))))
  
  
  marker1.marker2.lists <- marker1.marker2.lists %>% select(-dummy)
  colnames(marker1.marker2.lists) <- c(paste(cond1, 'cluster', sep = '.'), 
                                       paste(cond1, 'genes', sep = '.'),
                                       paste(cond2, 'cluster', sep = '.'), 
                                       paste(cond2, 'genes', sep = '.'),
                                       'common.genes',
                                       'num.common.genes',
                                       paste(cond1, 'only.genes', sep = '.'),
                                       paste('num', cond1, 'only.genes', sep = '.'),
                                       paste(cond2, 'only.genes', sep = '.'),
                                       paste('num', cond2, 'only.genes', sep = '.'))
  
  
  return(marker1.marker2.lists)
}

prod.desc <- read.csv('../Input/BdivCellCycle/ProductDescription/BDvi_Prod_desc.csv')
prod.desc <- prod.desc %>% transmute(GeneID = Gene.ID, Product.Description = Product.Description)

#saveRDS(prod.desc, '../Input/scClock/prod.desc.RData')

### b. divergense
input.dir.bdiv <- "../Input/scRNAseqBdivCS/"
count.files <- list.files(input.dir.bdiv)
num.total.files <- length(count.files)

file.info <- data.frame(time = rep(NA, num.total.files), 
                        reactivate = rep(NA, num.total.files), 
                        filename = rep(NA, num.total.files))

S.O.list <- list()

for(i in 1:num.total.files){
  file.info$filename[i] <- count.files[i]
  file.info$time[i] <- gsub('OUT', '', gsub('bd', '', strsplit(file.info$filename[i], split = '_')[[1]][1]))
  file.info$reactivate[i] = ifelse(grepl('OUT', file.info$filename[i]), 'Y', 'N')
  cat(paste('processing file', file.info$filename[i]))
  cat('\n')
  L <- processCount(input.dir.bdiv, file.info$filename[i], file.info$time[i], file.info$reactivate[i])
  S.O.list <- c(S.O.list, list(L))
}

file.info$spp <- paste('BDiv', file.info$time, file.info$reactivate, sep = '')

print(file.info)
