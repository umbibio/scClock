## functions
##
prep_S.O <- function(S.O, res = 0.1){
  S.O <- NormalizeData(S.O, normalization.method = "LogNormalize", scale.factor = 10000)
  S.O <- FindVariableFeatures(S.O, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(S.O)
  S.O <- ScaleData(S.O, features = all.genes)
  S.O <- RunPCA(S.O, features = VariableFeatures(object = S.O))
  S.O <- FindNeighbors(S.O, dims = 1:13)
  S.O <- FindClusters(S.O, resolution = res)
  S.O <- RunUMAP(S.O, dims = 1:13)
  return(S.O)
}


getNormExpr <- function(S.O){
  expr.norm <- as.matrix(S.O[["RNA"]]@data)
  rownames(expr.norm) <- gsub('-', '_',rownames(expr.norm))
  expr.norm <- as.data.frame(expr.norm) 
  expr.norm <- expr.norm %>% dplyr::mutate(GeneID = rownames(expr.norm))
  expr.norm <- expr.norm %>% gather(key = Sample, value = expr, -GeneID)
  expr.norm <- expr.norm  %>% group_by(GeneID) %>% 
    dplyr::mutate(quantile = ifelse(expr <= quantile(expr)[2], 1, 
                                    ifelse(expr <= quantile(expr)[3], 2, 
                                           ifelse(expr <= quantile(expr)[4], 3,4)))) %>% ungroup()
  return(expr.norm)
}


getUmap <- function(S.O){
  umap <- S.O[['umap']]@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) 
  #umap <- left_join(umap, S.O@meta.data, by = 'Sample')
  umap$cluster <- S.O$seurat_clusters
  
  return(umap)
}

getPCA <- function(S.O){
  pc <- S.O[['pca']]@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) 
  #pc <- left_join(pc, S.O@meta.data, by = 'Sample')
  pc$cluster <- S.O$seurat_clusters
  
  return(pc)
}

getPrinCurve <- function(pc.db){
  x <- as.matrix(pc.bd[,c(1,2)])
  fit <- principal_curve(x)
  pt <- fit$lambda
  
  ## reversing the order of time
  pt <- max(pt) - pt
  pt <- (pt - min(pt))/(max(pt) - min(pt))
  cell.ord <- fit$ord[seq(length(fit$ord), 1, by = -1)]
  
  
  s.data <- data.frame(Sample = rownames(fit$s), 
                       cell.ord = cell.ord,
                       pt = pt,
                       sc1 = fit$s[,1], 
                       sc2 = fit$s[,2])
  
  return(s.data)
}

getSlingShot <- function(S.O, method = 'pca'){
  sds <- slingshot(Embeddings(S.O, method), 
                   clusterLabels = S.O$seurat_clusters, 
                   start.clus = 0, end.clus = 3, stretch = 2)
  
  pt <- slingPseudotime(sds)
  
  s.data <- data.frame(cell.ord = sds@curves$curve1$ord, pt = pt, sds@curves$curve1$s[,1:2])
  s.data$Sample <- rownames(s.data)
  s.data <- s.data %>% transmute(Sample = Sample, cell.ord = cell.ord, 
                                 pt = curve1, sc1 = PC_1, sc2 = PC_2)
  
  #s.data <- left_join(s.data, S.O.combined@meta.data, by = 'Sample')
  
  return(s.data)
}


cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}


fitSingleSme <-function(my.tc, v){
  tryCatch(
    expr = {
      fit <- sme(my.tc[my.tc$variable==v,c("y","tme","ind")], criteria = 'AIC')
      return(fit)
    },
    error = function(v){
      message(paste('error:', v))
    },
    warning = function(w){
      message(paste('error:', v))
    },
    finally = {
      fit <- sme(my.tc[my.tc$variable==v,c("y","tme","ind")])
      return(fit)
    }
  )    
}

## spline the fitted values to get the means
splineSmeFits <- function(fits, variables){
  mus <- lapply(fits, function(x) spline(x = as.numeric(colnames(x$coefficients)),
                                         y = x$coefficients[1,], 
                                         n = as.numeric(colnames(x$coefficients))[ncol(x$coefficients)] - 
                                           as.numeric(colnames(x$coefficients))[1] + 1, 
                                         method = "natural")) 
  mus.y <- unlist(lapply(mus, `[[`, 2))
  mus.x <- unlist(lapply(mus, `[[`, 1))
  lens  <- unlist(lapply(lapply(mus, `[[`, 1), length))
  
  mus <- data.frame(variable = rep(variables, times = lens),
                    tme = mus.x,
                    y = mus.y)
  
  return(mus)
  
}

## spline the fitted values to get the means
smoothSplineSmeFits <- function(fits, variables, extend = F){
  ## Fitting the estimated kernel with smooth splines
  
  
  mus <- lapply(fits, function(x) 
    smooth.spline(x = as.numeric(colnames(x$coefficients)), 
                  y = x$coefficients[1,])) 
  if(extend){
    mus <- lapply(mus, function(x)
      predict(x, seq(0, 12, by = 1/3)))
  }
  
  mus.y <- unlist(lapply(mus, `[[`, 2))
  mus.x <- unlist(lapply(mus, `[[`, 1))
  lens  <- unlist(lapply(lapply(mus, `[[`, 1), length))
  
  mus <- data.frame(variable = rep(variables, times = lens),
                    tme = mus.x,
                    y = mus.y)
  
  return(mus)
  
}


plot.sme <-function(fit, v, conf = T){
  mu <- spline(x = as.numeric(colnames(fit$coefficients)), 
               y = fit$coefficients[1,], n = 500, 
               method = "natural")
  fs <- lapply(2:nrow(fit$coefficients), function(i) {
    spline(x = as.numeric(colnames(fit$coefficients)), 
           y = fit$coefficients[1,] + fit$coefficients[i, ], 
           method = "natural", 
           n = 500)})
  
  
  ylim <- range(fit$data$y, mu$y, sapply(fs, "[[", "y"))
  xlim <- range(as.numeric(colnames(fit$coefficients)))
  mu.variance <- diag(vcov(fit))
  if(conf){
    upper.band <- spline(x = as.numeric(colnames(fit$coefficients)), 
                         y = fit$coefficients[1, ] + 1.96 * sqrt(mu.variance), 
                         method = "natural", n = 500)
    lower.band <- spline(x = as.numeric(colnames(fit$coefficients)), 
                         y = fit$coefficients[1, ] - 1.96 * sqrt(mu.variance), 
                         method = "natural", n = 500)
    ylim <- range(ylim, upper.band$y, lower.band$y)
  }
  
  plot(x = fit$data$tme, y = fit$data$y, ylim = ylim, xlim = xlim, xaxt="none",
       xlab = '', ylab = '', col = 'black', cex = 1.2, main = '',
       cex.lab = 1.2, font = 2)
  grid()

  for (i in 1:length(fs)) {
    lines(fs[[i]], lty = "dashed", col = 'black', lwd = 0.8)
  }
  
  lines(mu, lwd = 2, col = 'red')
  
  col.meanCurve.rgb <- col2rgb('red')
  if(conf){
    polygon(x = c(upper.band$x, rev(lower.band$x)), 
            y = c(upper.band$y,rev(lower.band$y)), 
            col = rgb(col.meanCurve.rgb[1], 
                      col.meanCurve.rgb[2], col.meanCurve.rgb[3], alpha = 125, 
                      maxColorValue = 255), border = NA)
  }
 
  axis(1, seq(min(fit$data$tme), max(fit$data$tme), length = 13),
       labels = seq(0, 12), font=2)
  mtext(side=1, line=2, "Time (h)", col="black", font=2,cex=1.1)
  mtext(side=2, line=2, "log2(expr)", col="black", font=2,cex=1.1)
  title(main = v , cex.lab = 1.2, line = 0.5)
}



dtwClustCurves <- function(tc.mus, nclust = 6L){
  ## Calculate clusters in parallel
  num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
  
  workers <- makeCluster(num.cores)
  invisible(clusterEvalQ(workers, library("dtwclust")))
  registerDoParallel(workers)
  tc.mus <- lapply(2:ncol(tc.mus), function(i) c(as.numeric(tc.mus[,i])))
  hc_dtw <- tsclust(tc.mus, 
                    type = "h", 
                    k = nclust, 
                    distance = "dtw", 
                    control = hierarchical_control(method = "complete"),
                    centroid = shape_extraction, 
                    #preproc = NULL, 
                    preproc = zscore,
                    trace = T,
                    args = tsclust_args(dist = list(window.size = 4L))
  )
  
  stopCluster(workers)
  registerDoSEQ()
  
  return(hc_dtw)
}

withinCalssReOrder <- function(tc.mus){
  clust.ord <- tc.mus %>% dplyr::select(GeneID, cluster) %>% distinct() %>% 
    group_by(cluster) %>% summarise(GeneSet = list(GeneID)) 
  
  clusters <- unique(tc.mus$cluster)
  hc_eucledian.df <- {}
  for(i in 1:length(clusters)){
    class.i <- tc.mus %>% dplyr::filter(cluster == clusters[i]) %>% dplyr::select(GeneID, y, t) %>%
      pivot_wider(names_from = GeneID, values_from = y) %>% as.data.frame()
    
    
    ## Hierarchical clustering with Eucledian distance
    hc_eucledian <- hclust(dist(t(as.matrix(class.i[,-1] ))), method = "ward.D")
    hc_eucledian.df <- rbind(hc_eucledian.df, 
                             data.frame(GeneID = colnames(class.i[,-1]), 
                                        hc_eucledian.order = hc_eucledian$order,
                                        hc_eucledian.cluster = cutree(hc_eucledian,k = 10)))
  }
  
  return(hc_eucledian.df)
}

## Calculate th overlap of clusters with marker genes and re-order clusters accordingly
matchClustersToPhase <- function(hc_dtw.df, markers.sig){
  join.hc_dtw.df <- inner_join(hc_dtw.df, markers.sig, by = 'GeneID')
  totals <- join.hc_dtw.df %>% group_by(cluster.x) %>% summarise(totals = n())
  overlap <- join.hc_dtw.df %>% group_by(cluster.x, cluster.y) %>% summarise(overlap = n())
  overlap <- left_join(overlap, totals, by = 'cluster.x') %>% mutate(precent = overlap/totals)
  
  colnames(overlap) <- c('cluster', 'markers', 'counts', 'totals', 'percent')
  overlap <- overlap %>% dplyr::select(c('cluster','markers', 'percent')) %>% 
    pivot_wider(names_from = 'markers', values_from = 'percent')
  overlap[is.na(overlap)] <- 0
  overlap <- overlap %>% pivot_longer(-c('cluster'), names_to = 'markers', values_to = 'percent')
  
  overlap$markers[overlap$markers == 0] <- 'G1a'
  overlap$markers[overlap$markers == 3] <- 'G1b'
  overlap$markers[overlap$markers == 2] <- 'G1c'
  overlap$markers[overlap$markers == 1] <- 'S/M'
  
  
  overlap$cluster <- factor(overlap$cluster, levels = unique(sort(overlap$cluster)))
  overlap$markers <- factor(overlap$markers, levels = unique(sort(overlap$markers)))
  
  return(overlap)
  
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
    select(all_of(cond1), all_of(cond2) ,'num.comm.genes') %>% distinct()
  
  cluster.sim.marker1.marker2 <- cluster.sim.marker1.marker2 %>% 
    pivot_wider(names_from = !!cond2, values_from = 'num.comm.genes')
  
  cluster.sim.marker1.marker2[is.na(cluster.sim.marker1.marker2)] <- 0
  
  cluster.sim.marker1.marker2 <- cluster.sim.marker1.marker2 %>% 
    pivot_longer(-cond1, names_to = paste0(cond2), values_to = 'num.comm.genes')
  
  return(cluster.sim.marker1.marker2)
}


processCount <- function(input.dir, filename, tt, rr, down.sample = T){
  file.counts <- read.csv(paste(input.dir, filename, sep = ''))
  genes <- file.counts$X
  expr <- file.counts[,-1]
  rownames(expr) <- genes
  
  
  S.O <- CreateSeuratObject(counts = expr, min.cells = 10, min.features = 100)
  
  #VlnPlot(S.O, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
  #FeatureScatter(S.O, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  cutoffs <- quantile(S.O$nCount_RNA, probs = c(0.01, 0.9))
  print(cutoffs)
  S.O <- subset(S.O, subset = nFeature_RNA > cutoffs[1] & nFeature_RNA < cutoffs[2] )
  
  if(down.sample){
    S.O <- subset(x = S.O, downsample = 2000)
  }
  
  names(S.O$orig.ident)
  pheno <- data.frame(Sample = names(S.O$orig.ident))
  spp <- paste('BDiv', tt, rr, sep = '')
  pheno$spp <- spp
  pheno$time <- tt
  pheno$reactivate <- rr
  
  pheno$NAME <- paste(pheno$spp, pheno$Sample, sep = '_')
  
  
  L <- list(pheno = pheno, S.O = S.O)
  
  return(L)
}

mergeS.O <- function(L){
  num.objs <- length(L)
  
  phenos <- lapply(L, `[[`, 1)
  S.Os <-  lapply(L, `[[`, 2)
  
  all.samples <- bind_rows(phenos)
  
  S.O.merge <- merge(S.Os[[1]], y = S.Os[2:num.objs ], add.cell.ids = unique(all.samples$spp))
  
  rownames(all.samples) <- all.samples$NAME
  S.O.merge <- AddMetaData(S.O.merge, metadata = all.samples)  
  
  return(S.O.merge)
}

processeMergedS.O <- function(S.O.list, file.info, ref.ind, res = 0.2, SC = FALSE){
  ## Merge the S.O, add metadata, and re-split by spp and update the S.O.list
  S.O.merge <- mergeS.O(S.O.list[ref.ind])
  S.O.list <- SplitObject(S.O.merge, split.by = "spp")
  
  if(SC){
    S.O.list <- lapply(X = S.O.list, FUN = SCTransform)
    features <- SelectIntegrationFeatures(object.list = S.O.list, nfeatures = 3000)
    S.O.list <- PrepSCTIntegration(object.list = S.O.list, anchor.features = features)
    anchors <- FindIntegrationAnchors(object.list = S.O.list, normalization.method = "SCT", 
                                      anchor.features = features)
    S.O.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  }else{
    S.O.list <- lapply(X = S.O.list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })
    features <- SelectIntegrationFeatures(object.list = S.O.list)
    anchors <- FindIntegrationAnchors(object.list = S.O.list, anchor.features = features)
    S.O.integrated <- IntegrateData(anchorset = anchors)
  }
  
  
  # switch to integrated assay. Make sure to set to RNA for Differential Expression
  DefaultAssay(S.O.integrated) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  S.O.integrated <- ScaleData(S.O.integrated, verbose = FALSE)
  S.O.integrated <- RunPCA(S.O.integrated, npcs = 30, verbose = FALSE)
  S.O.integrated <- RunUMAP(S.O.integrated, reduction = "pca", dims = 1:30)
  S.O.integrated <- FindNeighbors(S.O.integrated, reduction = "pca", dims = 1:30)
  S.O.integrated <- FindClusters(S.O.integrated, resolution = res)
  
  S.O.integrated$phase.cond <- paste(S.O.integrated@meta.data$spp, 
                                     Idents(S.O.integrated), sep = "_")
  Idents(S.O.integrated) <- "phase.cond"
  
  return(S.O.integrated)
}


getCellCyclePhaseMarkers <- function(all.spp.list){
  all.markers.list <- mclapply(all.spp.list, function(x) FindAllMarkers(object = x, only.pos = TRUE))
  all.markers.list <- lapply(all.markers.list, function(x) {
    x$GeneID = gsub('-', '_', x$gene)
    x$glob.clust <- gsub('.*_', '', x$cluster)
    return(x)
  })
  
  
  all.markers.list.sig <- lapply(all.markers.list, function(x) {
    sig.marker = x %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)
    return(sig.marker)
  })
  
  L <- list(all.markers.list = all.markers.list, 
            all.markers.list.sig = all.markers.list.sig)
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
