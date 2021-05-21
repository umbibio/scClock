## functions
##
prep_S.O <- function(S.O, res = 0.1){
  S.O <- NormalizeData(S.O, normalization.method = "LogNormalize", scale.factor = 10000)
  S.O <- FindVariableFeatures(S.O, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(S.O)
  S.O <- ScaleData(S.O, features = all.genes)
  S.O <- RunPCA(S.O, features = VariableFeatures(object = S.O))
  S.O <- FindNeighbors(S.O, dims = 1:13)
  #S.O <- FindClusters(S.O, resolution = 0.1)
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
                                           ifelse(expr <= quantile(expr)[4], 3,4))))
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
