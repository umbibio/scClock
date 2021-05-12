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

getSlingShot <- function(S.O, method = 'pca'){
  sds <- slingshot(Embeddings(S.O, method), 
                   clusterLabels = S.O$seurat_clusters, 
                   start.clus = 1, stretch = 2)
  
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


plot.sme <-function(fit, v){
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
  upper.band <- spline(x = as.numeric(colnames(fit$coefficients)), 
                       y = fit$coefficients[1, ] + 1.96 * sqrt(mu.variance), 
                       method = "natural", n = 500)
  lower.band <- spline(x = as.numeric(colnames(fit$coefficients)), 
                       y = fit$coefficients[1, ] - 1.96 * sqrt(mu.variance), 
                       method = "natural", n = 500)
  ylim <- range(ylim, upper.band$y, lower.band$y)
  
  plot(x = fit$data$tme, y = fit$data$y, ylim = ylim, xlim = xlim, xaxt="none",
       xlab = '', ylab = '', col = 'black', cex = 1.2, main = '',
       cex.lab = 1.2, font = 2)
  grid()

  for (i in 1:length(fs)) {
    lines(fs[[i]], lty = "dashed", col = 'black', lwd = 0.8)
  }
  
  lines(mu, lwd = 2, col = 'red')
  
  col.meanCurve.rgb <- col2rgb('red')
  polygon(x = c(upper.band$x, rev(lower.band$x)), 
          y = c(upper.band$y,rev(lower.band$y)), 
          col = rgb(col.meanCurve.rgb[1], 
                    col.meanCurve.rgb[2], col.meanCurve.rgb[3], alpha = 125, 
                    maxColorValue = 255), border = NA)

  axis(1, seq(min(fit$data$tme), max(fit$data$tme), length = 13),
       labels = seq(0, 12), font=2)
  mtext(side=1, line=2, "Time (h)", col="black", font=2,cex=1.1)
  mtext(side=2, line=2, "log2(expr)", col="black", font=2,cex=1.1)
  title(main = v , cex.lab = 1.2, line = 0.5)
}

