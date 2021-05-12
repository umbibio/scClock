library(dtwclust)
library(doParallel)
library(tidyverse)

#### New method for fitting
priorCurveClass <- function(v, s.f, cell.cycle.int){
  
  ## Fitting the estimated kernel with smooth splines
  spline.fit <- smooth.spline(x = as.numeric(colnames(s.f$coefficients)), 
                              y = s.f$coefficients[1,])
  
  ## Compute the derivatives of the fitted splines
  s.0 <- predict(spline.fit, spline.fit$x, deriv=0)
  s.1 <- predict(spline.fit, spline.fit$x, deriv=1)
  s.derv <- data.frame(s0=s.0$y, s1=s.1$y)
  
  ## Get the location of the exterema
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
  
  ##take the max peak
  e <- entity.x[which.max(entity.y)]
  tm.idx <- cell.cycle.int$time.idx[which(cell.cycle.int$time.idx - e >=0)[1]]
  cc <- paste(unique(cell.cycle.int$cell.cycle[cell.cycle.int$time.idx == tm.idx]),
              collapse = '_')
  sub.cc <- paste(unique(cell.cycle.int$sub.cell.cycle[cell.cycle.int$time.idx == tm.idx]),
                  collapse = '_')
  
  L <- data.frame(GeneID = v,
                  entity.x = entity.x[which.max(entity.y)],
                  entity.y = entity.y[which.max(entity.y)],
                  class1 = cc,
                  class2 = sub.cc)
  
  
  return(L)
}

classifyCurves <- function(v, s.f, cell.cycle.int){
  
  ## Fitting the estimated kernel with smooth splines
  spline.fit <- smooth.spline(x = as.numeric(colnames(s.f$coefficients)), 
                              y = s.f$coefficients[1,])
  
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
  
  
  ind.e <- {}
  for(e in entity.x){
    ind.e <- c(ind.e, which(cell.cycle.int$time.idx - e >=0)[1])
  }
  
  if(length(ind.e) > 0){
    L <- data.frame(GeneID = rep(v, length(ind.e)),
                    entity.x = entity.x,
                    entity.y = entity.y,
                    class1 = cell.cycle.int$cell.cycle[ind.e],
                    class2 = cell.cycle.int$sub.cell.cycle[ind.e])
  }else{
    L <- data.frame(GeneID = v,
                    entity.x = NA,
                    entity.y = NA,
                    class1 = NA,
                    class2 = NA)
  }
  
  return(L)
}

dtwClustCurves <- function(tc.mus, nclust = 6L){
  ## Calculate clusters in parallel
  num.cores <- detectCores(all.tests = FALSE, logical = TRUE)
  
  workers <- makeCluster(num.cores)
  invisible(clusterEvalQ(workers, library("dtwclust")))
  registerDoParallel(workers)
  tc.mus <- lapply(2:ncol(tc.mus), function(i) c(as.numeric(sync.tc.mus[,i])))
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


## Get the sme mean splies
sync.tc.mus <- smoothSplineSmeFits(sync.tc.fits, unique(sync.tc.df$variable), extend = T)
colnames(sync.tc.mus) <- c('GeneID', 't', 'y')
sync.tc.mus <- sync.tc.mus %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  as.data.frame()

sc.tc.mus <- smoothSplineSmeFits(sc.tc.fits, unique(sc.tc.df$variable), extend = F)
colnames(sc.tc.mus) <- c('GeneID', 't', 'y')
sc.tc.mus <- sc.tc.mus %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  as.data.frame()

## Generate the clustes
num.clust <- 8L

sync.hc_dtw <- dtwClustCurves(sync.tc.mus, nclust = num.clust)
plot(sync.hc_dtw, type = 'sc')
plot(sync.hc_dtw, type = "series", clus = 1L)
plot(sync.hc_dtw, type = "centroids", clus = 1L)

sc.hc_dtw <- dtwClustCurves(sc.tc.mus, nclust = num.clust)
plot(sc.hc_dtw, type = 'sc')
plot(sc.hc_dtw, type = "series", clus = 1L)
plot(sc.hc_dtw, type = "centroids", clus = 1L)



## Scale the data for heatmaps
sync.tc.mus.scale <- sync.tc.mus
sync.tc.mus.scale[,2:ncol(sync.tc.mus.scale)] <- scale(sync.tc.mus.scale[,2:ncol(sync.tc.mus.scale)],
                                                       center = T,scale = T)
sync.tc.mus.scale <- sync.tc.mus.scale %>%  as.data.frame() %>% 
  pivot_longer(!starts_with('t'), names_to = 'GeneID', values_to = 'y')


sc.tc.mus.scale <- sc.tc.mus
sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)] <- scale(sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)],
                                                       center = T,scale = T)
sc.tc.mus.scale <- sc.tc.mus.scale %>%  as.data.frame() %>% 
  pivot_longer(!starts_with('t'), names_to = 'GeneID', values_to = 'y')


## Add curve cluster info
sync.hc_dtw.df <- data.frame(GeneID = unique(sync.tc.mus.scale$GeneID), 
                        order = as.numeric(sync.hc_dtw$order),
                        cluster = cutree(sync.hc_dtw,k = num.clust))
sync.tc.mus.scale <- left_join(sync.tc.mus.scale, sync.hc_dtw.df, by = 'GeneID')

sync.tc.mus.scale$GeneID <- factor(as.character(sync.tc.mus.scale$GeneID),
                                   levels = unique(sync.tc.mus.scale$GeneID[sync.tc.mus.scale$order]))


sc.hc_dtw.df <- data.frame(GeneID = unique(sc.tc.mus.scale$GeneID), 
                             order = as.numeric(sc.hc_dtw$order),
                             cluster = cutree(sc.hc_dtw,k = num.clust))

sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.hc_dtw.df, by = 'GeneID')

sc.tc.mus.scale$GeneID <- factor(as.character(sc.tc.mus.scale$GeneID),
                                   levels = unique(sc.tc.mus.scale$GeneID[sc.tc.mus.scale$order]))



## Must order the clusters and curves within cluster for better heatmap

p <- ggplot(sync.tc.mus.scale, aes(x = t, y = GeneID, fill = y)) + 
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
  facet_grid(cluster~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE)) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p)  

