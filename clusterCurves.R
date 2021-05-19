library(dtwclust)
library(doParallel)
library(tidyverse)
library(tidytext)

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


## Get the sme mean splines and filter to include marker genes only
sync.tc.mus <- smoothSplineSmeFits(sync.tc.fits, unique(sync.tc.df$variable), extend = T)
colnames(sync.tc.mus) <- c('GeneID', 't', 'y')
sync.tc.mus <- sync.tc.mus %>% dplyr::filter(GeneID %in% BD.markers.sig$GeneID) %>%
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  as.data.frame()

sc.tc.mus <- smoothSplineSmeFits(sc.tc.fits, unique(sc.tc.df$variable), extend = F)
colnames(sc.tc.mus) <- c('GeneID', 't', 'y')
sc.tc.mus <- sc.tc.mus %>% dplyr::filter(GeneID %in% BD.markers.sig$GeneID) %>%
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  as.data.frame()

## Generate the clusters
num.clust <- 4L

sync.hc_dtw <- dtwClustCurves(sync.tc.mus, nclust = num.clust)
plot(sync.hc_dtw, type = 'sc')
plot(sync.hc_dtw, type = "series", clus = 1L)
plot(sync.hc_dtw, type = "centroids", clus = 1L)

saveRDS(sync.hc_dtw, '../Input/scClock/sync.hc_dtw.RData')

sc.hc_dtw <- dtwClustCurves(sc.tc.mus, nclust = num.clust)
plot(sc.hc_dtw, type = 'sc')
plot(sc.hc_dtw, type = 'centroids')
plot(sc.hc_dtw, type = "series", clus = 2L)
plot(sc.hc_dtw, type = "centroids", clus = 2L)

saveRDS(sc.hc_dtw, '../Input/scClock/sc.hc_dtw.RData')

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



## Reorder the genes within each cluster.
hc_eucledian.df <- withinCalssReOrder(sync.tc.mus.scale) 
sync.tc.mus.scale <- left_join(sync.tc.mus.scale, hc_eucledian.df, by = 'GeneID')
sync.tc.mus.scale <- sync.tc.mus.scale %>% mutate(cluster = as.factor(cluster),
                                                  GeneID = reorder_within(GeneID, hc_eucledian.order, cluster)) 


hc_eucledian.df <- withinCalssReOrder(sc.tc.mus.scale) 
sc.tc.mus.scale <- left_join(sc.tc.mus.scale, hc_eucledian.df, by = 'GeneID')
sc.tc.mus.scale <- sc.tc.mus.scale %>% mutate(cluster = as.factor(cluster),
                                                  GeneID = reorder_within(GeneID, hc_eucledian.order, cluster)) 

## map the clusteres
sync.overlap <- matchClustersToPhase(sync.hc_dtw.df, BD.markers.sig)
sc.overlap <- matchClustersToPhase(sc.hc_dtw.df, BD.markers.sig)

sync.phase.match <- sync.overlap %>% group_by(cluster) %>% summarize(phase = markers[which.max(percent)])
sync.tc.mus.scale <- left_join(sync.tc.mus.scale, sync.phase.match, by = 'cluster')

p1 <- ggplot(sync.tc.mus.scale, aes(x = t, y = reorder_within(GeneID, order, phase), fill = y)) + 
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
  facet_grid(phase~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE)) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p1)  


ggsave(filename="../Output/scClockFigs/curve_cluster_heatmap_sync.pdf",
       plot=p1,
       width = 5, height = 5,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

sc.overlap <- matchClustersToPhase(sc.hc_dtw.df, BD.markers.sig)

sc.phase.match <- sc.overlap %>% group_by(cluster) %>% summarize(phase = markers[which.max(percent)])
sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.phase.match, by = 'cluster')

p2 <- ggplot(sc.tc.mus.scale, aes(x = t, y = reorder_within(GeneID, order, phase), fill = y)) + 
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
  facet_grid(phase~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE)) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p2)  


ggsave(filename="../Output/scClockFigs/curve_cluster_heatmap_sc.pdf",
       plot=p2,
       width = 5, height = 5,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

p3 <- ggplot(sc.overlap, aes(x = cluster, y = markers, fill = percent)) + 
  geom_tile() + 
  #scale_x_discrete(expand=c(0,0)) +
  ylab("markers") + xlab("clusters") + theme_bw() + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(1000)) +
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

