library(dtwclust)
library(doParallel)
library(tidyverse)
library(tidytext)


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

