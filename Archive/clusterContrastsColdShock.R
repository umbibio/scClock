## Merge the data 
## Compare 0h vs. 36h and 36h out
L <- processeMerged(S.O.list, file.info, c(1, 3, 4))

all.samples.integrated_0hr_36hrN_36hrY <- L$all.samples.integrated
all.spp.list_0hr_36hrN_36hrY <- L$all.spp.list

#saveRDS(all.samples.integrated_0hr_36hrN_36hrY, '../Input/scClock/all.samples.integrated_0hr_36hrN_36hrY.RData')
#saveRDS(all.spp.list_0hr_36hrN_36hrY, '../Input/scClock/all.spp.list_0hr_36hrN_36hrY.RData')



p1 <- DimPlot(all.samples.integrated_0hr_36hrN_36hrY, reduction = "pca", 
              #group.by = "cells", 
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


pc <- getPCA(all.samples.integrated_0hr_7dN_7dY)
pc$spp <- all.samples.integrated_0hr_7dN_7dY$spp
pc$orig.ident <- all.samples.integrated_0hr_7dN_7dY$orig.ident
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
p <- ggplot(cluster.sim.marker2.marker2 , aes(x = `7dN1`, y = `7dN2`, fill = num.comm.genes)) + 
  geom_tile() + theme_bw() + 
  geom_text(aes(label = num.comm.genes, size=6)) +
  scale_y_discrete(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  ylab("7dN2") + xlab("7dN1") + 
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
