library(ComplexHeatmap)
library(circlize)


## Get the average expression per timepoint
Idents(all.samples.integrated) <- 'spp'
global.avg.expr <- as.data.frame(log1p(AverageExpression(all.samples.integrated, verbose = FALSE)$RNA))
global.avg.expr$gene <- rownames(global.avg.expr)
global.avg.expr$GeneID <- gsub('-', '_', global.avg.expr$gene)

global.up.reg <- global.shock.markers.sig$GeneID[global.shock.markers.sig$avg_log2FC > 0]
global.down.reg <- global.shock.markers.sig$GeneID[global.shock.markers.sig$avg_log2FC < 0]

global.avg.expr.up <- global.avg.expr[global.avg.expr$GeneID %in% global.up.reg, ]
global.avg.expr.down <- global.avg.expr[global.avg.expr$GeneID %in% global.down.reg, ]


global.mat_up <- global.avg.expr.up[,1:(ncol(global.avg.expr.up) - 2)] %>% dplyr::select(!contains('Y'))
global.mat_down <- global.avg.expr.down[,1:(ncol(global.avg.expr.down) - 2)] %>% dplyr::select(!contains('Y'))


L1 <- makeGlobalHeatMap(global.mat_up, prod.desc, 7)
L2 <- makeGlobalHeatMap(global.mat_down, prod.desc, 7)

draw(L1$ht_list)
draw(L2$ht_list)

write.xlsx(L1$heat.clust, '../Output/scClockOut/global_heatmap_clusters_up.xlsx')
write.xlsx(L2$heat.clust, '../Output/scClockOut/global_heatmap_clusters_down.xlsx')



### Mathced
## Get the average expression per timepoint
Idents(all.samples.integrated) <- 'phase.cond'
matched.avg.expr <- as.data.frame(log1p(AverageExpression(all.samples.integrated, verbose = FALSE)$RNA))
matched.avg.expr$gene <- rownames(matched.avg.expr)
matched.avg.expr$GeneID <- gsub('-', '_', matched.avg.expr$gene)

matched.up.reg <- matched.DEGs.sig$GeneID[matched.DEGs.sig$avg_log2FC > 0]
matched.down.reg <- matched.DEGs.sig$GeneID[matched.DEGs.sig$avg_log2FC < 0]

matched.avg.expr.up <- matched.avg.expr[matched.avg.expr$GeneID %in% matched.up.reg, ]
matched.avg.expr.down <- matched.avg.expr[matched.avg.expr$GeneID %in% matched.down.reg, ]


matched.mat_up <- matched.avg.expr.up[,1:(ncol(matched.avg.expr.up) - 2)] %>% dplyr::select(!contains('Y'))
matched.mat_down <- matched.avg.expr.down[,1:(ncol(matched.avg.expr.down) - 2)] %>% dplyr::select(!contains('Y'))


L1 <- makeMatchedHeatMap(matched.mat_up, prod.desc, 7)
L2 <- makeMatchedHeatMap(matched.mat_down, prod.desc, 7)

draw(L1$ht_list)
draw(L2$ht_list)

write.xlsx(L1$heat.clust, '../Output/scClockOut/global_heatmap_clusters_up.xlsx')
write.xlsx(L2$heat.clust, '../Output/scClockOut/global_heatmap_clusters_down.xlsx')


