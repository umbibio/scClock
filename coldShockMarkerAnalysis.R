library(ggplot2)
library(cowplot)
library(patchwork)
library(tidytext)
library(ggVennDiagram)
library(UpSetR)


makeGlobalContrasts <- function(all.samples.integrated, c.mode  = c('full', 'WT')){
  
  c.mode <- match.arg(c.mode)
  
  objs <- as.character(unique(all.samples.integrated@meta.data$spp))
  
  contrasts <- data.frame(ref = objs, dummy = 1) %>% full_join( data.frame(query = objs, dummy = 1), by = 'dummy') %>% 
    mutate(ref.time = as.numeric(gsub("hr.*", "", gsub("7d", "168hr", gsub("BDiv", "", ref)))),
           query.time = as.numeric(gsub("hr.*", "", gsub("7d", "168hr", gsub("BDiv", "", query)))),
           ref.reactive = ifelse(grepl('Y', ref), 'Y', 'N'),
           query.reactive = ifelse(grepl('Y', query), 'Y', 'N'))
  if(c.mode == 'full'){
    my.contrasts <- contrasts %>%
      dplyr::filter( 
        ( (ref.time < query.time) ) | 
          ( (ref.time == query.time) & (ref.reactive != query.reactive) ) 
      ) %>% dplyr::filter(!(ref.reactive == 'N' & query.reactive == 'Y' & ref.time != 0)) %>%
      dplyr::filter(!(ref.reactive == 'Y' & query.reactive == 'Y')) %>% 
      dplyr::filter(!(ref.reactive == 'Y' & query.reactive == 'N' & ref.time != query.time ))
  }else if(c.mode == 'WT'){
    my.contrasts <- contrasts %>% dplyr::filter(ref.time == 0 & query.reactive != 'Y' & ref.time < query.time)
  }
  
  return(my.contrasts)
  
}

makeMatchedContrasts <- function(all.samples.integrated, c.mode  = c('full', 'WT')){
  
  c.mode <- match.arg(c.mode)
  tmp <- makeGlobalContrasts(all.samples.integrated, c.mode  = c.mode)
 
  clusters <- unique(all.samples.integrated@meta.data$seurat_clusters)
 
  
  my.contrasts <- data.frame(ref = paste(rep(tmp$ref, each = length(clusters)), rep(sort(clusters), nrow(tmp)), sep = '_'),
                             query = paste(rep(tmp$query, each = length(clusters)), rep(sort(clusters), nrow(tmp)), sep = '_'),
                             cluster = rep(sort(clusters), length(tmp$ref))) 
  
  return(my.contrasts)
  
}

## Integrate Samples
print(file.info)
ref.ind <- c(1,5,2,3,6,4,7) 
all.samples.integrated <- processeMergedS.O(S.O.list, file.info, ref.ind, res = 0.1, SC = F)
saveRDS(all.samples.integrated, '../Input/scClock/coldShock/all_samples_integrated.RData')

all.samples.integrated@meta.data$spp <- factor(all.samples.integrated@meta.data$spp, 
                                               levels = c("BDiv0hrN", "BDiv4hrN", "BDiv12hrN", 
                                                          "BDiv36hrN", "BDiv7dN", 
                                                          "BDiv36hrY", "BDiv7dY"))


#Idents(all.samples.integrated) <- "phase.cond"
p1 <- DimPlot(all.samples.integrated, reduction = "umap", 
              split.by = 'spp',
              pt.size = 1,
              #shape.by='spp',
              label = TRUE, label.size = 6) + NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p1)


# Differential Expression
# Identify global markers independent of cell cycle phase. Are there any global Cold Shock regulators?
# Two class comparison: fc = ident.1/ident.2

DefaultAssay(all.samples.integrated) <- "RNA"
Idents(all.samples.integrated) <- "spp"



c.mode <- 'WT'
contrasts <- makeGlobalContrasts(all.samples.integrated, 'WT')
global.DEGs <- mclapply(split(contrasts, seq(nrow(contrasts))), function(x){
  tmp <- FindMarkers(all.samples.integrated, ident.1 = x$query, ident.2 = x$ref, verbose = T)
  tmp$ref <- x$ref
  tmp$query <- x$query
  tmp$gene <- rownames(tmp)
  tmp$GeneID <- gsub('-', '_', tmp$gene)
  return(tmp)
})

global.shock.markers <- bind_rows(global.DEGs)
saveRDS(global.shock.markers, '../Input/scClock/coldShock/global_shock_markers.RData')

global.shock.markers.sig <- global.shock.markers %>% 
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > log2(1.5)) %>% arrange(desc(abs(avg_log2FC)))

global.shock.markers.overlap <- global.shock.markers.sig %>% group_by(GeneID) %>% 
  summarise(num.contrasts = length(query)) %>% arrange(desc(num.contrasts))

global.shock.markers.sig <- left_join(global.shock.markers.sig, global.shock.markers.overlap, by = 'GeneID')
global.shock.markers.stats <- global.shock.markers.sig %>%
  mutate(reg = ifelse(avg_log2FC > 0, 1, -1)) %>%
  group_by(ref, query, reg) %>% summarise(num.DEGs = n()) 

print(global.shock.markers.stats)

global.shock.markers.sig <- left_join(global.shock.markers.sig, prod.desc, by = 'GeneID')
write.xlsx(global.shock.markers.sig, '../Output/scClockOut/global_shock_markers_WT_base.xlsx')

tmp <- global.shock.markers.sig %>% dplyr::filter(avg_log2FC < 0) %>% group_by(ref, query) %>%
  summarise(genes = list(GeneID))

venn.list <- tmp$genes
names(venn.list) <- paste(tmp$ref, "_vs_", tmp$query, sep = '')

ggVennDiagram(venn.list)
upset(fromList(venn.list))


DefaultAssay(all.samples.integrated) <- 'RNA'
Idents(all.samples.integrated) <- 'phase.cond' 
p2 <- FeaturePlot(object = all.samples.integrated, 
                  shape.by = 'spp',
                  split.by = 'spp',
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = "Bdiv-033980",
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p2)



p3 <- DotPlot(all.samples.integrated, features = 'Bdiv-038820', 
              cols = c("blue", "red", "green", "yellow"),
              dot.scale = 8) + RotatedAxis()

plot(p3)

Idents(all.samples.integrated) <- 'spp'
DefaultAssay(all.samples.integrated) <- 'RNA'
VlnPlot(object = all.samples.integrated, features = 'Bdiv-033980')


### Differential expression analysis
## Cell cycle phase specific

all.spp.list <- SplitObject(all.samples.integrated, split.by = "spp")

## Re-normalize splitted data
for (i in 1:length(all.spp.list)) {
  Idents(all.spp.list[[i]]) <- 'phase.cond'
  all.spp.list[[i]] <- NormalizeData(all.spp.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
}


contrasts <- makeMatchedContrasts(all.samples.integrated, c.mode  = 'WT')

GM <- getCellCyclePhaseMarkers(all.spp.list)

## For each cluster, pool global markers.
pooled_markers <- bind_rows(GM$all.markers.list.sig) %>% group_by(glob.clust) %>%
  summarise(glob.markers = list(unique(gene)), num.markers = length(unique(gene)))

#ident.1 case, ident.2 is control
Idents(all.samples.integrated) <- 'phase.cond'
DefaultAssay(all.samples.integrated) <- 'RNA'
matched.DEGs <- mclapply(split(contrasts, seq(nrow(contrasts))), function(x){
  tmp <- FindMarkers(all.samples.integrated, ident.1 = x$query, ident.2 = x$ref, verbose = T)
  ind <- rownames(tmp) %in% unlist(pooled_markers$glob.markers[which(pooled_markers$glob.clust == x$cluster)]) 
  tmp$ref <- x$ref
  tmp$query <- x$query
  tmp$cluster <- x$cluster
  tmp$gene <- rownames(tmp)
  tmp$GeneID <- gsub('-', '_', tmp$gene)
  return(tmp[ind, ])
  #return(tmp)
})


matched.DEGs <- bind_rows(matched.DEGs)
saveRDS(matched.DEGs, '../Input/scClock/coldShock/matched_DEGs.Rdata')

matched.DEGs.sig <- matched.DEGs %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > log2(1.5)) %>% arrange(desc(abs(avg_log2FC)))

matched.DEGs.overlap <- matched.DEGs.sig %>% group_by(GeneID) %>% 
  summarise(num.contrasts = length(query)) %>% arrange(desc(num.contrasts))

matched.DEGs.sig <- left_join(matched.DEGs.sig, matched.DEGs.overlap, by = 'GeneID')

matched.DEGs.stats <- matched.DEGs.sig %>%
  mutate(reg = ifelse(avg_log2FC > 0, 1, -1)) %>%
  group_by(query, reg) %>% summarise(num.DEGs = n()) 

print(matched.DEGs.stats)

matched.DEGs.stats$cluster <- as.factor(as.numeric(gsub('.*_', '', matched.DEGs.stats$query)))


matched.DEGs.top <- matched.DEGs.sig %>% dplyr::filter(avg_log2FC > 0) %>% group_by(cluster) %>% top_n(1, abs(avg_log2FC))



p2 <- FeaturePlot(object = all.samples.integrated, 
                  shape.by = 'spp',
                  split.by = 'spp',
                  label = T, pt.size = 0.6, label.size = 3, 
                  features = matched.DEGs.top$gene,
                  cols = c("lightgrey", "red"), reduction = "pca") 

plot(p2)



VlnPlot(object = all.samples.integrated, features = 'Bdiv-015710')


tmp <- left_join(matched.DEGs.sig, prod.desc, by = 'GeneID')
write.xlsx(tmp, '../Output/scClockOut/matched_shock_markers.xlsx')


tmp <- matched.DEGs.sig %>% dplyr::filter(avg_log2FC > 0) %>% group_by(ref, query) %>%
  summarise(genes = list(GeneID))

venn.list <- tmp$genes
names(venn.list) <- paste(tmp$ref, "_vs_", tmp$query, sep = '')

ggVennDiagram(venn.list)
upset(fromList(venn.list))

