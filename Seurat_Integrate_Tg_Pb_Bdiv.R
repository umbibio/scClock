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
library(patchwork)


## Reading the data
## These files were prepared using the python script `get_h5ad_data_kz.ipynb`. 


input.dir <- '../Input/scRNAseqToxoAtlas/kz/'

rh384.expr.file <- 'rh384_expression_filtered.csv'
rh384.obs.file <- 'rh384_obs_filtered.csv'
rh384.var.file <- 'rh384_var_filtered.csv'


input.dir.10x <- "../Input/MalariaCellAtlas/Expression_Matrices/10X/pb10xIDC/"
plasmodium.10x.count.file <- 'pb10xIDC_counts.csv'
plasmodium.10x.pheno.file <- 'pb10xIDC_pheno.csv'

### bdivergense
input.dir.bdiv <- "../Input/scRNAseqBdiv/"
bdiv.count.file <- "bdiv.expr.csv"

## Gene product description
prod.disc <- read.xlsx('../Input/orthologs/ProductDescription.xlsx')
Ribosomal.RNA <- prod.disc$GeneID[grep('ribosom*',prod.disc$ProductDescription)]

## Reciprocal orthologs
GT1.Pberghei <- read.xlsx('../Input/orthologs/rec_GT1.vs.P.bergei.xlsx')
GT1.bdiv <- read.xlsx("../Input/orthologs/rec_GT1.vs.B.divergence.xlsx")
#bdiv.Pberghei <- read.xlsx("../Input/orthologs/rec_B.div_P.berg.xlsx")

GT1.PB <- GT1.Pberghei %>% dplyr::select('query_id', contains("subject_id"))
GT1.bdiv <- GT1.bdiv %>% dplyr::select('query_id', contains('subject_id'))

GT1.pb.bd <- inner_join(GT1.Pberghei, GT1.bdiv, by = 'query_id')
GT1.pb.bd <- GT1.pb.bd %>% dplyr::select('query_id' , contains('subject_id'))


colnames(GT1.pb.bd) <- c('GT1', 'PB', 'bdiv')
colnames(GT1.PB) <- c('GT1', 'PB')
colnames(GT1.bdiv) <- c('GT1', 'bdiv')

## Filter out Ribosomal RNAs
#GT1.pb.bd <- GT1.pb.bd %>% dplyr::filter(!(GT1 %in% Ribosomal.RNA))

#GT1.pb.bd <- left_join(GT1.pb.bd, prod.disc, by = c('GT1' = 'GeneID'))
#saveRDS(GT1.pb.bd, file = 'GT1_PB_BD_ortho.rds')

## Reading in the Toxo data
rh384.expr <- read.csv(paste(input.dir, rh384.expr.file, sep = ''))
rh384.var <- read.csv(paste(input.dir, rh384.var.file, sep = ''))
rh384.obs <- read.csv(paste(input.dir, rh384.obs.file, sep = ''))

## Reading Plasmodium data 10x pberghei
plasmodium.10x.count <- read.csv(paste(input.dir.10x, plasmodium.10x.count.file, sep = ''))
plasmodium.10x.pheno <- read.csv(paste(input.dir.10x, plasmodium.10x.pheno.file, sep = ''))


## Reading bdivergense data
bdiv.count <- read.csv(paste(input.dir.bdiv, bdiv.count.file, sep = ''))


## Processing the RH data
processCounts <- function(expr){
  cols <- expr[,1]
  expr <- t(expr[,-1])
  colnames(expr) <- cols
  return(expr)
}

rh384.expr <- processCounts(rh384.expr)

genes <- bdiv.count$X
bd.expr <- bdiv.count[,-1]
rownames(bd.expr) <- genes

genes.10x <- plasmodium.10x.count$X
pb.10x.expr <- plasmodium.10x.count[,-1]
rownames(pb.10x.expr) <- genes.10x

## sort out species names and conditions

rh384.obs$spp <- 'RH'

rh384.obs$cells <- paste('RH', gsub("\"", "",  as.character(rh384.obs$cell_cycle)), sep = '.')

rh384.obs$NAME <- paste(rh384.obs$cells, rh384.obs$index, sep = '_')

bdiv.pheno <- data.frame(X = colnames(bd.expr))
bdiv.pheno$spp <- 'BDiv' 


S.O.bd <- CreateSeuratObject(counts = bd.expr)
S.O.bd <- prep_S.O(S.O.bd)

bd.clust.info <- data.frame(X=as.character(names(S.O.bd$orig.ident)),
                            cluster=as.character(S.O.bd$seurat_clusters))
bdiv.pheno <- inner_join(bdiv.pheno, bd.clust.info, by = 'X')
bdiv.pheno$cells <- paste('BDiv', bdiv.pheno$cluster, sep = '')

bdiv.pheno$NAME <- paste(bdiv.pheno$cells, 
                         bdiv.pheno$X, sep = '_')



# Figure 3B
plasmodium.10x.pheno <- plasmodium.10x.pheno %>% 
  mutate(cells = case_when(absclust == 0 ~ "TrpE",
                           absclust == 1 ~ "TrpM",
                           absclust == 2 ~ "RngL",
                           absclust == 3 ~ "TrpL",
                           absclust == 4 ~ "SchE",
                           absclust == 5 ~ "SchL",
                           absclust == 6 ~ "RngE",
                           absclust == 7 ~ "SchM"))



plasmodium.10x.pheno$spp <- 'PBer'
plasmodium.10x.pheno$NAME <- paste(plasmodium.10x.pheno$cells, plasmodium.10x.pheno$X, sep = '_')

rh384.ind <- match(colnames(rh384.expr), rh384.obs$index)
pb.10x.ind <- match(colnames(pb.10x.expr), plasmodium.10x.pheno$X)
bd.ind <- match(colnames(bd.expr), bdiv.pheno$X)


colnames(rh384.expr) <- rh384.obs$NAME[rh384.ind]
colnames(pb.10x.expr) <- plasmodium.10x.pheno$NAME[pb.10x.ind]
colnames(bd.expr) <- bdiv.pheno$NAME[bd.ind]



all.inds <- c(rh384.ind, pb.10x.ind,bd.ind)

all.samples <- data.frame(Sample = rep(NA, length(all.inds)), 
                          spp = rep(NA, length(all.inds)), stringsAsFactors = F)

all.samples$Sample <- c(as.character(rh384.obs$index[rh384.ind]),
                        as.character(plasmodium.10x.pheno$X[pb.10x.ind]),
                        as.character(bdiv.pheno$X[bd.ind]))

all.samples$spp <- c(as.character(rh384.obs$spp[rh384.ind]),
                     as.character(plasmodium.10x.pheno$spp[pb.10x.ind]),
                     as.character(bdiv.pheno$spp[bd.ind]))

all.samples$cells <- c(as.character(rh384.obs$cells[rh384.ind]),
                       as.character(plasmodium.10x.pheno$cells[pb.10x.ind]),
                       as.character(bdiv.pheno$cells[bd.ind]))


#saveRDS(all.samples, 'all_samples.rds')

#### Create Individual Seurat objects
RH.cells <- gsub('_.*', '', colnames(rh384.expr))
meta.data.RH <- data.frame(spp = rep('RH', ncol(rh384.expr)), cells = RH.cells)
rownames(meta.data.RH) <- colnames(rh384.expr)
S.O.RH348 <- CreateSeuratObject(counts = rh384.expr, meta.data = meta.data.RH, 
                                min.cells = 3, min.features = 200)
S.O.RH348 <- subset(S.O.RH348, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 )

PB.cells <- gsub('_.*', '', colnames(pb.10x.expr))
meta.data.PB <- data.frame(spp = rep('PBer', ncol(pb.10x.expr)), cells = PB.cells)
rownames(meta.data.PB) <- colnames(pb.10x.expr)
S.O.PB.10x <- CreateSeuratObject(counts = pb.10x.expr, meta.data = meta.data.PB,
                                 min.cells = 3, min.features = 200)
S.O.PB.10x <- subset(S.O.PB.10x, subset = nFeature_RNA > 200 & nFeature_RNA < 2100 )

BD.cells <- gsub('_.*', '', colnames(bd.expr))
meta.data.BD <- data.frame(spp = rep('BDiv', ncol(bd.expr)), cells = BD.cells)
rownames(meta.data.BD) <- colnames(bd.expr)
S.O.bd <- CreateSeuratObject(counts = bd.expr, meta.data = meta.data.BD,
                             min.cells = 3, min.features = 200)
S.O.bd <- subset(S.O.bd, subset = nFeature_RNA > 200 & nFeature_RNA < 1000 )



#### Use orthologs 
ortho.ind <- rownames(rh384.expr) %in% GT1.pb.bd$GT1
rh384.expr.sub <- rh384.expr[ortho.ind, ]

ortho.ind <- rownames(pb.10x.expr) %in% GT1.pb.bd$PB
pb.10x.expr.sub <- pb.10x.expr[ortho.ind, ]

ortho.ind <- rownames(bd.expr) %in% GT1.pb.bd$bdiv
bd.expr.sub <- bd.expr[ortho.ind,]

## Write all the genes using GT1 id

rownames(pb.10x.expr.sub) <- GT1.pb.bd$GT1[match(rownames(pb.10x.expr.sub), GT1.pb.bd$PB)]
rownames(bd.expr.sub) <- GT1.pb.bd$GT1[match(rownames(bd.expr.sub), GT1.pb.bd$bdiv)]


## Combine the datasets
all.expr <- inner_join(inner_join(data.frame(GeneID=rownames(rh384.expr.sub), rh384.expr.sub), 
                                  data.frame(GeneID = rownames(pb.10x.expr.sub), pb.10x.expr.sub), 
                                  by = 'GeneID'),
                       data.frame(GeneID = rownames(bd.expr.sub), bd.expr.sub), 
                       by = 'GeneID')


rownames(all.expr) <- all.expr$GeneID
all.expr <- all.expr[,-1]

all.samples$NAME <- paste(all.samples$cells, all.samples$Sample, sep = '_')
inds <- all.samples$NAME %in% colnames(all.expr)
meta.data <- data.frame(spp = all.samples$spp[inds], cells = all.samples$cells[inds])
rownames(meta.data) <- all.samples$NAME[inds]

#Create Seurat Object
all.spp <- CreateSeuratObject(all.expr, meta.data = meta.data, min.cells = 3, min.features = 200)


## Split by spp
all.spp.list <- SplitObject(all.spp, split.by = "spp")

filt.vals <- c(2500, 2100, 1000)
for (i in 1:length(all.spp.list)) {
  all.spp.list[[i]] <- prep_S.O(all.spp.list[[i]])
  all.spp.list[[i]] <- subset(all.spp.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < filt.vals[i] )
}


reference.list <- all.spp.list[c("RH", "PBer", "BDiv")]
all.samples.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

all.samples.integrated <- IntegrateData(anchorset = all.samples.anchors, dims = 1:30)


# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(all.samples.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
all.samples.integrated <- ScaleData(all.samples.integrated, verbose = FALSE)
all.samples.integrated <- RunPCA(all.samples.integrated, npcs = 30, verbose = FALSE)
all.samples.integrated <- RunUMAP(all.samples.integrated, reduction = "pca", dims = 1:30)
all.samples.integrated  <- FindNeighbors(all.samples.integrated , reduction = "pca", dims = 1:30)
all.samples.integrated  <- FindClusters(all.samples.integrated , resolution = 0.1)


p1 <- DimPlot(all.samples.integrated, reduction = "pca", group.by = "cells", 
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

ggsave(filename="../Output/scClockFigs/Integrated_scRNAseq_split.pdf", 
       plot=p1,
       width =8, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)


p2 <- DimPlot(all.samples.integrated, reduction = "umap", 
              pt.size = 1,
              shape.by='spp',
              label = TRUE, label.size = 6) + NoLegend() + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p2)



ggsave(filename="~/work/ToxoPlasmaGondiiR/Output/scFigures/Integrated_scRNAseq_merge.pdf", 
       plot=p2,
       width = 9, height = 6, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)