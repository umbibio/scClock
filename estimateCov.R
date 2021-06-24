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
library(glassoFast)
library(igraph)
library(ggraph)
library(graphlayouts)




prod.desc <- read.csv('../Input/BdivCellCycle/ProductDescription/BDvi_Prod_desc.csv')
prod.desc <- prod.desc %>% transmute(GeneID = Gene.ID, Product.Description = Product.Description)

## Read the WT bdiv and create an initial Seurat object
file.counts <- read.csv("../Input/scRNAseqBdivCS/bd0hr_count.csv")
genes <- file.counts$X
expr <- file.counts[,-1]
rownames(expr) <- genes

S.O <- CreateSeuratObject(counts = expr, min.cells = 10, min.features = 100)


#VlnPlot(S.O, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#FeatureScatter(S.O, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

cutoffs <- quantile(S.O$nCount_RNA, probs = c(0.01, 0.9))
print(cutoffs)
S.O <- subset(S.O, subset = nFeature_RNA > cutoffs[1] & nFeature_RNA < cutoffs[2] )

## Downsample the cells
S.O.sub <- subset(x = S.O, downsample = 3000)

pheno <- data.frame(Sample = names(S.O.sub$orig.ident))
spp <- 'BDiv0hr'
pheno$spp <- spp
pheno$NAME <- paste(pheno$spp, pheno$Sample, sep = '_')
rownames(pheno) <- names(S.O.sub$orig.ident)
S.O.sub <- AddMetaData(S.O.sub, metadata = pheno) 

## Find variable Features and fit a KNN graph on PCA coordinates
S.O.sub <- NormalizeData(S.O.sub)
S.O.sub <- FindVariableFeatures(S.O.sub, selection.method = "vst", nfeatures = 800)
S.O.sub <- ScaleData(S.O.sub, verbose = FALSE)
S.O.sub <- RunPCA(S.O.sub, npcs = 30, verbose = FALSE)
S.O.sub <- RunUMAP(S.O.sub, reduction = "pca", dims = 1:30)
S.O.sub <- FindNeighbors(S.O.sub, reduction = "pca", dims = 1:10)
S.O.sub <- FindClusters(S.O.sub, resolution = 0.2)

DimPlot(S.O.sub, reduction = 'pca')


### Do a network smoothing

## Smoothing over KNN
AdjacencyMat <- as.matrix(S.O.sub@graphs$RNA_nn)

con <- colSums(AdjacencyMat)
con <- unlist(lapply(con, function(x) ifelse(x == 0, 0, 1/x)))

## Scaling adjacancy
for(i in 1:ncol(AdjacencyMat)){
  AdjacencyMat[,i] <- con[i] * AdjacencyMat[,i]
}

## Get the expression data
Idents(S.O.sub) <- 'Sample'
exprs <- AverageExpression(S.O.sub, features = VariableFeatures(object = S.O.sub))

## Smoothing for a few iterations
max.smoothing <- 6
alpha <- 0.2
scDat <- exprs$RNA
Ft <- scDat
for(i in 1:max.smoothing){
  Ft <- alpha * Ft %*% AdjacencyMat + (1 - alpha) * scDat
}

S.O.smooth <- S.O.sub
smooth_assay <- CreateAssayObject(counts = Ft)
S.O.smooth[["smooth"]] <- smooth_assay


Assays(S.O.smooth)
Idents(S.O.smooth) <- 'seurat_clusters'

DefaultAssay(S.O.smooth) <- "RNA"
FeaturePlot(S.O.smooth, "Bdiv-025600", reduction = 'pca')

DefaultAssay(S.O.smooth) <- "smooth"
FeaturePlot(S.O.smooth, "Bdiv-025600", reduction = 'pca')


sds <- apply(Ft, 1, sd)
rm.ind <- which(sds == 0)

if(length(rm.ind) > 0){
  Ft <- Ft[-rm.ind,]
}

genes <- rownames(Ft)
cells <- colnames(Ft)

Ft.scale <- scale(t(Ft))
S <- cov(Ft.scale)

lambda.max <- max(abs(S[col(S) != row(S)]))
lambda.min <- lambda.max * 0.1

RHO <- matrix(lambda.min, nrow = nrow(S), ncol = ncol(S))
diag(RHO ) <- 0

gl <- glassoFast(S, rho = RHO)

par.cor <- -cov2cor(gl$wi)
cutoff.cor <- par.cor %>% as.data.frame() %>% summarise(across(everything(), ~ quantile(abs(.x), prob = 0.25))) %>% 
  quantile(prob = 0.75)



recovered.network <- par.cor
recovered.network[which(abs(par.cor) <=  as.numeric(cutoff.cor))] <- 0
#recovered.network[which(abs(par.cor) >  as.numeric(cutoff.cor))] <- 1
diag(recovered.network) <- 0
recovered.network[upper.tri(recovered.network)] <- 0
colnames(recovered.network) <- genes

recovered.network.adj <- recovered.network
recovered.network.adj[which(abs(par.cor) >  as.numeric(cutoff.cor))] <- 1

node.con <- colSums(recovered.network.adj)
node.con <- data.frame(id = names(node.con), size = as.numeric(node.con))

edge.list <- recovered.network %>% as.data.frame() %>% mutate(from = genes) %>% 
  pivot_longer(-from, names_to = 'to', values_to = 'weight') %>% dplyr::filter(weight != 0) %>%
  transmute(from = from, to = to, weight = abs(weight))



nodes <- data.frame(GeneID = gsub("-", "_", genes), id = genes) %>% left_join(prod.desc, by = 'GeneID')

nodes.list <- nodes %>% dplyr::filter(id %in% unique(c(edge.list$from, edge.list$to))) %>%
  transmute(id = id, type = Product.Description) %>% left_join(node.con, by = 'id')

nodes.list$class <- 'Protein'
nodes.list$class[grepl('hypothetical', nodes.list$type)] <- 'hypo'
nodes.list$class[grepl('conserved', nodes.list$type)] <- 'cons'
nodes.list$class[grepl('histone', nodes.list$type)] <- 'hist'
nodes.list$class[grepl('domain', nodes.list$type)] <- 'domain'

bd_network <- graph_from_data_frame(d = edge.list, vertices = nodes.list, directed = F)


#plot(bd_network, vertex.cex = 1, mode = "circle")


# define a custom color palette
got_palette <- c("#1A5878", "#C44237", "#AD8941", "#E99093", "#50594B")

ggraph(bd_network,layout = "stress")+
  geom_edge_link0(aes(edge_width = weight),edge_colour = "grey66")+
  geom_node_point(aes(fill = class,size = size),shape=21)+
  geom_node_text(aes(filter = size >= 26, label = name),family="serif")+
  scale_fill_manual(values = got_palette)+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  theme_graph()+
  theme(legend.position = "right")


d <- degree(bd_network, mode="in")

pl <- power.law.fit(d+1, NULL)

d = degree(bd_network, mode = "all")
hist(d)

