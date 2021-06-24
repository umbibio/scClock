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


getNetwork <- function(gl, prod.desc){
  
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
  
  return(bd_network)
  
}

prod.desc <- read.csv('../Input/BdivCellCycle/ProductDescription/BDvi_Prod_desc.csv')
prod.desc <- prod.desc %>% transmute(GeneID = Gene.ID, Product.Description = Product.Description)

## Read Glasso output
input.dir <- "../Input/scClock/glasso/"

in.files <- list.files(input.dir)

file.info <- data.frame(rho = rep(0, length(in.files)), size = rep(0, length(in.files)), pval = rep(0, length(in.files)))
bd_networks <- list()
for(i in 1:length(in.files)){
  gl <- readRDS(paste(input.dir, in.files[i], sep = ''))
  bd_network <- getNetwork(gl, prod.desc)
  bd_networks <- c(bd_networks, list(bd_network))
  d <- degree(bd_network, mode="in")
  pl <- power.law.fit(d+1, NULL)
  file.info$rho[i] <- as.numeric(gsub('\\.RData', '', strsplit(in.files[i], split = '_')[[1]][2]))
  file.info$size[i] <- gsize(bd_network)
  file.info$pval[i] <- pl$KS.p
}

plot(file.info$rho, file.info$pval)
which(file.info$pval < 0.6)

#plot(bd_network, vertex.cex = 1, mode = "circle")


# define a custom color palette
got_palette <- c("#1A5878", "#C44237", "#AD8941", "#E99093", "#50594B")

ggraph(bd_networks[[498]],layout = "stress")+
  geom_edge_link0(aes(edge_width = weight),edge_colour = "grey66")+
  geom_node_point(aes(fill = class,size = size),shape=21)+
  geom_node_text(aes(filter = size >= 26, label = name),family="serif")+
  scale_fill_manual(values = got_palette)+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  theme_graph()+
  theme(legend.position = "right")


xx <- as.data.frame(do.call('cbind',vertex_attr(bd_networks[[45]], index = V(bd_networks[[45]]))))
d <- degree(bd_networks[[45]], mode="in")

xx[xx$name %in% names(d[order(d,decreasing=TRUE)][1:10]),]
pl <- power.law.fit(d+1, NULL)

d = degree(bd_network, mode = "all")
hist(d)

