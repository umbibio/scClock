library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(slingshot)
library(gam)
library(parallel)
#library(sctransform)


## Fit a pseudo-time curve and align using sync data

## down-sample the data to make it more manageable
S.O.bd.filt <- subset(x = S.O.bd, downsample = 800)

pc.bd <- getPCA(S.O.bd.filt)
sds.data <- getSlingShot(S.O.bd.filt, 'pca')
pc.sds.bd <- left_join(pc.bd, sds.data, by = "Sample")


p <- ggplot(pc.sds.bd, aes(x=PC_1,y=PC_2)) + 
  geom_point(aes(
    fill = cluster
  ), shape=21, size = 1.5)+
  geom_path(aes(x=sc1[cell.ord],y=sc2[cell.ord])) + 
  theme_bw(base_size = 14) + 
  theme(legend.position=c(1,1),legend.justification=c(1,1), 
        legend.title = element_blank(),
        legend.background = element_rect(fill=alpha('white', 0)),
        legend.direction="vertical") + 
  ylab('PC2') + xlab('PC1') + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=0.5, linetype="solid")) + 
  theme(
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  guides(color = FALSE)



plot(p)
ggsave(filename="../Output/scClockFigs/initial_pseudo_time_BD.pdf",
       plot=p,
       width = 5, height = 5,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)




## Pseudo-time analysis with SLingshot
## Identify genes that correlate with time

Y <- log2(S.O.bd.filt@assays$RNA@counts + 1)
var.genes <- names(sort(apply(Y, 1, var),decreasing = TRUE))#[1:1000] 
Y <- Y[var.genes, ]

pt <- sds.data$pt

## Map the pseudo-time to 0-12:20 hours 
t <- (12 + 1/3) * ((as.numeric(pt) - min(as.numeric(pt)))/(max(as.numeric(pt)) - min(as.numeric(pt))))
sds.data$t <- t

## time-index cells in 20 min intervals and identify cells in each partition
## They will be considered as replicates
time.breaks <- seq(1/3, 12 + 1/3, by = 1/3) 
time.idx <- rep(0, nrow(sds.data))

ind <- which(sds.data$t <= time.breaks[1])
time.idx[ind] <- 0

for(i in 2:length(time.breaks)){
  ind <- which(sds.data$t > time.breaks[(i-1)] & sds.data$t <= time.breaks[i])
  time.idx[ind] <- i - 1
}

sds.data$time.idx <- time.idx

sds.data <- sds.data %>%  
  group_by(time.idx) %>% mutate(rep = seq(1:n()))



## Run a GAM regression of expression on the pseudo-time
## Use parallel computation to speed things up. 16 cores
gam.pval <- mclapply(1:nrow(Y), function(z){
  d <- data.frame(z=as.numeric(Y[z,]), t=as.numeric(pt))
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
}, mc.cores = 16L)

gam.pval <- unlist(gam.pval)
names(gam.pval) <- rownames(Y)
## Remove the NA's and get the best fits
#gam.pval <- gam.pval[-which(is.na(gam.pval))]
gam.pval.adj <- p.adjust(gam.pval, method = 'fdr', n = length(gam.pval))
gam.pval.sig <- gam.pval[gam.pval.adj < 0.01] 
print(length(gam.pval.sig)) ## number of correlating genes

## Sort the cells on the pt
cell.ord <- sds.data$cell.ord

topgenes <- names(sort(gam.pval.sig, decreasing = FALSE))  
cell.cycle.genes.expr <- as.matrix(S.O.bd.filt@assays$RNA@data[topgenes, cell.ord])


cell.cycle.genes.df <- data.frame(GeneID = rownames(cell.cycle.genes.expr),
                                  cell.cycle.genes.expr) %>% 
  pivot_longer(-c(GeneID), names_to = 'Sample', values_to = 'log2.expr')


cell.cycle.genes.df$GeneID <- gsub('-', '_', cell.cycle.genes.df$GeneID)
cell.cycle.genes.df <- left_join(cell.cycle.genes.df, sds.data, by = 'Sample')
cell.cycle.genes.df$cluster <- S.O.bd.filt@meta.data$seurat_clusters[match(cell.cycle.genes.df$Sample, 
                                                                           rownames(S.O.bd.filt@meta.data))]
