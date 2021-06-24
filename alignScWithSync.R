library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(slingshot)
library(gam)
#library(sctransform)

## Generate the file cell.cycle.genes.df using fitPseudoTime.R and 
## tc.logCPM using ReadsyncTimeCourse.R 

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

## As a first pass, fit smoothing splines to both sync and sc and align the splines
sc.tc.df <- cell.cycle.genes.df %>% 
  transmute(y = log2.expr, tme = t, ind = rep, variable = GeneID)

sync.tc.df <- tc.logCPM %>% 
  transmute(y = as.numeric(expr), 
            tme = as.numeric(time_point), 
            ind = rep, variable = GeneID)


saveRDS(sync.tc.df, '../Input/scClock/sync.tc.df.RData')

## Get the common genes
comm.genes <- unique(sync.tc.df$variable)[which(unique(sync.tc.df$variable) %in%
                                                  unique(sc.tc.df$variable))]

saveRDS(comm.genes, '../Input/scClock/comm.genes.RData')
## Fit smoothing splines to both and sample at regular time points (every 20 min from 0 - 12h)
mu.sync.com <- mclapply(comm.genes, function(v){
  xx <- sync.tc.df[sync.tc.df$variable == v, c("y","tme","ind")]
  mu <-  smooth.spline(x = xx$tme, y = xx$y)
  mu
}, mc.cores = num.cores)

mu.sync.com.grid <- lapply(mu.sync.com, function(mu) predict(mu, seq(0, 12, by = 1/3)))

mu.sc.com <- mclapply(comm.genes, function(v){
  xx <- sc.tc.df[sc.tc.df$variable == v, c("y","tme","ind")]
  mu <-  smooth.spline(x = xx$tme, y = xx$y)
  mu
}, mc.cores = num.cores)

mu.sc.com.grid <- lapply(mu.sc.com, function(mu) predict(mu, seq(0, 12, by = 1/3)))

## Calculate the cross-correlation between the fitted smoothing splines
cc.sc.sync.genes <- mclapply(c(1:length(comm.genes)), function(i){
  ll <- ccf(mu.sc.com.grid[[i]]$y, mu.sync.com.grid[[i]]$y, plot = F, lag.max = length(mu.sc.com.grid[[i]]$y))
  ll <- ll$lag[which.max(ll$acf)]
  ll
}, mc.cores = num.cores)


saveRDS(mu.sync.com.grid, '../Input/scClock/mu.sync.com.grid.RData')
saveRDS(mu.sc.com.grid, '../Input/scClock/mu.sc.com.grid.RData')


# Histogram with density plot
dd <- data.frame(lag = unlist(cc.sc.sync.genes))
p <- ggplot(dd, aes(x=lag)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") + theme_minimal()

plot(p)  
ggsave(filename="../Output/scClockFigs/lag_time_dist.pdf",
       plot=p,
       width = 5, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

## calculate the optimal lag time
hh <- hist(dd$lag, nclass = 50)
lag.time <- hh$breaks[which.max(hh$counts)] + 2

adjusted.time <- (sds.data$time.idx * 1/3) -  sort(unique(sds.data$time.idx) * 1/3)[lag.time]
neg.ind <- ifelse(adjusted.time < 0, T, F)
adjusted.time[neg.ind] <- adjusted.time[neg.ind] + (12 + 1/3)
sds.data$adj.time <-  adjusted.time

## create fine-resolution 20 min cluster of cells
clusters <- paste('C', 1:length(unique(sds.data$adj.time)), sep = '')
sds.data$cluster <- clusters[as.integer((sds.data$adj.time) * 3 + 1)]


## Generate shifted curves
time.breaks <- seq(1/3, 12 + 1/3, by = 1/3) 
time.idx <- rep(0, nrow(sds.data))

ind <- which(sds.data$adj.time <= time.breaks[1])
time.idx[ind] <- 0

for(i in 2:length(time.breaks)){
  ind <- which(sds.data$adj.time > time.breaks[(i-1)] & sds.data$adj.time <= time.breaks[i])
  time.idx[ind] <- i - 1
}

sds.data$adj.time.idx <- time.idx



sds.data <- sds.data %>%  ungroup() %>%
  group_by(adj.time.idx) %>% mutate(rep = seq(1:n()))

sds.data <- as.data.frame(sds.data)
rownames(sds.data) <- sds.data$Sample

## Add the new clusters as meta-data
S.O.bd.filt <- AddMetaData(S.O.bd.filt, sds.data)

saveRDS(S.O.bd.filt, '../Input/scClock/S.O.bd.filt.RData')

pc <- S.O.bd.filt[['pca']]@cell.embeddings
pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc))
pc$cluster <- S.O.bd.filt$cluster



pc.sds.adj <- left_join(pc, sds.data, by = "Sample")

lvs <- paste('C', unique(sort(as.numeric(gsub('C', '', pc.sds.adj$cluster.y)))), sep = '')
pc.sds.adj$cluster.y <- factor(pc.sds.adj$cluster.y, levels = lvs)

saveRDS(pc.sds.adj, '../Input/scClock/pc.sds.adj.RData')

## Plot the new clusters and mark the time start point
p <- ggplot(pc.sds.adj, aes(x=PC_1,y=PC_2)) +
  geom_point(aes(
    fill = cluster.y
  ), shape=21, size = 1.5)+
  geom_path(aes(x=sc1[cell.ord],y=sc2[cell.ord])) +
  geom_point(aes(x=sc1[order(adj.time)][1],y=sc2[order(adj.time)][1]), col = 'red', shape=8, size = 3) +
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


## Too many clusters? Maybe reduce to 8? For now just remove the legend
p <- ggplot(pc.sds.adj, aes(x=PC_1,y=PC_2)) +
  geom_point(aes(
    fill = cluster.y
  ), shape=21, size = 1.5)+
  geom_path(aes(x=sc1[cell.ord],y=sc2[cell.ord])) +
  geom_point(aes(x=sc1[order(adj.time)][1],y=sc2[order(adj.time)][1]), col = 'black', shape=21, size = 5, stroke = 1.2)+
  geom_point(aes(x=sc1[order(adj.time)][1],y=sc2[order(adj.time)][1]), col = 'blue', shape=8, size = 4, stroke = 1.1)+
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
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
ggsave(filename="../Output/scClockFigs/aligned_pseudo_time_BD.pdf",
      plot=p,
      width = 5, height = 5,
      units = "in", # other options are "in", "cm", "mm"
      dpi = 300
)



cell.cycle.genes.df.adj <- left_join(cell.cycle.genes.df, sds.data[,c('Sample', 'adj.time', 
                                                                      'adj.time.idx', 'rep', 'cluster')], by = 'Sample')
saveRDS(cell.cycle.genes.df.adj, '../Input/scClock/cell.cycle.genes.df.adj.RData')

sc.tc.df.adj <- cell.cycle.genes.df.adj %>% 
  transmute(y = log2.expr, tme = adj.time, ind = rep.y, variable = GeneID)

saveRDS(sc.tc.df.adj, '../Input/scClock/sc.tc.df.adj.RData')


