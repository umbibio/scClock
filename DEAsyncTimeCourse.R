library(openxlsx)
library(tidyverse)
library(edgeR)
library(limma)
library(splines)
library(parallel)

## This version 
removeBatchNormalize2 <- function(count.file, expmnt.file){
  
  count <- read.xlsx(count.file)
  expmnt <- read.xlsx(expmnt.file)
  
  # cpm normalization
  CPM <- cpm(count[,2:ncol(count)]) 
  
  # filter low expressed genes
  keep <- rowSums(CPM > 2) >= 3     
  x <- count[keep, 2:ncol(count)]
  rownames(x) <- count$Gene[keep]
  treatment <- expmnt$treatment
  treatment <- factor(treatment, levels = unique(treatment))
  
  # Creating DGElist object 
  y <- DGEList(counts = x, group = treatment)  
  y <- calcNormFactors(y)
  # logarithmic scale
  logCPM <- cpm(y, log=TRUE, prior.count=3, normalized.lib.sizes=TRUE)  
  
  # Detecting batches/noisy samples
  clusters <- hclust(dist(t(logCPM)))
  clusterCut <- cutree(clusters, 4)
  CC <- data.frame(Sample = names(clusterCut), Batch = clusterCut, stringsAsFactors = F)
  expmnt <- left_join(expmnt, CC, by = "Sample")
  expmnt.cleaned <- expmnt %>% dplyr::filter(Batch != 3 & Batch !=4)
  
  x.cleaned <- x[, colnames(logCPM) %in% expmnt.cleaned$Sample]
  Samples <- data.frame(Sample = colnames(x.cleaned), 
                        time_point = expmnt.cleaned$time_point[match(colnames(x.cleaned), expmnt.cleaned$Sample)])
  Samples$ID <- gsub('CK2', 'CK', gsub('CK1', 'CK', gsub('BE2', 'BE', gsub('BE1', 'BE', Samples$Sample))))
 
  L <- list(x = x.cleaned, expmnt = expmnt.cleaned, Samples = Samples)
  return(L)
}


## Read in the raw data from synchronized time-course experiments
## Remove Batch-effects an log normalize the data
count.file <-  "../Input/BdivCellCycle/timeCourse/raw_counts_normal_growth.xlsx"
expmnt.file <- "../Input/BdivCellCycle/timeCourse/experiment_design.xlsx"

L <- removeBatchNormalize2(count.file, expmnt.file)

## Group the technical replicates and generate pool counts 
PooledCounts <- sumTechReps(L$x, ID=L$Samples$ID)

## Generate time-course data for edge-R DEA
## Note: an spline is fitted with 3 knots and F-test performed 
## to assess significant shift in trend at knot points.
Hours <- as.numeric(gsub('.*-', '', colnames(PooledCounts)))
Time <- paste0(Hours,"hrs")

y <- DGEList(counts=PooledCounts, group=Time)
y$genes <- data.frame(Symbol=rownames(y), stringsAsFactors=FALSE)

## Filter lowly expressed genes
keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
#plotMDS(y, labels=Hours)

## 3-knows natural splies
X <- poly(Hours, degree=3)
design <- model.matrix(~ X)
X <- ns(Hours, df=3)
design <- model.matrix(~X)
#design

y <- estimateDisp(y, design)
#sqrt(y$common.dispersion)
#plotBCV(y)  

fit <- glmQLFit(y, design, robust=TRUE)
#plotQLDisp(fit)

fit <- glmQLFTest(fit, coef=2:4)
tab <- as.data.frame(topTags(fit, n=nrow(y$counts)))
summary(decideTests(fit))

logCPM.obs <- cpm(y, log=TRUE, prior.count=fit$prior.count)
logCPM.fit <- cpm(fit, log=TRUE)

# par(mfrow=c(2,2))
# for(i in 1:4) {
#   GeneID <- row.names(tab)[i]
#   Symbol <- tab$Symbol[i]
#   logCPM.obs.i <- logCPM.obs[GeneID,]
#   logCPM.fit.i <- logCPM.fit[GeneID,]
#   plot(Hours, logCPM.obs.i, ylab="log-CPM", main=Symbol, pch=16)
#   lines(Hours, logCPM.fit.i, col="red", lwd=2)
# }

ind.sig <- tab$FDR[match(rownames(logCPM.obs), tab$Symbol)] < 0.01
print(sum(ind.sig)) ## number of DEGs
logCPM.obs.sig <- logCPM.obs[ind.sig,]
tc.logCPM <- logCPM.obs.sig
Samples <- L$Samples %>% transmute(Sample = ID, time_point = time_point) %>% distinct()


tc.logCPM <- tc.logCPM %>% as.data.frame() %>%
  mutate(GeneID = rownames(tc.logCPM)) %>%
  pivot_longer(-GeneID, names_to = "Sample", values_to = "expr")

tc.logCPM <- right_join(tc.logCPM, Samples, by = 'Sample')
tc.logCPM.rep <- tc.logCPM %>% group_by(GeneID, time_point) %>% summarise(rep = 1:n())
tc.logCPM$rep <- tc.logCPM.rep$rep
saveRDS(BD.markers.sig, '../Input/scClock/tc.logCPM.RData')
