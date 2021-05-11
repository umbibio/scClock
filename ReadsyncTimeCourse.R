library(openxlsx)
library(tidyverse)
library(edgeR)
library(limma)
library(stringr)
#library(TCseq)
library(fda)
library(sme)
library(dtwclust)
library(parallel)


removeBatchNormalize <- function(count.file, expmnt.file){
  
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
  
  logCPM.cleaned <- logCPM[, colnames(logCPM) %in% expmnt$Sample] 
  
  return(list(logCPM = logCPM.cleaned, expmnt = expmnt.cleaned))
  
}

## Read in the raw data from synchronized time-course experiments
## Remove Batch-effects an log normalize the data
count.file <-  "../Input/BdivCellCycle/timeCourse/raw_counts_normal_growth.xlsx"
expmnt.file <- "../Input/BdivCellCycle/timeCourse/experiment_design.xlsx"

L <- removeBatchNormalize(count.file, expmnt.file)

tc.logCPM <- L$logCPM
expmnt <- L$expmnt


tc.logCPM <- tc.logCPM %>% as.data.frame() %>%
  mutate(GeneID = rownames(tc.logCPM)) %>%
  pivot_longer(-GeneID, names_to = "Sample", values_to = "expr")

tc.logCPM <- right_join(tc.logCPM, expmnt, by = 'Sample')
tc.logCPM.rep <- tc.logCPM %>% group_by(GeneID, time_point) %>% summarise(rep = 1:n())
tc.logCPM$rep <- tc.logCPM.rep$rep

