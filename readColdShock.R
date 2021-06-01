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


source("./util_funcs.R")

set.seed(100)



prod.desc <- read.csv('../Input/BdivCellCycle/ProductDescription/BDvi_Prod_desc.csv')
prod.desc <- prod.desc %>% transmute(GeneID = Gene.ID, Product.Description = Product.Description)

#saveRDS(prod.desc, '../Input/scClock/prod.desc.RData')

### b. divergense
input.dir.bdiv <- "../Input/scRNAseqBdivCS/"
count.files <- list.files(input.dir.bdiv)
num.total.files <- length(count.files)

file.info <- data.frame(time = rep(NA, num.total.files), 
                        reactivate = rep(NA, num.total.files), 
                        filename = rep(NA, num.total.files))

S.O.list <- list()

for(i in 1:num.total.files){
  file.info$filename[i] <- count.files[i]
  file.info$time[i] <- gsub('OUT', '', gsub('bd', '', strsplit(file.info$filename[i], split = '_')[[1]][1]))
  file.info$reactivate[i] = ifelse(grepl('OUT', file.info$filename[i]), 'Y', 'N')
  cat(paste('processing file', file.info$filename[i]))
  cat('\n')
  L <- processCount(input.dir.bdiv, file.info$filename[i], file.info$time[i], file.info$reactivate[i])
  S.O.list <- c(S.O.list, list(L))
}

file.info$spp <- paste('BDiv', file.info$time, file.info$reactivate, sep = '')

print(file.info)


