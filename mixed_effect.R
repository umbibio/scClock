library(fda)
library(sme)
library(dtwclust)
library(parallel)

## After determining lag time with smoothing splines, fit (more expensive mix-effect) curves
## for better trend determination

## Run once or read from disk
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

sync.tc.fits <- mclapply(unique(sync.tc.df$variable), 
                         function(v) 
                           sme(sync.tc.df[sync.tc.df$variable==v,c("y","tme","ind")],
                               lambda.mu = 10, lambda.v = 10), mc.cores = num.cores)

saveRDS(object = sync.tc.fits ,file = "../Input/setClock/sme_fits_sync_tc_20min.RData")

#sync.tc.fits <- readRDS("../Input/setClock/sme_fits_sync_tc_20min.RData")


sc.tc.fits <- mclapply(unique(sc.tc.df.adj$variable),
                    function(v)
                      sme(sc.tc.df.adj[sc.tc.df.adj$variable==v,c("y","tme","ind")],
                          lambda.mu = 8, lambda.v = 8), mc.cores = num.cores)

saveRDS(object = sync.tc.fits ,file = "../Input/setClock/sme_fits_sc_tc_20min.RData")
#sc.tc.fits  <- readRDS("../Input/setClock/sme_fits_sc_tc_20min.RData")






## Plot a few curves to check the alignments

vs = unique(sc.tc.df.adj$variable)[1:16]


pdf(file = "../Output/scClockFigs/sme_fits_sc.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

par(mfrow = c(4,4))
for(v in vs){
  ind <- which(unique(sc.tc.df.adj$variable) == v)
  plot.sme(sc.tc.fits[[ind]], v)
}

dev.off()

pdf(file = "../Output/scClockFigs/sme_fits_sync.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

par(mfrow = c(4,4))
for(v in vs){
  ind <- which(unique(sync.tc.df$variable) == v)
  plot.sme(sync.tc.fits[[ind]], v)
}

dev.off()


