## Load the processed data
S.O.bd.filt <- readRDS('../Input/scClock/S.O.bd.filt.RData')
BD.markers.sig <- readRDS('../Input/scClock/BD.markers.sig.RData')
comm.genes <- readRDS('../Input/scClock/comm.genes.RData')
mu.sync.com.grid <- readRDS('../Input/scClock/mu.sync.com.grid.RData')
mu.sc.com.grid <- readRDS('../Input/scClock/mu.sc.com.grid.RData')
pc.sds.adj <- readRDS('../Input/scClock/pc.sds.adj.RData')
cell.cycle.genes.df.adj <- readRDS('../Input/scClock/cell.cycle.genes.df.adj.RData')
sync.tc.fits <- readRDS("../Input/scClock/sme_fits_sync_tc_20min.RData")
sc.tc.fits <- readRDS("../Input/scClock/sme_fits_sc_tc_20min.RData")
sync.tc.df <- readRDS('../Input/scClock/sync.tc.df.RData')
sc.tc.df.adj <- readRDS('../Input/scClock/sc.tc.df.adj.RData')
sc.hc_dtw <- readRDS('../Input/scClock/sc.hc_dtw.RData')
sync.hc_dtw <- readRDS('../Input/scClock/sync.hc_dtw.RData')
all.samples.integrated_0hr_36hrN_36hrY <- readRDS('../Input/scClock/all.samples.integrated_0hr_36hrN_36hrY.RData')
all.spp.list_0hr_36hrN_36hrY <- readRDS('../Input/scClock/all.spp.list_0hr_36hrN_36hrY.RData')
all.samples.integrated_0hr_7dN_7dY <- readRDS('../Input/scClock/all.samples.integrated_0hr_7dN_7dY.RData')
all.spp.list_0hr_7dN_7dY <- readRDS('../Input/scClock/all.spp.list_0hr_7dN_7dY.RData')
prod.desc <- readRDS('../Input/scClock/prod.desc.RData')
markers_0hr_36hrN_36hrY <- readRDS('../Input/scClock/markers_0hr_36hrN_36hrY.RData')
markers_0hr_7dN_7dY <- readRDS('../Input/scClock/markers_0hr_7dN_7dY.RData')

sc.tc.mus <- readRDS('../Input/scClock/sc.tc.mus.RData')
sync.tc.mus <- readRDS('../Input/scClock/sync.tc.mus.RData')
sc.tc.mus.scale <- readRDS('../Input/scClock/sc.tc.mus.scale.RData')
sync.tc.mus.scale <- readRDS('../Input/scClock/sync.tc.mus.scale.RData')




