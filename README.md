# scClock
## Repository for B. Divergence scRNA seq

1. `util_funcs.R` contains the functions that are needed for the analysis.
1. `ReadScRNAseqData_BD.R`  read in the 0h 37 degree WT data and generate the Seurat object. Note that the data is down-sampled to include 800 points per cluster. This is done for memory management. Markers of each cluster are identified for later use.
1. `DEAsyncTimeCourse.R` reads in the synchronized time-course data and performs a time-course DE analysis to include genes that show a significant change over time. Note that technical replicates are pooled, so there are 2 replicates per time-point. Outlier samples are excluded.
1. `fitPseudoTime.R` fits a principal curve to the first two PCA components of the sc data and diagonally projects the cells along the curve to order the cells and create a pseudo time. The pseudo time is scaled from 0 to 12 hours, binned every 20 min, and samples within each bin are considered as replicates. Genes are fitted against the pseudo time ina GAM model and genes that show no correlation with time are excluded.
1. `alignScWithSync.R` Fits smoothing splines to the gene curves and performs a cross-correlation analysis between the common genes in sc and sync data. The optimal shift is decided based on the distribution of the lag time between sc and sync common genes. The sc pseudo time is adjusted accordingly to re-map the start point.
1. `mixed_effect.R` Fits a mixed effect model to both sync and sc data to identify mean trends over time.
1. `clusterCurves.R` uses the mean-trend of Marker genes and performs a time-series clustering with dynamic time warping (dtw) to cluster the curves based on their similarity of shape, taking minor stretch and shifts into consideration. The curves are used to produce heatmaps.
1. `coldComp.R` performs a comparison between WT at 12h 4 degree. The comparison is performed in two ways. In the first samples are mixed and normalized with SCTransform. In the second approach samples are analyzed individually, clustered and Markers identified. Markers are then matched to assess similarity of identified clusters in 0h and 12h.

### Glasso
1. The Seurat objects for all three species are downsampled to 