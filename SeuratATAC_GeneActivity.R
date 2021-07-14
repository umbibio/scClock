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
library(Signac)
library(Seurat)
library(patchwork)
library(hdf5r)
library(GenomeInfoDbData)
library(GenomicRanges)
library(GenomicAlignments)
library(Biostrings)
library(rtracklayer)
library(GenomicFeatures)

ME49.fasta <- readDNAStringSet("../Input/ME49/ToxoDB-52_TgondiiME49_Genome.fasta")
chrs <- names(ME49.fasta)[grep("TGME49_chr", names(ME49.fasta))]

chr.len <- data.frame(chr = gsub(" ", "", unlist(lapply(strsplit(chrs, split = '\\|'), `[[`, 1))),
                      len = as.numeric(gsub('length=', '', unlist(lapply(strsplit(chrs, split = '\\|'), `[[`, 4)))))

txdb <- makeTxDbFromGFF(file="../Input/ME49/ToxoDB-52_TgondiiME49_filter.gtf",
                        dataSource="Toxodb",
                        organism="Toxoplasma")

trans_biotypes <- select(txdb, keys=keys(txdb, "TXID"), 
                         columns = "TXTYPE", keytype =  "TXID")


genome(txdb) <- 'ME49'
#tx_genes <- genes(txdb)

#tx_trans <- unlist(transcriptsBy(txdb, by = c("gene", "exon", "cds")))

tx_trans <- exonsBy(txdb, by = "tx", use.names = TRUE)
tx_names <- names(tx_trans)
num.exons <- lapply(tx_trans, function(x) length(x))
tx_names <- rep(tx_names, unlist(num.exons))
tx_trans <- unlist(tx_trans)
tx_trans$tx_id <- tx_names
tx_trans$gene_id <- gsub('-t.*', '', tx_trans$tx_id)
tx_trans$gene_name <- tx_trans$gene_id
tx_trans$type <- 'exon'
tx_trans$gene_biotype <- 'protein_coding'
tx_trans$exon_name <- tx_trans$exon_rank

tmp <- chr.len$len[match(names(seqlengths(txdb)), chr.len$chr)]
names(tmp) <- names(seqlengths(txdb))
#seqlengths(tx_genes) <- tmp
seqlengths(tx_trans) <- tmp
#seqlevels(tx_genes)
seqlevels(tx_trans)
#inds <- c(c(5,6,1,2,3,7,8,10,11,9,4,12,13,14) + 100, seq(1, 100, by = 1))
inds <- c(5,6,1,2,3,7,8,10,11,9,4,12,13,14) 

#seqlevels(tx_genes) <- gsub('TGME49_', '', seqlevels(tx_genes)[inds])
#seqlevels(tx_genes) <- seqlevels(tx_genes)[inds]
#isCircular(tx_genes) <- rep(F, length(isCircular(tx_genes)))

seqlevels(tx_trans) <- seqlevels(tx_trans)[inds]
isCircular(tx_trans) <- rep(F, length(isCircular(tx_trans)))

#seqinfo(tx_genes)
seqinfo(tx_trans)

counts <- Read10X_h5(filename = "../Input/scATAC/NovaSeq/Tg_ME49/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../Input/scATAC/NovaSeq/Tg_ME49/singlecell.csv",
  header = TRUE,
  row.names = 1
)

metadata.filt <- metadata
metadata.filt$Sample <- rownames(metadata.filt)
metadata.filt <- metadata.filt[metadata.filt$Sample %in% colnames(counts), ]
peak_anno <- read_tsv("../Input/scATAC/NovaSeq/Tg_ME49/filtered_peak_bc_matrix/peaks.bed")

#counts <- Read10X_h5(filename = "../Input/scATAC/ME49_cell_ranger/filtered_peak_bc_matrix.h5")
#peak_anno <- read_tsv("../Input/scATAC/ME49_cell_ranger/peak_annotation.tsv")

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = seqinfo(tx_trans),
  fragments = '../Input/scATAC/NovaSeq/Tg_ME49/fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

Tg_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata.filt
)


Tg_ATAC[['peaks']]
granges(Tg_ATAC)
#annotations <- tx_genes
annotations <- tx_trans

#annotations$gene_biotype <- 'protein_coding'
#annotations$gene_name <- annotations$gene_id
#annotations$gene_name <- gsub('-t.*', '', annotations$tx_name)
#annotations$gene_id <- annotations$gene_name
#annotations$type <- 'exon'

Annotation(Tg_ATAC) <- annotations
Tg_ATAC <- NucleosomeSignal(object = Tg_ATAC)
Tg_ATAC <- TSSEnrichment(object = Tg_ATAC, fast = FALSE)

peaks <- CallPeaks(
  object = Tg_ATAC,
  macs2.path = "/Users/kouroshz/miniconda3/envs/macs2/bin/macs2",
  extsize = 100,
  additional.args = "--nomodel -B --SPMR"
)


# fragments <- CreateFragmentObject(
#   path = "../Input/scATAC/ME49_cell_ranger/fragments.tsv.gz",
#   cells = colnames(Tg_ATAC),
#   validate.fragments = FALSE
# )
# #> Computing hash
# Fragments(Tg_ATAC) <- fragments





# add blacklist ratio and fraction of reads in peaks
Tg_ATAC$pct_reads_in_peaks <- Tg_ATAC$peak_region_fragments / Tg_ATAC$passed_filters * 100
Tg_ATAC$blacklist_ratio <- Tg_ATAC$blacklist_region_fragments / Tg_ATAC$peak_region_fragments

Tg_ATAC$high.tss <- ifelse(Tg_ATAC$TSS.enrichment > 3, 'High', 'Low')
TSSPlot(Tg_ATAC, group.by = 'high.tss') + NoLegend()

Tg_ATAC$nucleosome_group <- ifelse(Tg_ATAC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = Tg_ATAC, group.by = 'nucleosome_group')

VlnPlot(
  object = Tg_ATAC,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)


Tg_ATAC <- subset(
  x = Tg_ATAC,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 4000 &
    pct_reads_in_peaks > 40 &
    #blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
Tg_ATAC



Tg_ATAC <- RunTFIDF(Tg_ATAC)
Tg_ATAC <- FindTopFeatures(Tg_ATAC, min.cutoff = 'q0')
Tg_ATAC <- RunSVD(Tg_ATAC)

DepthCor(Tg_ATAC)

## Must remove highly correlating components
Tg_ATAC <- RunUMAP(object = Tg_ATAC, reduction = 'lsi', dims = seq(1:30)[-c(1,5)])
Tg_ATAC <- FindNeighbors(object = Tg_ATAC, reduction = 'lsi', dims = seq(1:30)[-c(1,5)])
Tg_ATAC <- FindClusters(object = Tg_ATAC, verbose = FALSE, algorithm = 3)


DimPlot(object = Tg_ATAC, label = TRUE, reduction = 'tsne') + NoLegend()

DefaultAssay(Tg_ATAC) <- "peaks"
gene.activities <- GeneActivity(Tg_ATAC, extend.upstream = 600,
                                extend.downstream = 200)

## gene.activity matrix created with other approaches can be passed here
Tg_ATAC[['RNA']] <- CreateAssayObject(counts = gene.activities)
Tg_ATAC <- NormalizeData(
  object = Tg_ATAC,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(Tg_ATAC$nCount_RNA)
)

DefaultAssay(Tg_ATAC) <- 'RNA'

FeaturePlot(
  object = Tg_ATAC,
  features = c('TGME49-208020'),
  pt.size = 0.5,
  max.cutoff = 'q95',
  ncol = 1
)



## For PCA
DefaultAssay(Tg_ATAC) <- 'RNA'
Tg_ATAC <- FindVariableFeatures(Tg_ATAC, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Tg_ATAC)
Tg_ATAC <- ScaleData(Tg_ATAC, features = all.genes)
Tg_ATAC <- RunPCA(Tg_ATAC, features = VariableFeatures(object = Tg_ATAC))
Tg_ATAC <- FindNeighbors(Tg_ATAC, dims = 1:10, reduction = 'pca')
Tg_ATAC <- FindClusters(Tg_ATAC, resolution = 0.2)
Tg_ATAC<- RunTSNE(object = Tg_ATAC,features = VariableFeatures(object = Tg_ATAC) )
DimPlot(object = Tg_ATAC, reduction = "tsne", label = TRUE) + NoLegend()

pc.ATAC <- getPCA(Tg_ATAC)
p1 <- ggplot(pc.ATAC, aes(x=PC_1,y=PC_3)) + 
  geom_point(aes(
    fill = cluster
  ), shape=21, size = 1.5)+ 
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
plot(p1)
##


## Coverage Browser
#levels(Tg_ATAC) <- c("CD4 Naive","CD4 Memory","CD8 Naive","CD8 Effector","DN T","NK CD56bright","NK CD56Dim","pre-B",'pro-B',"pDC","DC","CD14 Mono",'CD16 Mono')
Idents(Tg_ATAC) <- 'seurat_clusters'
DefaultAssay(Tg_ATAC) <- 'peaks'
CoveragePlot(
  object = Tg_ATAC,
  sep = c("-", "-"),
  #region = gsub('TGME49-', 'TGME49_', rownames(Tg_ATAC)[1:3]),
  region = "TGME49_chrIb-516129-529789",
  #region = gsub('TGME49-', 'TGME49_', VariableFeatures(Tg_ATAC)[20]),
  extend.upstream = 4000,
  extend.downstream = 2000
)

##Find Markers
DefaultAssay(Tg_ATAC) <- 'peaks'
da_peaks <- FindMarkers(
  object = Tg_ATAC,
  ident.1 = "0",
  ident.2 = "1",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

head(da_peaks)

plot1 <- VlnPlot(
  object = Tg_ATAC,
  features = rownames(da_peaks)[1],
  pt.size = 0.4,
  idents = c("0","1")
)

plot2 <- FeaturePlot(
  object = Tg_ATAC,
  features = rownames(da_peaks)[1],
  reduction = 'pca',
  pt.size = 0.4
)

plot1 | plot2

head(da_peaks)


region <- StringToGRanges(regions = gsub('TGME49-', 'TGME49_', rownames(da_peaks)[1]),  sep = c("-", "-"))

xx <- findOverlaps(region, tx_trans)
tx_trans[xx@to]$gene_id

my.gene <- tx_trans[xx@to]$gene_id[2]

DefaultAssay(Tg_ATAC) <- 'RNA'

p1 <- FeaturePlot(
  object = Tg_ATAC,
  features = gsub('_', '-', my.gene),
  pt.size = 0.4,
  max.cutoff = 'q0',
  ncol = 1,
  reduction = 'tsne'
)


ggsave(filename="../Output/scClockFigs/RON2_ATAC_activity.pdf",
       plot=p1,
       width = 5, height = 5,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


Idents(Tg_ATAC) <- 'seurat_clusters'
DefaultAssay(Tg_ATAC) <- 'peaks'
p2 <- CoveragePlot(
  object = Tg_ATAC,
  sep = c("-", "-"),
  #region = gsub('TGME49-', 'TGME49_', rownames(Tg_ATAC)[1:3]),
  region = region,
  #region = gsub('TGME49-', 'TGME49_', VariableFeatures(Tg_ATAC)[20]),
  extend.upstream = 4000,
  extend.downstream = 2000
)

ggsave(filename="../Output/scClockFigs/RON2_track_pileup.pdf",
       plot=p2,
       width = 5, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


p3 <- DimPlot(object = Tg_ATAC, label = TRUE, reduction = 'tsne') + NoLegend()

ggsave(filename="../Output/scClockFigs/tsne_clusters.pdf",
       plot=p3,
       width = 5, height = 5,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

