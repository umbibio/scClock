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
inds <- c(5,6,1,2,3,7,8,10,11,9,4,12,13,14)

#seqlevels(tx_genes) <- gsub('TGME49_', '', seqlevels(tx_genes)[inds])
#seqlevels(tx_genes) <- seqlevels(tx_genes)[inds]
#isCircular(tx_genes) <- rep(F, length(isCircular(tx_genes)))

seqlevels(tx_trans) <- seqlevels(tx_trans)[inds]
isCircular(tx_trans) <- rep(F, length(isCircular(tx_trans)))

#seqinfo(tx_genes)
seqinfo(tx_trans)

counts <- Read10X_h5(filename = "../Input/scATAC/ME49_cell_ranger/raw_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../Input/scATAC/ME49_cell_ranger/singlecell.csv",
  header = TRUE,
  row.names = 1
)

peak_anno <- read_tsv("../Input/scATAC/ME49_cell_ranger/raw_peak_bc_matrix/peaks.bed")

#counts <- Read10X_h5(filename = "../Input/scATAC/ME49_cell_ranger/filtered_peak_bc_matrix.h5")
#peak_anno <- read_tsv("../Input/scATAC/ME49_cell_ranger/peak_annotation.tsv")

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  #genome = seqinfo(tx_genes),
  genome = seqinfo(tx_trans),
  fragments = '../Input/scATAC/ME49_cell_ranger/fragments.tsv.gz',
  min.cells = 3,
  min.features = 10
)

Tg_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
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

Tg_ATAC$high.tss <- ifelse(Tg_ATAC$TSS.enrichment > 2, 'High', 'Low')
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
  subset = peak_region_fragments > 10 &
    peak_region_fragments < 100 &
    pct_reads_in_peaks > 15 &
    #blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
Tg_ATAC



Tg_ATAC <- RunTFIDF(Tg_ATAC)
Tg_ATAC <- FindTopFeatures(Tg_ATAC, min.cutoff = 'q0')
Tg_ATAC <- RunSVD(Tg_ATAC)

DepthCor(Tg_ATAC)


Tg_ATAC <- RunUMAP(object = Tg_ATAC, reduction = 'lsi', dims = 2:30)
Tg_ATAC <- FindNeighbors(object = Tg_ATAC, reduction = 'lsi', dims = 2:30)
Tg_ATAC <- FindClusters(object = Tg_ATAC, verbose = FALSE, algorithm = 3)
DimPlot(object = Tg_ATAC, label = TRUE) + NoLegend()

DefaultAssay(Tg_ATAC) <- "peaks"
gene.activities <- GeneActivity(Tg_ATAC,
                                extend.downstream = 100)

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
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 1
)




## Coverage Browser
#levels(Tg_ATAC) <- c("CD4 Naive","CD4 Memory","CD8 Naive","CD8 Effector","DN T","NK CD56bright","NK CD56Dim","pre-B",'pro-B',"pDC","DC","CD14 Mono",'CD16 Mono')
Idents(Tg_ATAC) <- 'seurat_clusters'
DefaultAssay(Tg_ATAC) <- 'peaks'
CoveragePlot(
  object = Tg_ATAC,
  sep = c("-", "-"),
  #region = gsub('TGME49-', 'TGME49_', rownames(Tg_ATAC)[1:3]),
  #region = "TGME49_chrIb-516129-529789",
  region = gsub('TGME49-', 'TGME49_', VariableFeatures(Tg_ATAC)[20]),
  extend.upstream = 4000,
  extend.downstream = 2000
)


region <- StringToGRanges(regions = gsub('TGME49-', 'TGME49_', VariableFeatures(Tg_ATAC)[20]),  sep = c("-", "-"))

xx <- findOverlaps(region, tx_trans)
tx_trans[xx@to]$gene_id

DefaultAssay(Tg_ATAC) <- 'RNA'
FeaturePlot(
  object = Tg_ATAC,
  features = gsub('_', '-', tx_trans[xx@to]$gene_id),
  pt.size = 0.1,
  max.cutoff = 'q0',
  ncol = 1
)

