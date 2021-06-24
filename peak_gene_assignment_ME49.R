library(tidyverse)
library(bedtoolsr)

annotation.file <- '../Input/ME49/ToxoDB-52_TgondiiME49.gtf'
narrow.peaks <- '../Input/scATAC/ME49_cell_ranger/MACS2_out/Tg_MACS2_peaks.narrowPeak'

anno <- read.table(annotation.file, header = F, sep = '\t', quote = NULL)
peaks <- read.table(narrow.peaks, header = F, sep = '\t')


anno.filt <- anno %>% dplyr::filter(!grepl('KE13', V1))
peaks.filt <- peaks %>% dplyr::filter(!grepl('KE13', V1))
## First and last exons contain 5' and 3', the rest are identical to CDS
## Remove peaks overlapping CDS
anno.CDS <- anno.filt %>% dplyr::filter(V3=='CDS')
anno.CDS.peaks.overlap <- bedtoolsr::bt.intersect(a = anno.CDS, b = peaks.filt, wa = T, wb = T)

peaks.intergenic <- peaks.filt %>% dplyr::filter(!(V4 %in% anno.CDS.peaks.overlap$V13))

## Assign the remaining peaks to nearest gene
anno.trans <- anno.filt %>% dplyr::filter(V3 == 'transcript')
anno.trans.sort <- anno.trans %>% arrange(V1, V4, V5)
peaks.intergenic.sort <- peaks.intergenic %>% arrange(V1, V2, V3)
peaks.genes <- bedtoolsr::bt.closest(a = anno.trans.sort, b = peaks.intergenic.sort, D = "a", k = 3)
peaks.genes.filt <- peaks.genes %>% dplyr::filter(((V7 == '+' & V20 >= 0) | (V7 == '-' & V20 <= 0)) & (abs(V20) < 3000)) %>%
  group_by(V13) %>% slice(which.min(abs(V20)))


parse.str <- strsplit(peaks.genes.filt$V9, split = ' ')
inds <- unlist(lapply(parse.str , function(x) which(grepl("gene_id", x)) + 1))
peaks.genes.filt$gene_name <- gsub(";", "", unlist(lapply(1:length(inds), function(i) parse.str[[i]][inds[[i]]])))

peaks.genes.filt <- peaks.genes.filt %>% ungroup() %>% transmute(Chr = V1, gene_start = V4, gene_stop = V5, 
                                                                 strand = V7, gene_id = gene_name, peak_id = V13, 
                                                                 peak_start = V11, peak_stop = V12, distance = V20)
