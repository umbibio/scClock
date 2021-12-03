library(openxlsx)
## Get the gene clusters
tc_clusts <- read.xlsx('../Output/scClockOut/sc_time_course_gene_clusters.xlsx')
tc_clusts <- tc_clusts %>% dplyr::select(GeneID, cluster, phase) %>% distinct()
tc_clusts$GeneID <- gsub('___.*', '', tc_clusts$GeneID)
clust.nums <- tc_clusts %>% group_by(phase) %>% summarise(total.genes = n())
orthologs.Bd.Tg <- read.xlsx('../Input/orthologs/GT1_BDiv.xlsx')

tc_clusts.orth <- inner_join(tc_clusts, orthologs.Bd.Tg, by = c('GeneID' = 'BDiv'))

tc_clusts.orth.gene.lists <- tc_clusts.orth %>% group_by(phase) %>%
  summarise(GT1.orth = list(GT1))

write.xlsx(tc_clusts.orth.gene.lists, '../Input/scClock/GO/tc_clusts_orth_gene_lists.xlsx')

## The genes in the above file were inputted into ToxoDB to perform geneset enrichment analysis
## on Cellular component, molecular function, and Biological Process

in.dir <- '../Input/scClock/GO/ToxoDB_GO/'
all.files <- list.files(in.dir)

all.clust.items <- list()
for(f in all.files){
  nn <- gsub('\\.tsv', '', f)
  GF <- strsplit(nn, split = '_')[[1]][1]
  phase <- strsplit(nn, split = '_')[[1]][2]
  tmp <- read_tsv(paste(in.dir, f, sep = ''))
  tmp$GF <- GF
  tmp$phase <- phase
  all.clust.items <- c(all.clust.items, list(tmp))
}


all.clust.items <- do.call(rbind, all.clust.items)

saveRDS(all.clust.items, '../Input/BdivCellCycle/RDS/all_clust_items.RData')

filtered.Go <- all.clust.items %>% arrange(phase, Benjamini) %>% distinct() %>%
  group_by(phase) %>% mutate(rank = row_number()) %>%
  dplyr::filter(Benjamini < 0.1 & rank < 30) %>% 
  arrange(phase, Benjamini) 


filtered.Go$phase <- factor(filtered.Go$phase)
filtered.Go$ID <- factor(filtered.Go$ID, level=unique(filtered.Go$ID))
filtered.Go$Name <- factor(filtered.Go$Name, level=unique(filtered.Go$Name))
## Category of contrasts
p <- ggplot(filtered.Go, aes(x = phase, y = Name)) + 
  geom_point(aes(colour = "red", size = -log(Benjamini))) +
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face="bold")) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face="bold")) +
  theme(legend.position="none") +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=0.5, linetype="solid")) + guides(color = FALSE)

plot(p)


ggsave(filename="../Output/scClockFigs/GO_gene_clusters.pdf", 
       plot=p,
       width = 8, height = 12, 
       units = "in", # other options are "in", "cm", "mm" 
       dpi = 300
)
