library(dplyr)
library(orthologr)
library(seqinr)
library(readxl)
library(openxlsx)

## This requires installing Blast from NCBI

rec_Bdiv.vs.Bmic <- blast_rec(query_file   = "../Input/rec_blast/PiroplasmaDB-52_Bdivergens1802A_AnnotatedProteins.fasta",
                              subject_file = "../Input/rec_blast/PiroplasmaDB-52_BmicrotiRI_AnnotatedProteins.fasta",
                              delete_corrupt_cds = T, seq_type = "protein",
                              format = "fasta", blast_algorithm = "blastp",
                              eval = 0.0001)

Bdiv.Bmic.orth <- rec_Bdiv.vs.Bmic %>% dplyr::select(query_id, subject_id)

Bdiv.Bmic.orth$query_id <- gsub('-.*', '', Bdiv.Bmic.orth$query_id)
Bdiv.Bmic.orth$subject_id <- gsub('-.*', '', Bdiv.Bmic.orth$subject_id)

colnames(Bdiv.Bmic.orth) <- c('Bd', 'Bm')
write.xlsx(Bdiv.Bmic.orth, "../Input/orthologs/Bdiv_Bmic_orth.xlsx")



