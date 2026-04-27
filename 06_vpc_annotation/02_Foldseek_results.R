library(tidyverse)

#Merge structure-based annotation results-----
#afdb-swissprot
swissprot_db <- read.delim('idmapping_afdb_swissprot.tsv') #download from Uniprot
swissprot_db_less <- swissprot_db[,c(1,4)]

swissprot_results <- read.delim('tophit_afdb_swissprot_new_annotation_1e-10_new.tsv')
swissprot_results_less <- swissprot_results[,c(1,13,3)]
swissprot_results_annotation <- merge(swissprot_results_less, swissprot_db_less, 
                                      by.x = 'uniprot_id', by.y = 'From',
                                      all.x = F)
swissprot_results_annotation <- swissprot_results_annotation[,c(2,3,1,4)]
colnames(swissprot_results_annotation) <- c('query', 'swissprot_evalue','swissprot_id', 'swissprot_protein_name')

#bfvd
bfvd_db_1 <- read.delim('idmapping_bfvd.tsv') #download from Uniprot
bfvd_db_less <- bfvd_db[,c(1,3)]

bfvd_results <- read.delim('tophit_bfvd_new_annotation_1e-10_new.tsv')
bfvd_results_less <- bfvd_results[,c(1,13,3)]
bfvd_results_annotation <- merge(bfvd_results_less, bfvd_db_less, 
                                      by.x = 'target_id', by.y = 'From',
                                      all.x = F)
bfvd_results_annotation <- bfvd_results_annotation[,c(2,3,1,4)]
colnames(bfvd_results_annotation) <- c('query', 'bfvd_evalue','bfvd_id', 'bfvd_protein_name')

#pdb100
pdb100_results <- read.delim("tophit_pdb100_new_annotation_1e-10_annotation.tsv")
pdb100_results_less <- pdb100_results[,c(1,3,2,13)]
colnames(pdb100_results_less) <- c('query', 'pdb100_evalue','pdb100_id', 'pdb100_annotation')

#phold
phold_results <- read.delim("CS_new_gene_phold_per_cds_predictions.tsv")
phold_results <- phold_results %>%
  filter(evalue <= 1e-10)
phold_results_less <- phold_results[,c(1,9,2:4)]
colnames(phold_results_less) <- c('query','phold_evalue','phrog','phold_function','phold_product')

#merge
phold_swissprot_pdb100_bfvd <- 
  merge(phold_results_less, swissprot_results_annotation, all = T, by = 'query') %>%
  merge(pdb100_results_less, all = T, by = 'query') %>%
  merge(bfvd_results_annotation, all = T, by = 'query')


#MCP------
keywords <- c("capsid", "major head protein")

# Filter rows containing any of the specified keywords
match_rows_all_rep <- apply(phold_swissprot_pdb100_bfvd, 1, function(row) {
  any(grepl(paste(keywords, collapse = "|"), row, ignore.case = TRUE))
})
filtered_all_rep <- phold_swissprot_pdb100_bfvd[match_rows_all_rep, ]

#SAVED
#The subsequent results were checked manually
write.csv(filtered_all_rep, 'MCP_new_rep.csv', row.names = F)


#TerL-------
keywords <- c("TerL", "terminase large subunit","large terminase", "Large subunit terminase")

# Filter rows containing any of the specified keywords
match_rows_all_rep <- apply(phold_swissprot_pdb100_bfvd, 1, function(row) {
  any(grepl(paste(keywords, collapse = "|"), row, ignore.case = TRUE))
})
filtered_all_rep <- phold_swissprot_pdb100_bfvd[match_rows_all_rep, ]

#SAVED
#The subsequent results were checked manually
write.csv(filtered_all_rep, 'TerL_new_rep.csv', row.names = F)

