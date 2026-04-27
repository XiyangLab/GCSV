##Update the information of host genomes
host_A <- read.delim("gtdbtk.ar53.summary.tsv")
host_A <- host_A[,1:2]
host_B <- read.delim("gtdbtk.bac120.summary.tsv")
host_B <- host_B[,1:2]
host_taxonomy_cs <- rbind(host_A, host_B) %>%
  separate(classification, 
           into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
           sep = ";")
gtdb_ar_r220 <- read.csv('ar53_taxonomy_r220.csv', header = T)
gtdb_bac_r220 <- read.csv('bac120_taxonomy_r220.csv', header = T)
gtdb_r220 <- rbind(gtdb_ar_r220,gtdb_bac_r220)
colnames(host_taxonomy_cs) <- colnames(gtdb_r220)
host_taxonomy <- rbind(gtdb_r220, host_taxonomy_cs)

#iphop_results
iphop_results <- read.csv("Host_prediction_to_genome_m90.csv")
iphop_results_2 <- iphop_results %>% select(-Host.taxonomy)
iphop_results_2 <- merge(iphop_results_2, host_taxonomy,
                              by.x = 'Host.genome', by.y = 'genome')
iphop_results_2 <- iphop_results_2[,c(2,1,6:ncol(iphop_results_2),3:5)]

#provirus_MAGs_matched
new_provirus_MAGs_matched <- read.delim("nr_provirus_all_final_host_match_results.tsv")
new_provirus_MAGs_matched_2 <- merge(new_provirus_MAGs_matched, host_taxonomy_cs, 
                                     by.x = "host_MAG", by.y = "genome") %>%
  mutate(Main.method = c("provirus_host"), Confidence.score = NA, Additional.methods = NA) 
new_provirus_MAGs_matched_2 <- new_provirus_MAGs_matched_2[,c(2,1,3:ncol(new_provirus_MAGs_matched_2))]
colnames(new_provirus_MAGs_matched_2) <- colnames(iphop_results_2)

Host_new <- rbind(iphop_results_2, new_provirus_MAGs_matched_2)

# Define the priority of the Main.method
method_priority <- c("provirus_host" = 1, "CRISPR" = 2, "blast" = 3, "iPHoP-RF" = 4)

Host_new_filtered <- Host_new %>%
  mutate(method_rank = method_priority[Main.method]) %>%
  group_by(Virus) %>%
  filter(method_rank == min(method_rank)) %>%              # Step 1: Filter by method priority
  filter(Main.method == "provirus_host" |                  # Step 2: Remain provirus_host
           Confidence.score == max(Confidence.score)) %>%  # Step 3: Retain the results with the highest score
  #  slice(1) %>%
  ungroup() %>%
  select(-method_rank)
