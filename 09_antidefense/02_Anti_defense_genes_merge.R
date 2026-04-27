library(tidyverse)

#dbAPIS_results_MERGE--------

#input
hmm_results_all <- read.delim('hmmscan_dbAPIS.out.parsed.tsv')
diamond_results_all <- read.delim('diamond_dbAPIS.out.parsed.tsv')

hmm_results_filter <- hmm_results_all %>% filter(Domain.c.evalue < 1E-10 & Domain.score > 50)
hmm_results_filter <- hmm_results_filter %>% filter (Hit.family != 'APIS028') # This family was deleted
diamond_results_filter <- diamond_results_all %>% filter(evalue < 1e-10 & pident >= 30 & qcov >= 0.8 & scov >= 0.8)

#Retain SEED
dbAPIS_SEED <- read.csv('dbAPIS_seed_info.csv')
dbAPIS_gene <- dbAPIS_SEED[,1:2]

hmm_results_SEED<-hmm_results_filter %>% filter(Defense.type != "")
diamond_results_SEED<-diamond_results_filter %>% filter(Defense.type != "")
diamond_results_SEED <- diamond_results_SEED %>% filter (famid != 'APIS063' & famid !='APIS051' & famid !='APIS067') 

hmm_results_tophit <- hmm_results_SEED %>%
  group_by(Query) %>%  
  filter(Domain.c.evalue == min(Domain.c.evalue)) %>%  
  ungroup()  

diamond_results_tophit <- diamond_results_SEED %>%
  group_by(qseqid) %>%  
  filter(evalue == min(evalue)) %>%  
  ungroup() 

#MERGE
hmm_tophit_2 <- hmm_results_tophit[,c(1,3,4,5,6,8)]
colnames(hmm_tophit_2) <- c('Query','hmm_family','hmm_defense_type','hmm_CLAN','hmm_CLAN_defense_type','hmm_evalue')
hmm_tophit_2 <- hmm_tophit_2 %>% distinct()
hmm_tophit_2 <- hmm_tophit_2 %>% 
  mutate(hmm_CLAN_defense_type = ifelse(hmm_CLAN_defense_type == "", hmm_defense_type, hmm_CLAN_defense_type))
hmm_tophit_2 <- hmm_tophit_2[,c(1,2,4,5,6)]

diamond_tophit_2 <- diamond_results_tophit[,c(1,2,3,4,5,9)]
colnames(diamond_tophit_2) <- c('Query','diamond_family','diamond_defense_type','diamond_CLAN','diamond_CLAN_defense_type','diamond_evalue')
diamond_tophit_2 <- diamond_tophit_2 %>% distinct()
diamond_tophit_2 <- diamond_tophit_2 %>% 
  mutate(diamond_CLAN_defense_type = ifelse(diamond_CLAN_defense_type == "", diamond_defense_type, diamond_CLAN_defense_type))
diamond_tophit_2 <- diamond_tophit_2[,c(1,2,4,5,6)]

dbAPIS_merge <- merge(hmm_tophit_2, diamond_tophit_2, by = 'Query') # Results supported by both methods
# Results with inconsistent annotations between the two methods
dbAPIS_merge_diff <- dbAPIS_merge %>%
  filter(hmm_family != diamond_family) 

# Prioritize HMM-based annotations
dbAPIS_merge_fixed <- dbAPIS_merge %>%
  mutate(
    source = ifelse(hmm_family == diamond_family, "dbAPIS_both", "dbAPIS_hmm")
  ) %>%
  transmute(
    Query,
    final_family = hmm_family,
    final_CLAN = hmm_CLAN,
    final_CLAN_defense_type = hmm_CLAN_defense_type,
    final_evalue = hmm_evalue,
    source
  ) %>% distinct()

# Retain sequences identified by only one method
only_hmm_tophit <- anti_join(hmm_tophit_2, dbAPIS_merge_mini, by = c("Query" = "Query"))
only_hmm_tophit <- only_hmm_tophit %>% mutate(source = 'dbAPIS_hmm')
only_diamond_tophit <- anti_join(diamond_tophit_2, dbAPIS_merge_mini, by = c("Query" = "Query"))
only_diamond_tophit <- only_diamond_tophit %>% mutate(source = 'dbAPIS_diamond')
# Merge the three tables
colnames(dbAPIS_merge_mini) <- c('Query', 'family','CLAN','defense_type','evalue','source')
colnames(only_hmm_tophit) <- c('Query', 'family','CLAN','defense_type','evalue','source')
colnames(only_diamond_tophit) <- c('Query', 'family','CLAN','defense_type','evalue','source')
dbAPIS_merge_final <- rbind(dbAPIS_merge_mini,only_diamond_tophit,only_hmm_tophit)

# Annotate gene names based on family assignments
dbAPIS_merge_parsed_SEED <- merge(dbAPIS_merge_final, dbAPIS_gene, 
                                  by.x = 'family', by.y = 'APIS.families',
                                  all.x = T)
dbAPIS_merge_parsed_SEED <- dbAPIS_merge_parsed_SEED %>%
  mutate(APIS.genes = ifelse(defense_type == 'CRISPR-Cas', family, APIS.genes))

# Extract rows with duplicated entries in the Query column for manual inspection
duplicate_row_SEED <- dbAPIS_merge_parsed_SEED %>%
  filter(duplicated(Query) | duplicated(Query, fromLast = TRUE))

write.csv(dbAPIS_merge_parsed_SEED, 'dbAPIS_merge_parsed_SEED.csv',row.names = F)
write.csv(duplicate_row_SEED, 'dbAPIS_merge_parsed_duplicate_row_SEED.csv', row.names = F)



#antidefensefinder_dbAPIS_merge------
# Input
dbAPIS <- read.csv('dbAPIS_merge_parsed_SEED_v2.csv')
dbAPIS_2 <- dbAPIS[,c(1,7,4,5,2,6)]
colnames(dbAPIS_2) <- c('Query','APIS_gene','APIS_defense_type','APIS_evalue','APIS_family','source')

antidf_df_all <- read.delim('coldseeps_vOTU_protein_20260323_defense_finder_genes.tsv')
antidf_df_rep <- antidf_df_all %>%
  filter(hit_id %in% cluster_mapping[,1])

antidf <- antidf_df_rep %>%
  filter(activity == "Antidefense")
antidf_2 <- antidf[,c(2,3,23,15,24)] 
colnames(antidf_2) <- c('Query','antidf_gene','antidf_defense_type','antidf_evalue','antidf_subtype')
antidf_2 <- antidf_2 %>% mutate(source = 'antidefensefinder')

# Standardize gene names
dbAPIS_2 <- dbAPIS_2 %>%
  mutate(
    APIS_gene = str_to_lower(APIS_gene),
    APIS_gene = str_replace_all(APIS_gene, "acriia7|tad2", "tad2_acriia7"),
    APIS_gene = str_replace_all(APIS_gene, "antidnd_p0020", "adnd_p0020"),
    APIS_gene = str_replace_all(APIS_gene, "antidnd_p0021", "adnd_p0021"),
    APIS_gene = str_replace_all(APIS_gene, "ntases", "ntase"),
    APIS_gene = str_replace_all(APIS_gene, "acriii-1", "acriii1")
  )
antidf_2 <- antidf_2 %>% 
  mutate(
    antidf_gene = str_to_lower(antidf_gene),
    antidf_gene = str_replace_all(antidf_gene, "acriia7|tad2", "tad2_acriia7")
  )

# Merge results
merge <- merge(dbAPIS_2, antidf_2, by = 'Query')
write.csv(merge, 'merge_antidf_dbAPIS_seed_0403.csv', row.names = F)
merge_diff <- merge %>% filter(APIS_gene != antidf_gene) # For manual inspection

merge_new <- read.csv('merge_antidf_dbAPIS_seed_0403_2.csv')
merge_new <- merge_new[,c(1,7,8,4,5,9,10)] %>% 
  mutate(source = 'both')

only_APIS <- anti_join(dbAPIS_2, merge_new, by = c("Query" = "Query"))
only_antidf <- anti_join(antidf_2, merge_new, by = c("Query" = "Query"))
colnames(merge_new) <- c('Query', 'gene','defense_type','APIS_evalue','APIS_family','antidf_evalue','antidf_subtype','source')
colnames(only_APIS) <- c('Query', 'gene','defense_type','APIS_evalue','APIS_family','source')
colnames(only_antidf) <- c('Query', 'gene','defense_type','antidf_evalue','antidf_subtype','source')
merge_all <- bind_rows(merge_new, only_antidf,only_APIS)

# Unify tad2 and acriia7 as tad2
merge_all <- merge_all %>%
  mutate(
    gene = ifelse(gene == "tad2_acriia7", "tad2", gene),
    defense_type = ifelse(gene == "tad2", "Anti_Thoeris", defense_type)
  )

# Standardize defense_type
merge_all_rename <- merge_all %>% mutate(
  defense_type = str_replace_all(defense_type, "CRISPR-Cas evasion by DNA repair", 
                                 "Anti_CRISPR"),
  defense_type = str_replace_all(defense_type, "restriction-modification \\(RM\\)\\;bacteriophage exclusion \\(BREX\\)", 
                                 "Anti_RM"),
  defense_type = str_replace_all(defense_type, "CBASS\\, Pycsar and CRISPR-Cas \\(type III\\)", 
                                 "Anti_CBASS"),
  defense_type = str_replace_all(defense_type, "CRISPR-Cas", 
                                 "Anti_CRISPR"),
  defense_type = str_replace_all(defense_type, "pyrimidine cyclase system for antiphage resistance \\(Pycsar\\)", 
                                 "Anti_Pycsar"),
  defense_type = str_replace_all(defense_type, "restriction-modification \\(RM\\)", 
                                 "Anti_RM"),
  defense_type = str_replace_all(defense_type, "cyclic oligonucleotide-based antiphage signaling system \\(CBASS\\)", 
                                 "Anti_CBASS"),
  defense_type = str_replace_all(defense_type, "\\bGabija\\b", 
                                 "Anti_Gabija"),
  defense_type = str_replace_all(defense_type, "\\bThoeris\\b", 
                                 "Anti_Thoeris"),
  defense_type = str_replace_all(defense_type, "\\bRecBCD\\b", 
                                 "Anti_RecBCD"),
  defense_type = str_replace_all(defense_type, "\\bRetron\\b", 
                                 "Anti_Retron"),
  defense_type = str_replace_all(defense_type, "\\bDnd\\b", 
                                 "Anti_Dnd"),
  defense_type = str_replace_all(defense_type, "toxin-antitoxin \\(TA\\)", 
                                 "Anti_TA"),
  defense_type = str_replace_all(defense_type, "NAD\\+ reconstitution pathway \\(NARP\\)", 
                                 "NADP"),
  defense_type = str_replace_all(defense_type, "\\bOther\\b", 
                                 "broad-spectrum counter-defense")
)

#viral_defense_genes-----
antidf_df_all <- read.delim('defense_finder_genes.tsv')

virus_defense <- antidf_df_all %>%
  filter(activity == "Defense")
virus_defense$gene_name <- str_replace(virus_defense$gene_name, ".*__", "")%>% 
  str_replace("\\..*", "")

count_virus_df_type <-  virus_defense %>% count(type) %>%
  arrange(desc(n))
count_virus_df_gene <-  virus_defense %>% count(type,gene_name) %>%
  arrange(desc(n))


# top 15 systems
top_systems <- count_virus_df_type %>%
  slice_head(n = 15) %>%
  pull(type)
count_virus_df_type_plot <- virus_defense %>%
  mutate(type_grouped = fct_other(type, keep = top_systems)) %>%
  count(type_grouped) %>%
  arrange(desc(n))
#plot
library(ggplot2)

color_system <- c('RM' = '#2f7f9c',        
                  'Other' = '#b3b3b3',    
                  'Cas' = '#6ba36c',       
                  'CBASS' = '#f1c40f',     
                  'SoFIC' = '#f39c12',    
                  'Toutatis' = '#E9633B',  
                  'AbiE' = '#e25d6e',       
                  'HEC-09' = '#e06e6e',     
                  'PD-Lambda-5' = '#f77b7a',
                  'Retron' = '#a8e6cf',     
                  'Gabija' = '#f9b69c',     
                  'DarTG' = '#fad6a6',      
                  'BREX' = '#63AA9C',       
                  'Gao_Iet' = '#8758a3',    
                  'Septu' = '#dcedc1',      
                  'Shedu' = '#ffab91')      


ds_plot <- ggplot(count_virus_df_type_plot, aes(x = reorder(type_grouped, -n), y = n, fill = type_grouped)) +
  geom_bar(stat = "identity") +
  labs(x = "Defense systems", y = "Count") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black")) +
  scale_fill_manual(values = color_system ) +
  labs(fill = "Defense systems")

ds_plot

#top 20 gene
top_gene <- count_virus_df_gene %>%
  slice_head(n = 20) 

dg_plot <- ggplot(top_gene, aes(x = reorder(gene_name, -n), y = n, fill = type)) +
  geom_bar(stat = "identity") +
  labs(x = "Defense genes", y = "Count") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black")) +
  scale_fill_manual(values = color_system ) +
  guides(fill = "none")
dg_plot


#viral_antidf_df------
antidf_df_all <- read.delim('coldseeps_vOTU_protein_20260323_defense_finder_genes.tsv')
virus_defense <- antidf_df_all %>%
  filter(activity == "Defense")
virus_defense_2 <- virus_defense[,c(1,2,3,23,25)]
colnames(virus_defense_2) <- c('vOTU','Query','gene','defense_type','activity')

virus_antidefense <- read.csv('merge_antidf_dbAPIS_seed_all_rename_0403.csv')
virus_antidefense <- virus_antidefense %>% 
  mutate(activity = 'antidefense') %>%
  mutate(vOTU = sub("_\\d+$", "", Query))
virus_antidefense <- virus_antidefense[,c(10,1,2,3,9)]
colnames(virus_antidefense) <- colnames(virus_defense_2)

virus_defense_antidefense <- rbind(virus_antidefense,virus_defense_2)
write.csv(virus_defense_antidefense, 'virus_defense_antidefense.csv', row.names = F)

qualified_vOTUs <- virus_defense_antidefense %>%
  distinct(vOTU, activity) %>%
  group_by(vOTU) %>%
  summarise(n = n(), .groups = "drop") %>%  
  filter(n == 2) %>%  
  pull(vOTU)  

virus_defense_antidefense_both <- virus_defense_antidefense %>% 
  filter(vOTU %in% qualified_vOTUs) %>%
  arrange(Query)
virus_defense_antidefense_both_d <- virus_defense_antidefense_both %>%
  filter(duplicated(Query) | duplicated(Query, fromLast = TRUE))
write.csv(virus_defense_antidefense_both, 'virus_defense_antidefense_both.csv', row.names = F)

#viral_antidf_host_df_match------
##input------
virus_antidefense_genes <- read.csv('merge_antidf_dbAPIS_seed_all_rename_0403.csv')
host <- read.delim("Host.tsv")
host <- host[,1:2]
colnames(host) <- c("vOTU","Host_genome")
virus_antidefense_genes <- virus_antidefense_genes %>%
  mutate(vOTU = sub("_\\d+$", "", Query))
virus_antidefense_genes_host <- merge(virus_antidefense_genes[,c(9,1:3)], host[,1:2], by = "vOTU")

host_genome_defense_genes <- read.delim("host_genome_all_protein_defense_finder_genes.tsv")
host_genome_defense_genes_2 <- host_genome_defense_genes[,c(2:4,6,24,25,27)]
colnames(host_genome_defense_genes_2) <- c("Host_genome", "Gene_id", "Gene_name",
                                           "Model", "type", "Subtype")

virus_host_merge <- merge(virus_antidefense_genes_host, host_genome_defense_genes_2, 
                          by = 'Host_genome')

##count: virus_antidefense/host_defense------
#Anti-Pycsar/Pycsar 
virus_host_Pycsar <- virus_host_merge %>% 
  filter(defense_type == 'Anti_Pycsar') %>% 
  filter(type == 'Pycsar') %>%
  count(Host_genome, type, vOTU, defense_type)

#Anti-CBASS/CBASS 
virus_host_CBASS <- virus_host_merge %>% 
  filter(defense_type == 'Anti_CBASS') %>% 
  filter(type == 'CBASS') %>%
  count(Host_genome, type, vOTU, defense_type)

#Anti-Gabija/Gabija 
virus_host_gabija <- virus_host_merge %>% 
  filter(defense_type == 'Anti_Gabija') %>% 
  filter(type == 'Gabija') %>%
  count(Host_genome, type, vOTU, defense_type)

#Anti-RM/RM 
virus_host_RM <- virus_host_merge %>% 
  filter(defense_type == 'Anti_RM') %>% 
  filter(type == 'RM') %>%
  count(Host_genome, type, vOTU, defense_type)

#Anti-CRISPR/Cas 
virus_host_CRISPR <- virus_host_merge %>% 
  filter(defense_type == 'Anti_CRISPR') %>% 
  filter(type == 'Cas') %>%
  count(Host_genome, type, vOTU, defense_type)

#Anti_Thoeris/Thoeris 
virus_host_Thoeris <- virus_host_merge %>% 
  filter(defense_type == 'Anti_Thoeris') %>% 
  filter(type == 'Thoeris') %>%
  count(Host_genome, type, vOTU, defense_type)

#Anti_TA/PsyrTA, SanaTA 
virus_host_TA <- virus_host_merge %>% 
  filter(defense_type == 'Anti_TA') %>% 
  filter(grepl("TA$", type)) %>%
  count(Host_genome, type, vOTU, defense_type)

#Anti_Retron/Retron  
virus_host_Retron <- virus_host_merge %>% 
  filter(defense_type == 'Anti_Retron') %>% 
  filter(type == 'Retron') %>%
  count(Host_genome, type, vOTU, defense_type)

#NADP/Thoeris, SEFIR, Dsr 
#Phages encoding NARP1 can overcome a diverse set of defence systems, including Thoeris, DSR1, DSR2, SIR2–HerA and SEFIR, all of which deplete NAD+ as part of their defensive mechanism.
virus_host_NADP <- virus_host_merge %>% 
  filter(defense_type == 'NADP') %>% 
  filter(type %in% c("Thoeris", "Dsr", "SEFIR") | Subtype == 'Gao_Her_SIR') %>%
  count(Host_genome, type, vOTU, defense_type)

#orf148 0
virus_host_Orf <- virus_host_merge %>% 
  filter(defense_type == 'broad-spectrum counter-defense')
virus_host_RecBCD <- virus_host_merge %>% 
  filter(defense_type == 'Anti_RecBCD')
