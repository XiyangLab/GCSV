virus_id <- read.delim("vOTU_v2.id", header = F)
colnames(virus_id) <- "vOTU"

#genomad-------
genomad_virus <- read.delim('genomad_coldseeps_vOTU_260319_virus_summary.tsv')
genomad_rename <- read.csv("genomad_ICTV_MSL40_rename.csv")
genomad_virus_rename <- merge(genomad_virus[,c(1,11)], genomad_rename,
                              by = "taxonomy", all.x = T) %>%
  mutate(seq_name = sub("\\|.*", "", seq_name)) %>%
  filter(Viruses != "Unclassified") %>%
  select(-Viruses, -taxonomy) %>%
  mutate(across(2:7, ~na_if(., ""))) %>%
  filter(!if_all(2:7, is.na))

#vcontact3--------
vcontact3_virus <- read.csv('vcontact3_final_assignments_combined.csv')
vcontact3_virus <- vcontact3_virus[,c(2,5,8,10,12,14,16,20,21)] %>%
  filter(Reference == "FALSE") %>%
  filter(realm..prediction. != "No Realm" & realm..prediction. != "No prediction") %>%
  mutate(across(everything(), ~ ifelse(grepl("novel_", .), NA, .))) %>%
  select(-Reference) %>%
  unite("taxonomy", 2:7, sep = ";", remove = TRUE)
mapping_ICTV_name <- read.csv("vcontact3_mapping_ICTV_name.csv")
vcontact3_virus_rename <- merge(vcontact3_virus, mapping_ICTV_name,
                                by = "taxonomy", all.x = T) %>%
  select(-taxonomy) %>%
  filter(!if_all(3:9, is.na))

#VITAP-------
VITAP_virus <- read.delim('VITAP_best_determined_lineages.tsv')
VITAP_virus <- VITAP_virus %>%
  separate(lineage, 
           into = c("species", "genus", "family", "order", "class", "phylum", "kingdom", "realm"),
           sep = ";")
VITAP_virus <- VITAP_virus[,c(1,9,8,7,6,5,4,3,10,11)] %>%
  mutate(across(2:8, ~ ifelse(grepl("-|\\[", .), NA, .)))

#Merge------
tax_ranks <- c("realm", "kingdom", "phylum", "class", "order", "family", "genus")

##Step 1: Merge the three result tables
colnames(genomad_virus_rename)<- c("vOTU",tax_ranks)
vcontact3_virus_rename_2 <- vcontact3_virus_rename[,c(1,3:9)]
colnames(vcontact3_virus_rename_2)<- c("vOTU",tax_ranks)
VITAP_virus_2 <- VITAP_virus[,1:8]
colnames(VITAP_virus_2)<- c("vOTU",tax_ranks)

combined <- bind_rows(
  genomad_virus_rename %>% mutate(software = "genomad"),
  vcontact3_virus_rename_2 %>% mutate(software = "vcontact3"),
  VITAP_virus_2 %>% mutate(software = "VITAP")
)

##Sept 2: For each taxonomic level of each vOTU_id, retain only if ≥2 tools/software agree
consensus <- combined %>%
  group_by(vOTU) %>%
  summarise(across(all_of(tax_ranks), function(x) {
    x_clean <- x[!is.na(x) & x != ""]
    if (length(x_clean) == 0) return(NA_character_)
    tab <- table(x_clean)
    best <- names(tab)[which.max(tab)]
    if (max(tab) >= 2) best else NA_character_
  }), .groups = "drop")

## Sept 3: Join with the vOTU_id table
virus_taxonomy <- virus_id %>%
  left_join(consensus, by = "vOTU")

#Plot------
#Stacked Bar Chart
plot_data <- virus_taxonomy %>%
  select(realm, kingdom, phylum, class, order, family) %>%
  pivot_longer(everything(), names_to = "level", values_to = "value") %>%
  mutate(is_annotated = ifelse(is.na(value), "Unannotated", "Annotated")) %>%
  count(level, is_annotated) %>%
  group_by(level) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(
    level = factor(level, levels = c("realm", "kingdom", "phylum", "class", "order", "family")),
    is_annotated = factor(is_annotated, levels = c( "Unannotated","Annotated"))
  )

ggplot(plot_data, aes(x = percentage, y = level, fill = is_annotated)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = c("Annotated" = "#4A90A4", "Unannotated" = "#C0C0C0")) +
  scale_x_continuous(
    labels = function(x) paste0(abs(x), "%"),
    expand = c(0, 0),
    breaks = seq(0, 100, 25)
  ) +
  labs(
    x = "Relative Proportion (%)",
    y = NULL,
    fill = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none", 
    plot.margin = margin(10, 10, 10, 10)
  ) +
  scale_y_discrete(limits = rev(levels(plot_data$level)))

#Donut Chart 1
donut_data <- virus_provirus_taxonomy %>%
  select(realm, kingdom, phylum, class, order, family) %>%
  mutate(
    is_all_na = if_all(everything(), is.na),
    status = ifelse(is_all_na, "Unannotated", "Annotated")
  ) %>%
  count(status) %>%
  mutate(
    percentage = n / sum(n) * 100,
    status = factor(status, levels = c("Unannotated", "Annotated")),
    label = sprintf("%s\n%.1f%%", status, percentage)
  )

total_vOTUs <- sum(donut_data$n)

ggplot(donut_data, aes(x = 2, y = percentage, fill = status)) +
  geom_col(width = 1, color = "white", size = 0.5) +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3.5, color = "black") +
  scale_fill_manual(values = c("Unannotated" = "#C0C0C0", "Annotated" = "#4A90A4")) +
  xlim(0.5, 2.5) +
  annotate("text", x = 0.5, y = 0, 
           label = paste0(format(annotated_total, big.mark = ","), "\nvOTUs"), 
           size = 5, fontface = "bold") +  labs(fill = NULL) +
  theme_void() +
  theme(legend.position = "none")

#Donut Chart 2
family_data <- virus_provirus_taxonomy %>%
  filter(!is.na(family)) %>%
  count(family, name = "count") %>%
  arrange(desc(count))
top9 <- family_data %>% 
  slice_head(n = 9)
others <- family_data %>% 
  slice_tail(n = nrow(family_data) - 9) %>%
  summarise(family = "Other viruses", count = sum(count))
final_data <- bind_rows(top9, others) %>%
  mutate(
    percentage = count / sum(count) * 100,
    label = ifelse(percentage > 2, sprintf("%.1f%%", percentage), ""),
    family = factor(family, levels = c(family[family != "Other viruses"], "Other viruses"))
  )

total_annotated_family <- sum(final_data$count)
colors <- c(
  setNames(c("#A8D0E6", "#B8E6B8", "#4CAF50", "#F48FB1", "#D32F2F", 
             "#FFB74D", "#FF9800", "#CE93D8", "#7B1FA2"), 
           setdiff(final_data$family, c("Other viruses"))),
  "Other viruses" = "grey"
)

ggplot(final_data, aes(x = 2, y = count, fill = family)) +
  geom_col(width = 1, color = "white", size = 0.3) +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), 
            size = 3, color = "black") +
  scale_fill_manual(values = colors, name = "Family") +
  xlim(0.2, 2.5) +
  annotate("text", x = 0.2, y = 0, 
           label = paste0(format(total_annotated_family, big.mark = ","), " vOTUs"), 
           size = 4, fontface = "bold") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  ) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE))
