# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggpubr)

# Load data ---------------------------------------------------------------

multi <-
  read.csv(
    here::here("data", "Full_multimap_parent_genes_data_combined_with_NAs.csv"))

GBA_GBAP1_multi <-
  read.csv(
    here::here("data", "Multimapping_results.txt"))


# Main --------------------------------------------------------------------

# Change name column to match between the two files
multi$sample <- sub(".[^.]+$", "", multi$sample)

# Diseaase group per sample
disease_group <- 
  multi %>% 
  dplyr::select(sample, Disease_Group) %>% 
  unique()

# Parent gene and protein coding gene multimapping

multimapping_data_to_plot <-
  multi %>% 
  na.omit() %>% 
  dplyr::filter(total_counts > 100) %>% 
  dplyr::mutate(percent_unique_counts = percent_unique_counts * 100,
                gene_category = case_when(gene_category == "protein_coding" ~ "Protein coding",
                                          gene_category %in% c("Parent_processed", "Parent_unprocessed") ~ "Parent gene",
                                          gene_name == "GBA" ~ "Parent gene"),
                facet = case_when(percent_unique_counts >= 0 & percent_unique_counts < 10 ~ "0-10%",
                                  percent_unique_counts >= 10 & percent_unique_counts < 20 ~ "10-20%",
                                  percent_unique_counts >= 20 & percent_unique_counts < 30 ~ "20-30%",
                                  percent_unique_counts >= 30 & percent_unique_counts < 40 ~ "30-40%",
                                  percent_unique_counts >= 40 & percent_unique_counts < 50 ~ "40-50%",
                                  percent_unique_counts >= 50 & percent_unique_counts < 60 ~ "50-60%",
                                  percent_unique_counts >= 60 & percent_unique_counts < 70 ~ "60-70%",
                                  percent_unique_counts >= 70 & percent_unique_counts < 80 ~ "70-80%",
                                  percent_unique_counts >= 80 & percent_unique_counts < 90 ~ "80-90%",
                                  percent_unique_counts >= 80 & percent_unique_counts <= 100 ~ "90-100%")) %>% 
  count(facet, gene_category, sample) %>% 
  group_by(gene_category, sample) %>% 
  dplyr::mutate(Percentage = (n/sum(n) *100)) 
  
# groupwise comparison 
stats <- 
  multimapping_data_to_plot %>% 
  group_by(facet) %>% 
  t_test(Percentage ~ gene_category) %>% 
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>% 
  add_xy_position(fun = "mean_sd", x = "facet", dodge = 0.8)
  
# Plot
multimapping_plot <-
ggpubr::ggbarplot(data = multimapping_data_to_plot,
                  x = "facet",
                  y = "Percentage",
                  fill = "gene_category",
                  add = c("mean_se"),
                  position = position_dodge(0.8),
                  label = TRUE, 
                  lab.nb.digits = 1) +
  stat_compare_means(label.x.npc = "center") +
  stat_pvalue_manual(stats,  
                     label = "p.adj.signif", 
                     tip.length = 0,
                     bracket.nudge.y = -2) +
  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "Unique mapping rate (%)",
       y = "Fraction of genes (%)") +
  theme(legend.title=element_blank(),
        legend.position = c(0.6, 0.8),
        legend.text = element_text(size = 12),
        legend.background = element_blank(),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(face = "bold",
                                   size = 10),
        axis.text.y = element_text(face = "bold",
                                   size = 10),
        axis.title.y = element_text(size = 14))

# GBA multimapping
GBA_multimapping_plot <-
  multi %>% 
  dplyr::filter(gene_name == "GBA",
                Disease_Group != "DLB") %>% 
  dplyr::mutate(percent_unique_counts = percent_unique_counts * 100,
                Disease_Group = case_when(Disease_Group == "Control" ~ "Control",
                                          Disease_Group != "Control" ~ "Parkinson's disease")) %>% 
  
  ggplot(aes(x = Disease_Group,
             y = percent_unique_counts)) +
  geom_violin(fill = "grey83", show.legend = F) +
  geom_boxplot(show.legend = F, 
               width=0.1) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), 
                     limits = c(0, 100)) +
  labs(title = "Short RNA-sequencing reads that uniquely map to GBA",
       x = "",
       y = "Reads that uniquely map (%)") +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position = c(0.8, 0.8),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 16))

# GBA multimapping to GBAP1
GBA_to_GBAP1_multi_plot <-
  GBA_GBAP1_multi %>% 
  left_join(., disease_group, by = c("name" = "sample")) %>% 
  dplyr::filter(Disease_Group != "DLB") %>% 
  dplyr::mutate(multi_gene_in_pseudo.vs.gene_ratio = multi_gene_in_pseudo.vs.gene_ratio * 100,
                Disease_Group = case_when(Disease_Group == "Control" ~ "Control",
                                          Disease_Group != "Control" ~ "Parkinson's disease")) %>% 
  ggplot(aes(x = Disease_Group,
             y = multi_gene_in_pseudo.vs.gene_ratio)) +
  geom_violin(fill = "grey83", show.legend = F) +
  geom_boxplot(show.legend = F, 
               width=0.1) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), 
                     limits = c(0, 100)) +
  labs(title = "GBA multimapping short RNA-sequencing reads\nthat also map to GBAP1",
       x = "",
       y = "GBA multimapping reads mapping\nto GBAP1 (%)") +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position = c(0.8, 0.8),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 16))

# Save data ---------------------------------------------------------------
ggsave(plot = multimapping_plot, 
       filename = "multimapping_plot.svg", 
       path = here::here("results", "multimapping"), 
       width = 8, 
       height = 5, 
       dpi = 600
)


ggsave(plot = GBA_multimapping_plot, 
       filename = "GBA_multimapping_plot.svg", 
       path = here::here("results", "multimapping"), 
       width = 6, 
       height = 4, 
       dpi = 600
)

ggsave(plot = GBA_to_GBAP1_multi_plot, 
       filename = "GBA_to_GBAP1_multi_plot.svg", 
       path = here::here("results", "multimapping"), 
       width = 6, 
       height = 4, 
       dpi = 600
)
