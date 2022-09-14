
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)

# Load data ---------------------------------------------------------------

df_lncRNA_genes_IDB_tidy <- 
  readRDS(
    here::here("results", "intron_db", "prop_genes_IDB_threshold_5.rds"))

# Main --------------------------------------------------------------------

prop_genes_threshold_5_IDB_plot <- 
  df_lncRNA_genes_IDB_tidy %>% 
  
  ggplot(aes(x = gene_category,
             y = proportion)) + 
  geom_violin(aes(fill = gene_category)) +
  geom_boxplot(width = 0.5) + 
  scale_x_discrete(name = "",
                   labels = c("parent gene" = "Parent genes", "protein coding" = "Protein coding\ngenes")) + 
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"),
                     name = "Proportion of genes with novel junctions\n misspliced in at least 5 % of individuals.") +
  ggpubr::stat_compare_means(label.x = 1.2,
                             label.y = 0.78,
                             size = 5) +
  scale_fill_manual(guide = "none",
                    values = c("#5ab4ac", "#ffffff")) + 
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))

# Save data ---------------------------------------------------------------

ggsave(
  plot = prop_genes_threshold_5_IDB_plot, 
  filename = "prop_genes_threshold_5_IDB_plot.svg", 
  path = here::here("results", "intron_db"), 
  width = 6, 
  height = 5, 
  dpi = 600
)
