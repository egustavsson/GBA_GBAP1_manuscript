
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
  scale_x_discrete(name = "Gene category") + 
  scale_y_continuous(name = paste0(
    "Proportion of genes with novel junctions misspliced in\nat least ",
    "5 % of individuals.")) +
  ggpubr::stat_compare_means(label.x.npc = "center",
                             label.y.npc = "top",
                             size = 6) +
  scale_fill_discrete(guide = "none") + 
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))

# Save data ---------------------------------------------------------------

ggsave(
  plot = prop_genes_threshold_5_IDB_plot, 
  filename = "prop_genes_threshold_5_IDB_plot.png", 
  path = here::here("results", "intron_db"), 
  width = 6, 
  height = 5, 
  dpi = 600
)
