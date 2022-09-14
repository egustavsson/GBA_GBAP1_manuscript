
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)

# Load data ---------------------------------------------------------------

prop_ers_per_tissue <- readr::read_csv(
  here::here("results", "ers", "prop_ers_per_tissue.csv")
)

# Main --------------------------------------------------------------------

prop_ers_per_tissue_plot <- 
  prop_ers_per_tissue %>% 
  ggplot(aes(
    x = gene_type, 
    y = prop_ers
  )) + 
  geom_violin(
    aes(fill = gene_type)
  ) + 
  geom_boxplot(
    width = 0.5
  ) + 
  scale_x_discrete(name = "",
                   labels = c("Parent gene" = "Parent genes", "Protein Coding" = "Protein coding\ngenes")) + 
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"),
                     name = "Proportion of genes with\nevidence of incomplete annotation",
                     limits = c(0, 0.18)) +
  ggpubr::stat_compare_means(label.x = 1.2,
                             label.y = 0.17,
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
  plot = prop_ers_per_tissue_plot, 
  filename = "prop_ers_per_tissue.svg", 
  path = here::here("results", "ers"), 
  width = 6, 
  height = 5, 
  dpi = 600
)
