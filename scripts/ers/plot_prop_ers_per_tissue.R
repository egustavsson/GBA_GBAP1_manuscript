
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
  scale_x_discrete(name = "Gene category") + 
  scale_y_continuous(name = paste0(
    "Proportion of genes with\n", 
    "evidence of incomplete annotation"
    )) + 
  scale_fill_discrete(guide = "none") + 
  theme_bw()

# Save data ---------------------------------------------------------------

ggsave(
  plot = prop_ers_per_tissue_plot, 
  filename = "prop_ers_per_tissue.png", 
  path = here::here("results", "ers"), 
  width = 4, 
  height = 4, 
  dpi = 600
)