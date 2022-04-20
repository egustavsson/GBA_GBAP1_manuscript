# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggpubr)

# Load data ---------------------------------------------------------------

multi <-
  read.csv(
    here::here("data", "GBA_GBAP1_multimapping.csv"))

# Main --------------------------------------------------------------------

GBA_multimapping_plot <-
  multi %>% 
  
  ggplot(
    aes(
      x = factor(Disease_Group, 
                 levels = c("Control", "PD", "PDD", "DLB")),
      y = unique_mapping)) +
  geom_violin(show.legend = F, 
              colour = "Black", 
              aes(fill = Disease_Group)) +
  geom_boxplot(show.legend = F, 
               width=0.1) +
  ggpubr::stat_compare_means(label.x.npc = "center") +
  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), 
                     limits = c(0, 100)) +
  labs(x = "",
       y = "Reads that uniquely map") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold",
                                   size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 16))

# Save data ---------------------------------------------------------------

ggsave(plot = GBA_multimapping_plot, 
       filename = "GBA_multimapping_plot.png", 
       path = here::here("results", "multimapping"), 
       width = 6, 
       height = 4, 
       dpi = 600
)
