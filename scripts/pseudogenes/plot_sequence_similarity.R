# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)

# Load data ---------------------------------------------------------------

parent_pseudo_homology <- read_tsv(here::here("data", "similarity.dat"), col_names = TRUE, show_col_types = FALSE)

# Main --------------------------------------------------------------------

sequence_similarity_plot <-
  parent_pseudo_homology %>% 
  ggplot(aes(x = body, fill = ..x..)) +
  geom_histogram(colour = "black", 
                 bins = 50) +
  geom_segment(aes(x = 0.96, y = 0, xend = 0.96, yend = 620), colour = "black") + 
  geom_text(aes(x=0.96, y= 680, label = "GBA1-GBAP1"), 
            colour = "black", 
            vjust = 1, 
            size = 6, fontface = "italic") +
  scale_fill_gradient(low = "#FFFFFF", high = "#ef6548", limits = c(0.35, 1)) +
  ylim(0, 700) +
  labs(x = "Sequence similarity (CDS)",
       y = "Parent-pseudogene pairs") +
  scale_x_continuous(n.breaks = 10, 
                     labels = scales::percent_format(accuracy = 1)) +
  theme_classic() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))

# Save data ---------------------------------------------------------------
ggsave(plot = sequence_similarity_plot, 
       filename = "sequence_similarity_plot.svg", 
       path = here::here("results", "pseudogenes"), 
       width = 9, 
       height = 6, 
       dpi = 600
)
