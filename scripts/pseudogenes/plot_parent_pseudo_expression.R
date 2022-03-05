# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)

# Load data ---------------------------------------------------------------

load("/home/egust/Projects/pseudogenes/results/pseudogene_annotated.rda")
load("/home/egust/Projects/pseudogenes/results/parentgene_annotated.rda")

# Main --------------------------------------------------------------------

parent_pseudo_to_plot <-
  
  rbind(
    parentgene_annotated %>%
      dplyr::filter(associated_pseudo %in% c("Unprocessed", "Processed")) %>% 
      dplyr::select(no_tissues_expressed,
                    type = associated_pseudo) %>%
      dplyr::mutate(gene = "Parent genes"),
    
    pseudogene_annotated %>%
      dplyr::filter(pseudogene_type %in% c("Unprocessed", "Processed")) %>% 
      dplyr::select(no_tissues_expressed,
                    type = pseudogene_type) %>% 
      dplyr::mutate(gene = "Pseudogenes")
  ) %>% 
  
  
  ggplot(aes(x = no_tissues_expressed, 
             fill = type))+
  geom_histogram(aes(y=..density..),
                 position = "dodge", 
                 colour = "black", 
                 binwidth = 1) + 
  theme_classic() +
  labs(x = "Number of tissues", y = "Fraction of pseudo genes") +
  scale_fill_manual(values = c("#7570B3", "#E6AB02")) +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(face = "bold",
                                   size = 8),
        axis.text.y = element_text(face = "bold",
                                   size = 12),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(face = "bold",
                                    size = 12),
        strip.background =element_rect(fill="gray63"),
        legend.title = element_blank()) +
  facet_wrap(~gene)

# Save data ---------------------------------------------------------------

ggsave(plot = pseudogenes_expression_to_plot, 
       filename = "pseudogenes_expression_to_plot.png", 
       path = here::here("results", "pseudogenes"), 
       width = 4, 
       height = 3, 
       dpi = 600
)