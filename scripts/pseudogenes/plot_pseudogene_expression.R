# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggalluvial)

# Load data ---------------------------------------------------------------

load("/home/egust/Projects/pseudogenes/results/pseudogene_annotated.rda")

# Main --------------------------------------------------------------------

pseudogenes_expression_to_plot <-
  pseudogene_annotated %>%
  dplyr::filter(pseudogene_type %in% c("Unprocessed", "Processed"),
                is_expressed != "NA") %>% 
  
  ggplot(aes(axis1 = pseudogene_type,
             axis2 = is_expressed,
             fill = pseudogene_type)) +
  scale_x_discrete(limits = c("pseudogene_type", "is.expressed"),
                   labels = c("Pseudogene\ntype", "Expressed in\nGTEx"),
                   expand = c(.2, .05)) +
  geom_flow() +
  geom_stratum() +
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum))) +
  scale_fill_manual(values = c("#7570B3", "#E6AB02")) +
  theme_void() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12, 
                                   face = "bold"))

# Save data ---------------------------------------------------------------

ggsave(plot = pseudogenes_expression_to_plot, 
       filename = "pseudogenes_expression_to_plot.png", 
       path = here::here("results", "pseudogenes"), 
       width = 4, 
       height = 3, 
       dpi = 600
)
