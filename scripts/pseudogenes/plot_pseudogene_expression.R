# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)

# Load data ---------------------------------------------------------------

load(here::here("results", "pseudogenes", "pseudogene_annotated.rda"))

# Main --------------------------------------------------------------------

pseudogene_expression_to_plot <- 
  pseudogene_annotated %>%
      dplyr::filter(pseudogene_type %in% c("Unprocessed", "Processed")) %>%
      dplyr::select(no_tissues_expressed,
                    type = pseudogene_type) %>%
      dplyr::mutate(gene = "Pseudogenes") %>% 
      na.omit() %>% 
  
  # Plot
  ggplot(aes(x = no_tissues_expressed))+
  geom_histogram(aes(y=..density..),
                 position = "dodge", 
                 colour = "black", 
                 binwidth = 2,
                 fill = "lightblue") + 
  labs(x = "Number of tissues expressed (GTEx)", y = "Percentage of processed and unprocessed\n pseudogenes expressed") +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%"),
                     n.breaks = 5,
                     limits = c(0, 0.2)) +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_blank())
  

# Save data ---------------------------------------------------------------

ggsave(plot = pseudogene_expression_to_plot, 
       filename = "pseudo_expression_plot.svg", 
       path = here::here("results", "pseudogenes"), 
       width = 7, 
       height = 7, 
       dpi = 600
)
