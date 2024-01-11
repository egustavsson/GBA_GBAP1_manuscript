# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)

# Load data ---------------------------------------------------------------

# Parent genes
parentgenes <- readRDS(here::here("results", "pseudogenes", "parentgene_annotated.rds"))

# Main --------------------------------------------------------------------

parent_genes_OMIM_plot <-
  parentgenes %>% 
  dplyr::filter(associated_pseudo %in% c("Unprocessed", "Processed")) %>% 
  dplyr::select(omim_morbid) %>% 
  table() %>%
  prop.table() %>% 
  data.frame() %>% 
  dplyr::rename("OMIM morbid" = ".") %>% 
  
  ggplot(aes(x = "", y = Freq, fill = `OMIM morbid`)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  scale_fill_manual(values=c("white", "black"),
                    labels = c("Not in OMIM morbid", "In OMIM morbid"),
                    name = "Parent genes") +
  coord_polar("y", start = 20)+
  geom_text(aes(label = paste0(round(Freq*100), "%")), position = position_stack(vjust = 0.5), color = c("black", "white"), size = 8) +
  
  theme_void() +
  theme(legend.position="top", 
        legend.justification="top",
        legend.direction="horizontal",
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_legend(nrow = 1, 
                             byrow = T, 
                             title.position = "top",
                             title.hjust = 0.5))

# Save data ---------------------------------------------------------------

ggsave(plot = parent_genes_OMIM_plot, 
       filename = "parent_genes_OMIM_plot.svg", 
       path = here::here("results", "pseudogenes"), 
       width = 6, 
       height = 5, 
       dpi = 600, 
       bg = "white"
)

