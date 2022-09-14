# Load libraries ----------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
library(plyranges)
library(here)
library(ggrepel)

# Load data ---------------------------------------------------------------

load(here::here("results", "pseudogenes", "pseudogene_annotated.rda"))

# Main --------------------------------------------------------------------

GENCODEv38_pseudogenes_to_plot <-
  pseudogene_annotated %>% 
  dplyr::select(gene_id, pseudogene_type) %>% 
  unique() %>% # have to unique() since some genes are duplicated when left_join() with GTEx data. This is due to PAR1 region.
  group_by(pseudogene_type) %>% 
  tally() %>% 
  dplyr::mutate(csum = rev(cumsum(rev(n))),
                pos = n/2 + lead(csum, 1),
                pos = ifelse(is.na(pos), n/2, pos)) %>% 
  
  ggplot(aes(x = "",
             y = n,
             fill = pseudogene_type
  )) +
  geom_bar(width = 1,
           stat = "identity",
           color = "black") +
  geom_label_repel(aes(label = n, 
                       y = pos, size=3), 
                   show.legend = F, 
                   nudge_x = 1,
                   size = 6)+
  coord_polar("y", start = 0) +
  scale_fill_manual("Human pseudogenes", 
                    values = c("#7570B3", "#E6AB02", "#FFFFFF"),
                    breaks = c("Processed", "Unprocessed", "Other")) +
  theme_classic() +
  theme(legend.position="top", 
        legend.justification="top",
        legend.direction="horizontal",
        legend.title=element_text(size = 16,
                                  face = "bold"),
        legend.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank()) + 
  guides(fill = guide_legend(nrow = 1, 
                             byrow = T, 
                             title.position = "top",
                             title.hjust = 0.5))

# Save data ---------------------------------------------------------------

ggsave(plot = GENCODEv38_pseudogenes_to_plot, 
       filename = "pseudogenes_per_type.svg", 
       path = here::here("results", "pseudogenes"), 
       width = 4, 
       height = 4, 
       dpi = 600
)
