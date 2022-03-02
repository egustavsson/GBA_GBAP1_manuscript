# Load libraries ----------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
library(plyranges)
library(here)
library(ggrepel)

# Load data ---------------------------------------------------------------

GENCODEv38 <- rtracklayer::import("/home/egust/Projects/pseudogenes/raw_data/GENCODE_annotations/gencode.v38.annotation.gtf.gz")

# Main --------------------------------------------------------------------

GENCODEv38_pseudogenes_to_plot <-
  GENCODEv38 %>% 
  plyranges::filter(str_detect(gene_type, "pseudogene") & type == "gene") %>%
  plyranges::mutate(pseudogene_type = case_when(gene_type %in% c("IG_V_pseudogene", 
                                                                 "IG_C_pseudogene",
                                                                 "IG_J_pseudogene",
                                                                 "IG_pseudogene") ~ "IG",
                                                gene_type %in% c("TR_J_pseudogene",
                                                                 "TR_V_pseudogene") ~ "TR",
                                                gene_type %in% c("transcribed_unprocessed_pseudogene",
                                                                 "unprocessed_pseudogene",
                                                                 "translated_unprocessed_pseudogene") ~ "Unprocessed",
                                                gene_type %in% c("processed_pseudogene",
                                                                 "transcribed_processed_pseudogene",
                                                                 "translated_processed_pseudogene") ~ "Processed",
                                                gene_type %in% c("transcribed_unitary_pseudogene",
                                                                 "unitary_pseudogene") ~ "Unitary",
                                                gene_type %in% c("polymorphic_pseudogene",
                                                                 "rRNA_pseudogene",
                                                                 "pseudogene") ~ "Other")) %>% 
  data.frame() %>% 
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
                   nudge_x = 1)+
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  theme(legend.position="top", 
        legend.justification="top",
        legend.direction="horizontal",
        legend.title=element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),
        legend.background = element_rect(fill="gray93",
                                         size=0.8, linetype="solid", 
                                         colour ="black")) + 
  guides(fill = guide_legend(nrow = 2, 
                             byrow = T))

# Save data ---------------------------------------------------------------

ggsave(plot = GENCODEv38_pseudogenes_to_plot, 
       filename = "pseudogenes_per_type.png", 
       path = here::here("results", "pseudogenes"), 
       width = 4, 
       height = 4, 
       dpi = 600
)
