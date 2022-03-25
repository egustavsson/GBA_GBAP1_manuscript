# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)

# Load data ---------------------------------------------------------------

Transcripts <-
  list(
    GBA = read_tsv("/home/egust/Projects/GBA_GBAP1/raw_data/GBA_classification_filtered.txt", show_col_types = F),
    GBAP1 = read_tsv("/home/egust/Projects/GBA_GBAP1/raw_data/GBAP1_classification_filtered.txt", show_col_types = F)
  )

# Main --------------------------------------------------------------------

## Data to plot ##

Transcripts_to_plot <-
  bind_rows(!!!Transcripts) %>% 
  dplyr::select(Isoform_class, 
                associated_gene, 
                starts_with("NFLR.")) %>% 
  rename_with(~paste0("", gsub("NFLR.Clontech_5p..|_3p", "", .)),
              starts_with('NFLR.Clontech_5p..')) %>% 
  
  pivot_longer(!c(Isoform_class, associated_gene), 
               names_to = "region", 
               values_to = "count") %>% 
  aggregate(count ~ Isoform_class + region + associated_gene,
            data = .,
            FUN = "sum") %>% 
  dplyr::mutate(region = str_replace_all(region, "_", " "),
                associated_gene = str_replace_all(associated_gene, c("ENSG00000177628.16" = "GBA",
                                                                     "ENSG00000160766.14" = "GBAP1"))) 

## Function to define order ##

plot_order <- function(data, category, gene){
  
    data %>%
    group_by(associated_gene, region) %>% 
    dplyr::mutate(freq = count / sum(count)) %>% 
    dplyr::filter(associated_gene == gene,
                  Isoform_class == category) %>%
    arrange(freq) %>% 
    pull(region)
}

## Plot ##

Transcripts_per_tissue_plot <-
  Transcripts_to_plot %>% 
  
  ggplot(aes(x = factor(region,
                        plot_order(data = Transcripts_to_plot,
                                   category = "Coding known (alternate 3/5 end)",
                                   gene = "GBA")),
             y = count,
             fill = Isoform_class)) +
  geom_col(colour = "black", 
           position = "fill") +
  labs(x = "Tissue",
       y = "Expression per transcript category") +
  facet_wrap(~factor(associated_gene)) +
  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(face = "bold",
                                   size = 8),
        axis.text.y = element_text(face = "bold",
                                   size = 12),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(face = "bold",
                                    size = 12),
        strip.background =element_rect(fill="grey"),
        legend.title = element_blank())

# Save data ---------------------------------------------------------------

ggsave(plot = Transcripts_per_tissue_plot, 
       filename = "Transcripts_per_tissue.png", 
       path = here::here("results", "transcripts"), 
       width = 10, 
       height = 4, 
       dpi = 600
)
