# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(rstatix)

# Load data ---------------------------------------------------------------

GCase_activity <-
  read_tsv(
    here::here("results", "biochem", "GCase_activity_formatted.txt"), show_col_types = F)

# Main --------------------------------------------------------------------

## Rename constructs ##

GCase_activity_to_plot <-
  GCase_activity %>%
  na.omit() %>%
  dplyr::mutate(construct = case_when(construct == "pcDNA3.1(+)" ~ "pcDNA3.1\n(empty vector)",
                                      construct == "Construct 1" ~ "PB.845.525\n(GBAP1)",
                                      construct == "Construct 2" ~ "PB.845.2627\n(GBA affecting GH30 and SP)",
                                      construct == "Construct 3" ~ "PB.845.2629\n(GBA affecting GH30 and SP)",
                                      construct == "Construct 4" ~ "PB.845.1693\n(GBAP1)",
                                      construct == "Construct 5" ~ "ENST00000368373.8\n(GBA MANE select)",
                                      construct == "Construct 6" ~ "PB.845.2954\n(GBA affecting SP)"))

## Define order of columns ##

GCase_activity_to_plot$construct <- 
  factor(GCase_activity_to_plot$construct,
         levels = c("pcDNA3.1\n(empty vector)",
                    "ENST00000368373.8\n(GBA MANE select)",
                    "PB.845.2954\n(GBA affecting SP)",
                    "PB.845.2627\n(GBA affecting GH30 and SP)",
                    "PB.845.2629\n(GBA affecting GH30 and SP)",
                    "PB.845.525\n(GBAP1)",
                    "PB.845.1693\n(GBAP1)"))

GCase_activity_to_plot$Genotype <- 
  factor(GCase_activity_to_plot$Genotype,
         levels = c("H4 Parental",
                    "H4 GBA KO"))

## Tuke post hoc ##

tukey_data <- 
  GCase_activity_to_plot %>% 
  rstatix::group_by(Genotype) %>% 
  tukey_hsd(`Enzyme activity` ~ construct) %>% 
  dplyr::filter(group1 == "pcDNA3.1\n(empty vector)" | group2 == "pcDNA3.1\n(empty vector)")

## ggboxplot ##

Enzyme_activity_plot <-
  
  ggpubr::ggbarplot(GCase_activity_to_plot, 
                  x = "construct",
                  y = "Enzyme activity",
                  fill = "construct", 
                  add = c("mean_sd", "jitter"),
                  facet.by = c("Genotype"),
                  width = 0.9) +
  ggpubr::stat_pvalue_manual(tukey_data,
                             hide.ns = TRUE,
                             y.position = 9,
                             size = 5,
                             bracket.size = 0.5,
                             tip.length = 0.01) +
  ylim(0, 10) +
  labs(x = "",
       y = "Enzyme activity (normalised to UTC)") +
  scale_fill_brewer(palette = "Dark2") +
  
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text.x = element_text(face = "bold",
                                   size = 10, angle = 45, 
                                   hjust = 1),
        axis.text.y = element_text(face = "bold",
                                   size = 10),
        axis.title.y = element_text(size = 14,
                                    face = "bold"),
        strip.text.x = element_text(face = "bold",
                                    size = 12),
        strip.background =element_rect(fill = "gray80"))

# Save data ---------------------------------------------------------------

ggsave(plot = Enzyme_activity_plot, 
       filename = "Enzyme_activity_plot.png", 
       path = here::here("results", "biochem"), 
       width = 8, 
       height = 6, 
       dpi = 600
)
