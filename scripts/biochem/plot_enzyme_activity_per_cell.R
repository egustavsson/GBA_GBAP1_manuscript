# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(rstatix)
library(ggpubr)

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

GCase_activity_to_plot$`Cell line` <- 
  factor(GCase_activity_to_plot$`Cell line`,
         levels = c("H4 GBA KO",
                    "H4 Parental"))

## Tuke post hoc ##

tukey_data <- 
  GCase_activity_to_plot %>% 
  rstatix::group_by(`Cell line`) %>% 
  tukey_hsd(`Enzyme activity` ~ construct) %>% 
  add_y_position(scales = "free_y") %>%
  dplyr::filter(group1 == "pcDNA3.1\n(empty vector)" | group2 == "pcDNA3.1\n(empty vector)")

## ggboxplot ##

for (i in c("H4 GBA KO", "H4 Parental")) {
  
  data_to_plot <-
    list(activity = dplyr::filter(GCase_activity_to_plot, `Cell line` == i),
         tukey = dplyr::filter(tukey_data, `Cell line` == i))
  
  
  
  ggpubr::ggboxplot(data_to_plot$activity, 
                    x = "construct",
                    y = "Enzyme activity",
                    fill = "grey83", 
                    add = c("jitter"),
                    notch = F,
                    width = 0.9, title = i) +
    ggpubr::stat_pvalue_manual(data_to_plot$tukey,
                               hide.ns = TRUE,
                               size = 5,
                               bracket.size = 0.5,
                               tip.length = 0.01) +
    labs(x = "",
         y = "Enzyme activity (normalised to UTC)") +
    scale_fill_brewer(palette = "Dark2") +
    ylim(0, 9) +
    labs(title = i) +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(color="black", size=20, face="bold", vjust = 0.5, hjust = 0.5),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(face = "bold",
                                     size = 10, angle = 45, 
                                     hjust = 1),
          axis.text.y = element_text(face = "bold",
                                     size = 10),
          axis.title.y = element_text(size = 14,
                                      face = "bold"))
  
  # Save data ---------------------------------------------------------------
  
  ggsave(filename = paste0(i, "_Enzyme_activity_plot.svg"), 
         path = here::here("results", "biochem"), 
         width = 8, 
         height = 7, 
         dpi = 600
  )  
  
}

