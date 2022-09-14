# Load libraries ----------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
library(plyranges)
library(here)
library(rtracklayer)
library(patchwork)

# Load data ---------------------------------------------------------------

load(here::here("results", "pseudogenes", "pseudogene_annotated.rda"))

# Functions ---------------------------------------------------------------

is.brain <- 
  function(x) {
    
    grepl("Brain", x)
  }

pseudogene_expression <-
  
  function(pseudogenes) {
    
    pseudo_tidy <-
      pseudogenes %>%
      dplyr::filter(pseudogene_type %in% c("Unprocessed", "Processed")) %>%
      dplyr::select(pseudogene_type, 34:86) %>% 
      dplyr::select(-c("Bladder",
                       "Breast - Mammary Tissue",
                       "Cells - Cultured fibroblasts",
                       "Cells - EBV-transformed lymphocytes",
                       "Cervix - Endocervix",
                       "Fallopian Tube",
                       "Ovary",
                       "Prostate",
                       "Testis",
                       "Uterus",
                       "Vagina",
                       "Brain - Cerebellum",
                       "Whole Blood",
                       "Skin - Sun Exposed (Lower leg)",
                       "Skin - Not Sun Exposed (Suprapubic)",
                       "Brain - Cerebellar Hemisphere")) %>% 
      na.omit() %>% 
      dplyr::mutate_if(is.numeric, ~1 * (. > 0)) %>% 
      dplyr::group_by(pseudogene_type) %>% 
      dplyr::summarise_if(is.numeric, sum) %>%
      pivot_longer(!pseudogene_type,
                   values_to = "Pseudogenes",
                   names_to = "Tissue") %>% 
      dplyr::mutate(Pseudogenes = ifelse(pseudogene_type == "Processed",
                                         Pseudogenes / sum(pseudogene_annotated$pseudogene_type == "Processed"),
                                         Pseudogenes / sum(pseudogene_annotated$pseudogene_type == "Unprocessed")))
    
    return(pseudo_tidy)
  }

Mean_expression_per_tissue <-
  
  function(data) {
    
    means <-
      dplyr::left_join(data[data$pseudogene_type == "Unprocessed", ],
                       data[data$pseudogene_type == "Processed", ],
                       by = c("Tissue" = "Tissue")) %>% 
    
      dplyr::mutate(exp = rowMeans(.[, c("Pseudogenes.x", "Pseudogenes.y")]) *100,
                    Brain = sapply(Tissue, is.brain))
    
    return(means)
  }

# Main --------------------------------------------------------------------

pseudo_exp <- 
  pseudogene_expression(pseudogenes = pseudogene_annotated)

# Pseudogenes expressed per tissue

pseudogene_expression_per_tissue <-
  Mean_expression_per_tissue(data = pseudo_exp) %>% 
  dplyr::filter(Tissue != "MeanTPM") %>% 
  
  ggplot(aes(x = reorder(Tissue, -exp), y = exp)) +
  geom_col(fill = "grey83", color = "Black", size = 0.5, width = 0.7) +
  scale_colour_manual(labels = c("Non-brain", "Brain"),
                      values = c("grey", "darkturquoise")) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "",
       y = "Pseudogenes expressed (%)") +
  theme_classic() +
  coord_flip() +
  theme(legend.title=element_blank(),
        legend.position = c(0.8, 0.5),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 12,
                                 face = "bold"),
        axis.text.x = element_text(size = 12,
                                   face = "bold"))

# Save data ---------------------------------------------------------------

ggsave(
  plot = pseudogene_expression_per_tissue, 
  filename = "pseudogene_expression_per_tissue_plot.svg", 
  path = here::here("results", "pseudogenes"), 
  width = 9, 
  height = 7, 
  dpi = 600
)
