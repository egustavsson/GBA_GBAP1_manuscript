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

Ratio_expression_per_tissue <-
  
  function(data) {
    
    ratio <-
      dplyr::left_join(data[data$pseudogene_type == "Unprocessed", ],
                       data[data$pseudogene_type == "Processed", ],
                       by = c("Tissue" = "Tissue")) %>% 
      dplyr::mutate(`Unprocessed/Processed` = log2(Pseudogenes.x/Pseudogenes.y),
                    Brain = sapply(Tissue, is.brain))
  }

# Main --------------------------------------------------------------------

pseudo_exp <- 
  pseudogene_expression(pseudogenes = pseudogene_annotated)

# Pseudogenes expressed per tissue

pseudogene_expression_per_tissue <-
  Mean_expression_per_tissue(data = pseudo_exp) %>% 
  
  ggplot(aes(x = reorder(Tissue, +exp), y = exp)) +
  geom_point(aes(color = Brain), size = 2.5) +
  scale_colour_manual(labels = c("Non-brain", "Brain"),
                      values = c("grey", "darkturquoise")) +
  coord_flip() +
  labs(x = "Tissue",
       y = "Pseudogenes expressed (%)") +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position = c(0.8, 0.5),
        legend.text = element_text(size = 10))

# Pseudogenes fractions expressed per tissue

ratio_pseudogene_expression_per_tissue <-
  Ratio_expression_per_tissue(data = pseudo_exp) %>% 
  
  ggplot(aes(x = reorder(Tissue, +`Unprocessed/Processed`), y = `Unprocessed/Processed`)) +
  geom_point(aes(color = Brain), size = 2.5) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             color = "lightgrey") +
  scale_colour_manual(labels = c("Non-brain", "Brain"),
                      values = c("grey", "darkturquoise")) +
  coord_flip() +
  labs(x = "Tissue",
       y = bquote(log[2]~fold~change~(Unprocessed/Processed))) +
  theme_classic() +
  theme(legend.title=element_blank(),
        legend.position = c(0.8, 0.5),
        legend.text = element_text(size = 10))

# Save data ---------------------------------------------------------------

ggsave(
  plot = pseudogene_expression_per_tissue / ratio_pseudogene_expression_per_tissue, 
  filename = "pseudogene_expression_per_tissue_plot.png", 
  path = here::here("results", "pseudogenes"), 
  width = 7, 
  height = 9, 
  dpi = 600
)
