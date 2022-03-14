# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)

# Load data ---------------------------------------------------------------

gtex_filtered <- 
  readRDS("/home/egust/Projects/GBA_GBAP1/results/gtexv8_GBA_GBAP1_tpm.Rds") %>% 
  dplyr::mutate(Organ = gsub("Brain - ", "" , Organ))

# Main --------------------------------------------------------------------

## Plot GBA and GBAP1 TPM ##

Expression_to_plot <- 
  gtex_filtered %>%
  dplyr::select(Description, 
                tpm, 
                Organ,
                `Brain tissue?`) %>%
  dplyr::group_by(Description, 
                  Organ) %>% 
  dplyr::filter(Organ %in% c("Frontal Cortex (BA9)", "Anterior cingulate cortex (BA24)")) %>%
  
  ggplot(aes(x = Description, 
             y = tpm)) +
  geom_violin(aes(fill = `Brain tissue?`), 
              scale = TRUE, 
              show.legend = F) +
  geom_boxplot(width=0.05, 
               outlier.shape = NA, 
               show.legend = F) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             colour = "Black", 
                             size = 4) +
  facet_wrap(. ~Organ, 
             ncol = 4) +
  labs(x = "", 
       y = "TPM") +
  scale_fill_manual(values = c("#888888", "#00BFC4")) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(face = "bold",
                                   size = 10),
        axis.text.y = element_text(face = "bold",
                                   size = 12),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(face = "bold",
                                    size = 12),
        strip.background =element_rect(fill="gray83"),
        legend.title = element_blank())

## Plot fold change ##




Fold_change_to_plot <-
  gtex_filtered %>% 
  split(.$Description) %>% 
  purrr::reduce(left_join, by = "SAMPID") %>% 
  dplyr::mutate(FC = log2(tpm.x / tpm.y)) %>% 
  dplyr::select(FC, 
                Organ.x) %>% 
  dplyr::filter(Organ.x %in% c("Frontal Cortex (BA9)", "Anterior cingulate cortex (BA24)")) %>% 
  dplyr::group_by(Organ.x) %>% 
  
  ggplot(aes(x = FC, 
             y = reorder(Organ.x, -FC))) +
  geom_violin(aes(fill = "gray83"),
              scale = TRUE, 
              show.legend = F) +
  geom_boxplot(width=0.1, 
               outlier.shape = NA) +
  #geom_vline(aes(xintercept = mean(delta)), linetype = "dashed", colour = "red") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", colour = "grey") +
  labs(y = "", 
       x = bquote(log[2]~fold~change~(GBA/GBAP1))) +
  #xlim(-2, 10) +
  scale_fill_manual(values = c("#888888", "#00BFC4")) +
  theme_classic() +
  coord_flip() +
  facet_wrap(. ~Organ.x, scales = "free_x") +
  theme(legend.position = "top", 
        text = element_text(size = 12),
        strip.text.x = element_text(face = "bold",
                                    size = 12),
        strip.background =element_rect(fill = "gray83"),
        axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "bold",
                                   size = 10),
        axis.title.y = element_text(size = 14))
