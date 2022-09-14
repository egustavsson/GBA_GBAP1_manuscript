# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)

# Load data ---------------------------------------------------------------

gtex_filtered <- 
  readRDS(here::here("results", "GTEx", "gtexv8_GBA_GBAP1_tpm.Rds")) %>% 
  dplyr::mutate(Organ = gsub("Brain - ", "" , Organ))

# Main --------------------------------------------------------------------

## Plot GBA and GBAP1 expression in all tissues ##

# TPM

GBA_GBA_GTEx_expression_to_plot <- 
  gtex_filtered %>%
  dplyr::select(Description, 
                tpm, 
                Organ,
                `Brain tissue?`) %>%
  dplyr::group_by(Description, 
                  Organ) %>% 
  ggplot(aes(x = Description, 
             y = tpm)) +
  geom_violin(aes(fill = Description), 
              scale = TRUE, 
              show.legend = F) +
  geom_boxplot(width=0.05, 
               outlier.shape = NA, 
               show.legend = F) +
  scale_fill_manual(values = c("#5ab4ac", "#ffffff"),
                    breaks = c("GBA", "GBAP1")) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             colour = "Black", 
                             size = 2, vjust = 0.7) +
  facet_wrap(. ~Organ, 
             ncol = 6) +
  labs(x = "", 
       y = "TPM") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(face = "bold",
                                   size = 10),
        axis.text.y = element_text(face = "bold",
                                   size = 10),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(face = "bold",
                                    size = 6),
        strip.background =element_rect(fill="gray83"),
        legend.title = element_blank())

# Fold change

GBA_GBA_GTEx_fold_change_to_plot <-
  gtex_filtered %>% 
  split(.$Description) %>% 
  purrr::reduce(left_join, by = "SAMPID") %>% 
  dplyr::mutate(FC = log2(tpm.x / tpm.y)) %>% 
  dplyr::select(FC, 
                Organ.x,
                `Brain tissue?.x`) %>% 
  dplyr::group_by(Organ.x) %>% 
  
  ggplot(aes(x = FC, 
             y = reorder(Organ.x, -FC))) +
  geom_violin(aes(fill = `Brain tissue?.x`),
              scale = TRUE, 
              show.legend = T) +
  geom_boxplot(width=0.1, 
               outlier.shape = NA) +
  scale_fill_manual("Brain tissue?",
                    breaks = c("TRUE", "FALSE"),
                    values = c("#5ab4ac", "#ffffff")) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", colour = "black") +
  labs(y = "", 
       x = bquote(log[2]~fold~change~(GBA/GBAP1))) +
  xlim(-1, 7) +
  theme_classic() +
  theme(legend.position = "top",
        text = element_text(size = 12),
        strip.text.x = element_text(face = "bold",
                                    size = 12),
        strip.background =element_rect(fill = "gray83"),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "bold",
                                   size = 10),
        axis.title.y = element_text(size = 14))

## Plot GBA and GBAP1 expression in select tissues ##

# TPM

select_tissue_expression_to_plot <- 
  gtex_filtered %>%
  dplyr::select(Description, 
                tpm, 
                Organ,
                `Brain tissue?`) %>%
  dplyr::group_by(Description, 
                  Organ) %>% 
  dplyr::filter(Organ %in% c("Frontal Cortex (BA9)", "Hippocampus")) %>%
  
  ggplot(aes(x = Description, 
             y = tpm)) +
  geom_violin(aes(fill = Description),
              scale = TRUE, 
              show.legend = F) +
  geom_boxplot(width=0.2, 
               outlier.shape = NA, 
               show.legend = F) +
  scale_fill_manual(values = c("#5ab4ac", "#ffffff"), 
                    breaks = c("GBA", "GBAP1")) +
  ggpubr::stat_compare_means(method = "wilcox.test", 
                             colour = "Black", 
                             size = 6) +
  facet_wrap(. ~Organ, 
             ncol = 4) +
  labs(x = "", 
       y = "Short-read RNA seq expression (TPM)") +
  theme_classic() +
  theme(axis.title = element_text(size = 22),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18,
                                   face = "bold.italic"),
        strip.text.x = element_text(face = "bold",
                                    size = 18),
        strip.background = element_rect(fill="gray83"),
        legend.title = element_blank())


# Save data ---------------------------------------------------------------

# All tissues
ggsave(plot = GBA_GBA_GTEx_expression_to_plot, 
      filename = "GBA_GBA_GTEx_expression_to_plot.svg", 
      path = here::here("results", "GTEx"), 
      width = 11, 
      height = 10, 
      dpi = 600
)

ggsave(plot = GBA_GBA_GTEx_fold_change_to_plot, 
       filename = "GBA_GBA_GTEx_fold_change_to_plot.svg", 
       path = here::here("results", "GTEx"), 
       width = 8, 
       height = 8, 
       dpi = 600
)

# Select tissues
ggsave(plot = select_tissue_expression_to_plot, 
       filename = "select_tissue_expression_to_plot.svg", 
       path = here::here("results", "GTEx"), 
       width = 10, 
       height = 6, 
       dpi = 600
)
