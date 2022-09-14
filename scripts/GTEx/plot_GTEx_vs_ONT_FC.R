# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)

# Load data ---------------------------------------------------------------

# GTEx gene level expression data
gtex_filtered <- 
  readRDS("/home/egust/Projects/GBA_GBAP1/results/gtexv8_GBA_GBAP1_tpm.Rds") %>% 
  dplyr::mutate(Organ = gsub("Brain - ", "" , Organ))

# ONT direct cDNA sequencing data from Frontal lobe and Hippocampus
ONT <- 
  list(FL = read_tsv(here::here("data", "ONT_long_read", "gene_abundance_FL.tab"), show_col_types = FALSE) %>% 
         dplyr::mutate(Organ = "Frontal lobe"),
       
       Hipp = read_tsv(here::here("data", "ONT_long_read", "gene_abundance_Hipp.tab"), show_col_types = FALSE) %>% 
         dplyr::mutate(Organ = "Hippocampus"))

# Main --------------------------------------------------------------------

# Calculate log2 FC from ONT data
LR_FC <-
  lapply(ONT, function(x){
    
    data.frame(Organ = unique(x$Organ),
               FC = log2(x$TPM[x$`Gene Name` == "GBA"] / x$TPM[x$`Gene Name` == "GBAP1"]))
  }) %>% 
  bind_rows(!!!.)
 
# Change to matching tissue name
LR_FC[LR_FC$Organ == "Frontal lobe", ]$Organ = "Frontal Cortex (BA9)"

# Calculate log2 FC from GTEx data
GTEx_GBA_GBAP1 <-
  gtex_filtered %>% 
  dplyr::select(SAMPID, Description, Organ, tpm) %>% 
  split(., list(.$Description))
  
GTEx_FC <-
  dplyr::left_join(GTEx_GBA_GBAP1$GBA,
                   GTEx_GBA_GBAP1$GBAP1,
                   by = "SAMPID") %>% 
  dplyr::mutate(FC = log2(tpm.x / tpm.y)) %>% 
  dplyr::select(Organ = Organ.x, FC) %>% 
  dplyr::filter(Organ %in% c("Frontal Cortex (BA9)", "Hippocampus")) 


# add tha stats and a dotted line to plot #

# Frontal lobe
grubbs_frontal <-
  GTEx_FC %>% 
  dplyr::filter(Organ == "Frontal Cortex (BA9)") %>% 
  bind_rows(., LR_FC[LR_FC$Organ == "Frontal Cortex (BA9)", ]) %>% .$FC %>% 
  outliers::grubbs.test(x = ., two.sided = F)

# Hippocampus
grubbs_hipp <-
  GTEx_FC %>% 
  dplyr::filter(Organ == "Hippocampus") %>% 
  bind_rows(., LR_FC[LR_FC$Organ == "Hippocampus", ]) %>% .$FC %>% 
  outliers::grubbs.test(x = ., two.sided = F)

stats <- data.frame(Organ = c("Frontal Cortex (BA9)", "Hippocampus"),
                    p_val = c(paste0("P = ", round(grubbs_frontal$p.value, 3)),
                              paste0("P = ", round(grubbs_hipp$p.value, 3))))

# Plot
expression_GTEx_ONT_plot <-
  GTEx_FC %>% 
  dplyr::group_by(Organ) %>% 
  ggplot(aes(x = FC)) +
  geom_density(fill = "gray63", show.legend = F, adjust = 0.8, size = 1) +
  
  # GTEx mean FC
  geom_vline(data = filter(GTEx_FC, Organ == "Frontal Cortex (BA9)"), aes(xintercept = mean(FC)), size = 1, linetype = "dashed") +
  geom_vline(data = filter(GTEx_FC, Organ == "Hippocampus"), aes(xintercept = mean(FC)), size = 1, linetype = "dashed") +
  
  # ONT FC
  geom_vline(data = LR_FC[LR_FC$Organ == "Frontal Cortex (BA9)", ], aes(xintercept = FC), size = 1, linetype = "dashed", colour = "red") +
  geom_vline(data = LR_FC[LR_FC$Organ == "Hippocampus", ], aes(xintercept = FC), size = 1, linetype = "dashed", colour = "red") +
  
  # grubbs p value
  geom_text(data = filter(stats, Organ == "Frontal Cortex (BA9)"),
            label = stats[stats$Organ == "Frontal Cortex (BA9)", ]$p_val, size = 8, aes(x = 1.5, y = 0.75)) +
  geom_text(data = filter(stats, Organ == "Hippocampus"),
            label = stats[stats$Organ == "Hippocampus", ]$p_val, size = 8, aes(x = 1.75, y = 0.75)) +
  
  labs(y = "Density", x = bquote(log[2]~fold~change~(GBA/GBAP1))) +
  facet_wrap(. ~Organ) +
  theme_classic() +
  theme(legend.position = "top", 
        text = element_text(size = 14),
        strip.text.x = element_text(face = "bold",
                                    size = 18),
        strip.background = element_rect(fill="gray83"),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))

# Save data ---------------------------------------------------------------

ggsave(plot = expression_GTEx_ONT_plot, 
       filename = "expression_GTEx_ONT_plot.svg", 
       path = here::here("results", "GTEx"), 
       width = 10, 
       height = 6, 
       dpi = 600
)
