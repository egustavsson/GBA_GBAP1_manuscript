# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggforce)

# Load data ---------------------------------------------------------------

sqanti_class <-
  list(
    GBA = read.table("/home/egust/Projects/GBA_GBAP1/raw_data/GBA_classification_processed.txt", header = TRUE, sep="\t"),
    GBAP1 = read.table("/home/egust/Projects/GBA_GBAP1/raw_data/GBAP1_classification_processed.txt", header = TRUE, sep="\t")
  )

# Main --------------------------------------------------------------------

######### Format data per gene #########

number_of_transcripts_per_gene_plot <-
  lapply(sqanti_class, 
         
         function(x){
           
           n_samples <- x %>%
             dplyr::select(starts_with("FL.")) %>% 
             ncol()
           
           counts.long <- x %>%
             dplyr::select(isoform, 
                           starts_with("NFLR")) %>%
             pivot_longer(cols = c(-isoform), 
                          names_to = "sample", 
                          values_to = "read_count")
           
           counts.long$sample <- counts.long$sample %>% 
             str_replace("NFLR.*\\.", "") %>% 
             str_replace("_3p", "")
           
           threshs <- seq(from = 0, to = max(x$NFLR_mean), by = 0.1)
           
           res <- data.frame()
           
           for(i in threshs){
             transcripts_per_sample <- counts.long %>%
               dplyr::filter(read_count >=i) %>%
               dplyr::group_by(sample) %>%
               dplyr::summarise(n_transcripts = dplyr::n(), 
                         .groups = "keep") %>%
               dplyr::mutate(NFLR_threshold = i)
             
             res <- rbind(res, transcripts_per_sample)
             
             rm(transcripts_per_sample)
             
           }
           
           df <- res %>%
             dplyr::filter(sample %in% "NFLR_mean")
           
           sd <- res %>%
             dplyr::filter(!sample %in% "NFLR_mean") %>%
             dplyr::group_by(NFLR_threshold) %>%
             dplyr::summarise(sd = sd(n_transcripts), 
                       .groups = "keep")
           
           df <- left_join(df, 
                           sd, 
                           by = "NFLR_threshold")
           
    return(df)
           
         }
  
) %>% 
  
  bind_rows(.id = "gene")

######### Plot data per gene #########

for (i in unique(number_of_transcripts_per_gene_plot$gene)) {
  
  number_of_transcripts_per_gene_plot %>% dplyr::filter(gene == i) %>%
    
    ggplot(aes(x=NFLR_threshold, 
               y = n_transcripts)) +
    geom_line(group=1) +
    geom_ribbon(aes(NFLR_threshold, 
                    ymax = n_transcripts + sd, 
                    ymin = n_transcripts - sd), 
                alpha = 0.2, 
                fill = "blue", group=1) +
    ggtitle(i) +
    scale_x_continuous(name = expression(paste("Normalised expression (", NFLR[T], ")"))) +
    scale_y_continuous(name = "No. of transcripts detected") +
    guides(fill="none") +
    facet_zoom(xlim = c(0, 5),
             ylim = c(0, 50),
             horizontal = F) +
    
    theme(plot.title = element_text(face = "bold",
                                    size = 16,
                                    hjust = 0.5),
          panel.border = element_blank(),
          axis.line = element_line(color = "black", size=0.4),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(face = "bold",
                                     size = 10),
          axis.text.y = element_text(face = "bold",
                                     size = 10),
          axis.title.y = element_text(size = 14),
          strip.text.x = element_text(face = "bold",
                                      size = 12),
          legend.title = element_blank())
  
  ggsave(filename = paste0(i, "_transcripts_detected_plot.svg"), 
         path = here::here("results", "transcripts"), 
         width = 6, 
         height = 4, 
         dpi = 600
  )
  
}  
