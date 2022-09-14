# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggpubr)

# Load data ---------------------------------------------------------------

Transcripts <-
  list(
    GBA = read_tsv(here::here("data", "Brain", "GBA_classification_filtered.txt"), show_col_types = F),
    GBAP1 = read_tsv(here::here("data", "Brain", "GBAP1_classification_filtered.txt"), show_col_types = F)
  )

# Functions ---------------------------------------------------------------

plot_transcripts_per_gene <-
  function(data, gene, gene_name, labelling) {
    
    # fill colour to use
    fill_colour <- c("Coding known (alternate 3/5 end)" = "#74add1",
                     "Coding novel" = "#4d9221",
                     "NMD novel" = "#d53e4f",
                     "Non-coding known" = "#b2abd2",
                     "Non-coding novel" = "#d8daeb")
    
    
    ## Plot number of transcripts per category ##
    Transcripts_per_category <-
      data[[gene_name]] %>% 
      dplyr::select(Isoform_class, 
                    associated_gene) %>% 
      dplyr::mutate(associated_gene = gene_name) %>% 
      dplyr::count(associated_gene, Isoform_class)
    
    # if categories are missing populate df
    if(length(c("Coding known (alternate 3/5 end)",
                 "Coding novel",
                 "NMD novel",
                 "Non-coding known",
                 "Non-coding novel")) != length(Transcripts_per_category$Isoform_class)) {
      
      missing <-
      data.frame(associated_gene = gene_name,
                 Isoform_class = setdiff(c("Coding known (alternate 3/5 end)",
                                           "Coding novel",
                                           "NMD novel",
                                           "Non-coding known",
                                           "Non-coding novel"), 
                                         Transcripts_per_category$Isoform_class),
                 n = 0)
      
      Transcripts_per_category <- bind_rows(Transcripts_per_category, missing)
    }
    
    
    Transcripts_per_category_plot <-
      Transcripts_per_category %>% 
      ggplot(aes(x = factor(Isoform_class, 
                            levels = c("Coding known (alternate 3/5 end)",
                                       "Coding novel",
                                       "NMD novel",
                                       "Non-coding known",
                                       "Non-coding novel")), 
                 y = n, 
                 fill = Isoform_class)) +
      geom_col(show.legend = F, colour = "Black") +
      scale_fill_manual(values = fill_colour) +
      scale_x_discrete(labels = c("Coding known (alternate 3/5 end)" = "Coding known\n(alternate 3/5 end)",
                                  "Coding known (complete match)" = "Coding known\n(complete match)")) +
      labs(y = "No. unique transcripts", x = "Transcript category") +
      theme_classic() +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.title = element_text(size = 14),
            axis.text.y = element_text(size = 12),
            axis.text.x = element_text(angle = 45, 
                                       hjust=1,
                                       size = 12))
    
    
    ## Plot expression per transcript ##
    Expression_per_transcript <-
      data[[gene_name]] %>%
      dplyr::select(isoform, Isoform_class, starts_with("NFLR.")) %>%
      pivot_longer(cols = starts_with("NFLR."), names_to = "Sample", values_to = "NFLR") %>%
      group_by(isoform, Isoform_class) %>%
      summarise(NFLR_mean = mean(NFLR), NFLR_sd = sd(NFLR), .groups = "keep") %>%
      arrange(desc(NFLR_mean)) %>%
      tibble::rowid_to_column(., "isoform_index")
    
    Expression_per_transcript_plot <-
      Expression_per_transcript %>% 
      ggplot(aes(x=isoform_index, fill = Isoform_class)) +
      geom_bar(aes(y = NFLR_mean), color=NA, size=0.3, width=0.8, stat="identity") +
      geom_errorbar(aes(ymin = NFLR_mean-NFLR_sd, ymax = NFLR_mean+NFLR_sd), width = 0.2) +
      scale_fill_manual(values = fill_colour) +
      scale_y_continuous(name = "Transcript expression\nacross samples (%)", limits = c(0, 48)) +
      scale_x_continuous(name = "Transcripts ranked by expression", limits = c(0, 48)) +
      guides(fill=guide_legend(title="Transcript category")) +
      theme_classic() +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            axis.line = element_line(color = "black", size=0.4))
    
    ## Plot expression per category ##
    Expression_per_category <-
      data[[gene_name]] %>%
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
    
    Expression_per_category_plot <-
      Expression_per_category %>% 
      ggplot(aes(x = region,
                 y = count,
                 fill = Isoform_class)) +
      geom_col(colour = "black", 
               position = "fill") +
      labs(x = "",
           y = "Expression per transcript category") +
      scale_fill_manual(values = fill_colour) +
      scale_y_continuous(labels = scales::percent) +
      coord_flip() +
      theme_classic() +
      theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.title = element_text(size = 16),
            axis.text.x = element_text(face = "bold",
                                       size = 12),
            axis.text.y = element_text(face = "bold",
                                       size = 12),
            axis.title.y = element_text(size = 14),
            legend.title = element_blank())
    
    plot <- ggpubr::ggarrange(Transcripts_per_category_plot, Expression_per_transcript_plot, Expression_per_category_plot, 
                      nrow = 1, 
                      common.legend = T, align = "h",
                      labels = labelling, 
                      font.label = list(size = 24),
                      widths = c(1, 1.5, 2))
    
    return(annotate_figure(plot, top = text_grob(paste0(gene_name, " transcripts"), 
                                             color = "black", face = "bold", size = 20))) # return plots as ggarrange
  }

# Main --------------------------------------------------------------------

# Plot GBA

GBA_transcript_plot <-
  plot_transcripts_per_gene(data = Transcripts, gene = "ENSG00000177628.16", gene_name = "GBA", labelling = c("a", "b", "c"))

# Plot GBAP1

GBAP1_transcript_plot <-
  plot_transcripts_per_gene(data = Transcripts, gene = "ENSG00000160766.14", gene_name = "GBAP1", labelling = c("d", "e", "f"))


# Save data ---------------------------------------------------------------
ggsave(plot = GBA_transcript_plot, 
       filename = "GBA_transcript_plot.svg", 
       path = here::here("results", "transcripts"), 
       width = 16, 
       height = 6, 
       dpi = 600, 
       bg = "white"
)

ggsave(plot = GBAP1_transcript_plot, 
       filename = "GBAP1_transcript_plot.svg", 
       path = here::here("results", "transcripts"), 
       width = 16, 
       height = 6, 
       dpi = 600,
       bg = "white"
)
