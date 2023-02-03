# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(plyranges)
library(ggtranscript)
library(cowplot)

# Load data ---------------------------------------------------------------

# snRNAseq
snRNA_files = list.files(here::here("data", "5prime_CUTRUN"), 
                             full.names = TRUE, 
                             pattern = glob2rx("PFCTX*rev*"))

# Long-read transcripts
long_read <- rtracklayer::import(
  here::here("data", "Brain", "Brain_regions_corrected.gtf.cds.gff")) %>% 
  plyranges::filter(transcript_id %in% c(
    "PB.845.2888",
    "PB.845.2786", 
    "PB.845.2627", 
    "PB.845.2629",
    "PB.845.2848",
    "PB.845.2954",
    "PB.845.525", 
    "PB.845.1693",
    "PB.845.565"
  )) 


# Functions ---------------------------------------------------------------

plot_suppotive_data_per_gene <-
  function(seqnames, start, end, strand, gene_id) {
    
    # Filter long read gtf/gff for gene of interest
    long_read <- long_read[long_read$gene_id == gene_id]
    
    
    # loci used to filter data
    locus_subset <- 
      GRanges(seqnames = seqnames,
              ranges = IRanges(start = start, 
                               end = end),
              strand = strand)
    
    # 5' single nuclear RNAseq data
    snRNA_data <- lapply(snRNA_files, function(x) {
      
      rtracklayer::import.bw(x) %>%
        subsetByOverlaps(., locus_subset) %>%
        plyranges::mutate(
          sample = stringi::stri_replace_all_regex(x,
                                                   paste0(here::here("data", "5prime_CUTRUN"), "/"), 
                                                   "", 
                                                   vectorize_all = F),
          cell_type = case_when(grepl("Excitatoryneurons", sample) ~ "Excitatory neurons",
                                grepl("Inhibitoryneurons", sample) ~ "Inhibitory neurons",
                                grepl("Astrocytes", sample) ~ "Astrocytes",
                                grepl("Microglia", sample) ~ "Microglia",
                                grepl("OPC", sample) ~ "Oligodendrocyte progenitor cells",
                                grepl("Oligodendrocytes", sample) ~ "Oligodendrocytes")) %>% 
        data.frame()
      
    }
    ) %>% 
      bind_rows(!!!.)
    
    # Plot GBA and GBAP1 transcripts
    exons <- data.frame(long_read) %>% dplyr::filter(type == "exon")
    introns <- exons %>% to_intron(group_var = "transcript_id")
    CDS <- data.frame(long_read) %>% dplyr::filter(type == "CDS")
    
    GBA_GBAP1_plot <-
      exons %>%
      ggplot(
        aes(xstart = start, xend = end, y = factor(transcript_id, levels = c("PB.845.2888",
                                                                             "PB.845.2786", 
                                                                             "PB.845.2627", 
                                                                             "PB.845.2629",
                                                                             "PB.845.2848",
                                                                             "PB.845.2954",
                                                                             "PB.845.525", 
                                                                             "PB.845.1693",
                                                                             "PB.845.565")))) +
      geom_range(fill = "white",
                 height = 0.25) +
      geom_range(data = CDS) +
      geom_intron(data = introns, 
                  arrow.min.intron.length = 500, 
                  arrow = grid::arrow(ends = "first", length = grid::unit(0.1, "inches"))) +
      labs(y = "Transcript name",
           x = "") +
      xlim(start(locus_subset), end(locus_subset)) +
      theme_classic() +
      theme(axis.title = element_text(size = 14),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 10),
            axis.title.y = element_text(size = 14),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    
    # Plot 5' single nuclear RNAseq
    
    # Define colors per cell type
    fill_colors <- c("Excitatory neurons" = "#a6cee3", 
                     "Inhibitory neurons" = "#ab9c6d", 
                     "Astrocytes" = "#7fc564",
                     "Microglia" = "#fdbf6f",
                     "Oligodendrocyte progenitor cells" = "#3c8cab",
                     "Oligodendrocytes" = "#e73334")
    
    snRNA_plot <-
      snRNA_data %>% 
      ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = score, colour = cell_type, fill = cell_type)) +
      geom_rect(show.legend = F, alpha=0.8) +
      xlim(start(locus_subset), end(locus_subset)) +
      labs(y = "Expression levels (RPKM)",
           x = "") +
      scale_fill_manual(values = fill_colors) +
      scale_color_manual(values = fill_colors) +
      facet_wrap(vars(cell_type), 
                 ncol = 1, 
                 strip.position = "right",
                 labeller = label_wrap_gen(18)) +
      theme_classic() +
      theme(axis.title = element_text(size = 12),
            axis.text.x = element_text(size = 14),
            axis.ticks.x = element_line(),
            axis.text.y = element_text(size = 10),
            strip.text.y = element_text(face = "bold",
                                        size = 10),
            strip.background = element_rect(fill ="lightgrey"),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            panel.spacing = unit(0, "cm"))
    
    # Final plot
    GBA_GBAP1_transcript_supportive_data_plot <-
      plot_grid(GBA_GBAP1_plot, 
                snRNA_plot, 
                ncol = 1, 
                align = "hv", 
                rel_heights = c(1, 3.5), 
                axis = "lr")
    
    return(GBA_GBAP1_transcript_supportive_data_plot)
  }


# Main --------------------------------------------------------------------

# Plot GBA
GBA_snRNAseq_data_plot <-
  plot_suppotive_data_per_gene(seqnames = "chr1", start = 155233000, end = 155246000, strand = "-", gene_id = "ENSG00000177628.16")


# Plot GBAP1
GBAP1_snRNAseq_data_plot <-
  plot_suppotive_data_per_gene(seqnames = "chr1", start = 155213587, end = 155219638, strand = "-", gene_id = "ENSG00000160766.14")

# Save data ---------------------------------------------------------------

ggsave(plot = GBA_snRNAseq_data_plot, 
       filename = "GBA_snRNAseq_data_plot.svg", 
       path = here::here("results", "transcripts"), 
       width = 10, 
       height = 11, 
       dpi = 600
)

ggsave(plot = GBAP1_snRNAseq_data_plot, 
       filename = "GBAP1_snRNAseq_data_plot.svg", 
       path = here::here("results", "transcripts"), 
       width = 10, 
       height = 11, 
       dpi = 600
)
