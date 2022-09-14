# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(plyranges)
library(ggtranscript)
library(cowplot)

# Arguments ---------------------------------------------------------------

args <-
  list(
    snRNA_files = list.files(here::here("data", "5prime_CUTRUN"), 
                              full.names = TRUE, 
                              pattern = glob2rx("PFCTX*rev*")),
    CUTRUN_files = list.files(here::here("data", "5prime_CUTRUN"), 
                               full.names = TRUE, 
                               pattern = "^CR116")
  )

# Load data ---------------------------------------------------------------

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
  


# GBA and GBAP1 loci to use when filtering data
locus_subset <- 
  GRanges(seqnames = "chr1",
          ranges = IRanges(start = 155213078, 
                           end = 155245198),
          strand = "-")

# 5' single nuclear RNAseq data
snRNA_data <- lapply(args$snRNA_files, function(x) {
  
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

# CUT&RUN seq data
CUT_RUN <- lapply(args$CUTRUN_files, function(x) {
  
  rtracklayer::import.bw(x) %>%
    subsetByOverlaps(., locus_subset) %>% 
    plyranges::mutate(
      sample = stringi::stri_replace_all_regex(x,
                                               paste0(here::here("data", "5prime_CUTRUN"), "/"), 
                                                  "", 
                                                  vectorize_all = F),
      cell_type = case_when(grepl("IgG", sample) ~ "IgG",
                            grepl("H3K4me3", sample) ~ "H3K4me3")) %>% 
    data.frame()
  
  }
  ) %>% 
  bind_rows(!!!.)

# FANTOM5 CAGE seq data
CAGE_seq <- 
  rtracklayer::import.bed(here::here("data", "hg38.cage_peak_phase1and2combined_coord.bed.gz")) %>% 
  subsetByOverlaps(., locus_subset)

# Main --------------------------------------------------------------------

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
  annotate("text", 
           x=c(155217078,155237078), 
           y=c("PB.845.2848", "PB.845.1693"), 
           label=c("GBAP1", "GBA"), 
           size=4, 
           fontface="bold.italic") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        panel.border = element_rect(colour = "black", fill=NA, size=1))



# Plot CUT&RUN
CUT_RUN_plot <- 
  CUT_RUN %>% 
  ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = score)) +
  geom_rect(colour = "skyblue2", fill = "skyblue2", show.legend = F, alpha=0.8) +
  xlim(start(locus_subset), end(locus_subset)) +
  labs(y = "CUT&RUN read counts",
       x = "") +
  facet_wrap(vars(cell_type), 
             ncol = 1, 
             strip.position = "right") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10),
        strip.text.y = element_text(face = "bold",
                                    size = 10),
        strip.background = element_rect(fill ="lightgrey"),
        strip.placement = "outside",
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.spacing = unit(0, "cm"))

# Plot 5' single nuclear RNAseq
snRNA_plot <-
  snRNA_data %>% 
  ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = score)) +
  geom_rect(colour = "lightsalmon", fill = "lightsalmon", show.legend = F, alpha=0.8) +
  xlim(start(locus_subset), end(locus_subset)) +
  labs(y = "Expression levels (RPKM)",
       x = "") +
  facet_wrap(vars(cell_type), 
             ncol = 1, 
             scales = "free_y",
             strip.position = "right",
             labeller = label_wrap_gen(18)) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10),
        strip.text.y = element_text(face = "bold",
                                    size = 10),
        strip.background = element_rect(fill ="lightgrey"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.spacing = unit(0, "cm"))

# Plot CAGE seq coords (hg38)
CAGE_plot <-
  CAGE_seq %>% 
  data.frame() %>% 
  ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = 1)) +
  geom_rect(colour = "black", fill = "black", show.legend = F, alpha=0.8) +
  xlim(start(locus_subset), end(locus_subset)) +
  labs(y = "CAGE peaks",
       x = "") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 10,
                                   colour = "white"),
        axis.ticks.y = element_line(colour = "white"),
        strip.text.y = element_text(face = "bold",
                                    size = 10),
        strip.background = element_rect(fill ="lightgrey"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.spacing = unit(0, "cm"))

# Final plot
GBA_GBAP1_transcript_supportive_data_plot <-
  plot_grid(GBA_GBAP1_plot, 
            CUT_RUN_plot, 
            snRNA_plot, 
            CAGE_plot,
            ncol = 1, 
            align = "hv", 
            rel_heights = c(1, 1, 3.2, 0.5), 
            axis = "lr",
            labels = c("a", "b", "c", "d"), 
            label_size = 18)

# Save data ---------------------------------------------------------------

ggsave(plot = GBA_GBAP1_transcript_supportive_data_plot, 
       filename = "GBA_GBAP1_H3k_snRNAseq_plot.png", 
       path = here::here("results", "transcripts"), 
       width = 10, 
       height = 14, 
       dpi = 800
)
