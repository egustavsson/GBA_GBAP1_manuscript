# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggtranscript)
library(patchwork)
library(rtracklayer)
library(plyranges)
library(here)

# Load data ---------------------------------------------------------------

Filtered_tx <-
  list(
    GBA = read.table(here::here("data", "Brain", "GBA_classification_filtered.txt"), header = TRUE, sep="\t"),
    GBAP1 = read.table(here::here("data", "Brain", "GBAP1_classification_filtered.txt"), header = TRUE, sep="\t")
  )

LR_gff <- 
  list(
    GBA = rtracklayer::import(here::here("data", "Brain", "GBA_corrected.gtf.cds.gff")),
    GBAP1 = rtracklayer::import(here::here("data", "Brain", "GBAP1_corrected.gtf.cds.gff"))
  )

# Functions ---------------------------------------------------------------

get_TOI <- 
  
  function(gff, pb_ids) {
    
    gff <- gff[gff$transcript_id %in% pb_ids]
    gff_exons <- gff[gff$type == "exon"]
    gff_cds <- gff[gff$type == "CDS"]
    
    # return as list as we need both exons and cds
    gff_exons_cds <- 
      list(
        exons = data.frame(gff_exons),
        cds = data.frame(gff_cds)
      )
    
    return(gff_exons_cds)
  }


get_MANE <- 
  
  function(ref, MANE_id) {
    
    # remove any NA transcript ids (i.e. type == "gene")
    MANE <- ref[!is.na(ref$transcript_id)]
    
    # remove to .XX after ENST
    GenomicRanges::mcols(MANE)[["transcript_id"]] <- 
      GenomicRanges::mcols(MANE)[["transcript_id"]] %>%
      stringr::str_remove("\\..*")
    
    MANE <- MANE[MANE$transcript_id == mane_id, ]
    MANE_exons <- MANE[MANE$type == "exon"]
    MANE_cds <- MANE[MANE$type == "CDS"]
    
    MANE_exons_cds <- 
      list(
        exons = MANE_exons,
        cds = MANE_cds
      )
    
    return(MANE_exons_cds)
  }

annotate_TX <- 
  
  function(tx, gene) {
    
    sapply(tx, function(x){
      tmp <- Filtered_tx[[gene]]
      x$transcript_biotype <- tmp[match(x$transcript_id, tmp$isoform),]$Isoform_class
      
      return(x)
    }, simplify=F)
    
  }

plot <-
  function(gene_name, filter) {
    
    TOI <- 
      get_TOI(
        gff = LR_gff$GBA,
        pb_ids = Filtered_tx$GBA$isoform) %>% 
      annotate_TX(tx = ., gene = gene_name) 
    
    if(filter == TRUE){
      
      TOI <-
        TOI %>% 
        lapply(., function(x){dplyr::filter(x, transcript_biotype == "Coding novel")})
      
    }
    
    exons <- TOI$exons
    introns <- exons %>% to_intron(group_var = "transcript_id")
    cds <- TOI$cds
    
    final_plot <-
      exons %>%
      ggplot(aes(
        xstart = start,
        xend = end,
        y = transcript_id
      )) +
      geom_range(height = 0.25) +
      geom_range(
        data = cds
      ) +
      geom_intron(
        data = introns,
        arrow.min.intron.length = 3500,
        strand = "-"
      ) +
      labs(
        y = "",
        x = "Genomic position"
      ) +
      ggtitle(gene_name)+ 
      scale_fill_brewer(palette = "Dark2") +
      theme_bw() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    
    return(final_plot)
  }

# Main --------------------------------------------------------------------

# All filtered GBA transcripts #

GBA_filtered <- 
  get_TOI(
    gff = LR_gff$GBA,
    pb_ids = Filtered_tx$GBA$isoform) %>% 
  annotate_TX(tx = ., gene = "GBA") %>% 
  lapply(., function(x){dplyr::filter(x, transcript_biotype == "Coding novel")})


GBA_exons <- GBA_filtered$exons
GBA_introns <- GBA_exons %>% to_intron(group_var = "transcript_id")
GBA_cds <- GBA_filtered$cds

GBA_base_plot <-
  GBA_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range(height = 0.25) +
  geom_range(
    data = GBA_cds
  ) +
  geom_intron(
    data = GBA_introns,
    arrow.min.intron.length = 3500,
    strand = "-"
  ) +
  labs(
    y = "",
    x = "Genomic position"
  ) +
  ggtitle("GBA")+ 
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# All filtered GBAP1 transcripts #

GBAP1_filtered <- 
  get_TOI(
    gff = LR_gff$GBAP1,
    pb_ids = Filtered_tx$GBAP1$isoform) %>% 
  annotate_TX(tx = ., gene = "GBAP1") %>% 
  lapply(., function(x){dplyr::filter(x, transcript_biotype == "Coding novel")})


GBAP1_exons <- GBAP1_filtered$exons
GBAP1_introns <- GBAP1_exons %>% to_intron(group_var = "transcript_id")
GBAP1_cds <- GBAP1_filtered$cds

GBAP1_base_plot <-
  GBAP1_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range(height = 0.25) +
  geom_range(
    data = GBAP1_cds
  ) +
  geom_intron(
    data = GBAP1_introns,
    arrow.min.intron.length = 3500,
    strand = "-"
  ) +
  labs(
    y = "",
    x = "Genomic position"
  ) +
  ggtitle("GBAP1")+
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

GBA_GBAP1_base <- ggpubr::ggarrange(GBA_base_plot, GBAP1_base_plot)

# Save data ---------------------------------------------------------------

ggsave(
  plot = GBA_GBAP1_base, 
  filename = "GBA_GBAP1_base.png", 
  path = here::here("results", "transcripts"), 
  width = 12, 
  height = 5, 
  dpi = 600, 
  bg = "white"
)

