# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(GenomicRanges)
library(plyranges)

# Load data ---------------------------------------------------------------

# PacBio filtered transcripts from Brain (SQANTI output)
Classification <-
  list(
    GBA = read_tsv(here::here("data", "Brain", "GBA_classification_filtered.txt"), show_col_types = F),
    GBAP1 = read_tsv(here::here("data", "Brain", "GBAP1_classification_filtered.txt"), show_col_types = F)
)  

# GTF
Structures <-
  GRangesList(
    GBA = rtracklayer::import(here::here("data", "Brain", "GBA_corrected.gtf.cds.gff")),
    GBAP1 = rtracklayer::import(here::here("data", "Brain", "GBAP1_corrected.gtf.cds.gff"))
)

# GBA and GBAP1 loci to use when filtering data
locus_subset <- 
  GRanges(seqnames = "chr1",
          ranges = IRanges(start = 155213078, 
                           end = 155245198),
          strand = "-")

# FANTOM5 CAGE seq data
CAGE_seq <- 
  rtracklayer::import.bed(here::here("data", "hg38.cage_peak_phase1and2combined_coord.bed.gz")) %>% 
  subsetByOverlaps(., locus_subset)


# Functions ---------------------------------------------------------------

dist_to_CAGE <- function(data, CAGE) {
  
  data_transcript <- 
    data %>% plyranges::filter(type == "transcript")
  
  table <- data.frame(Transcript = data_transcript$transcript_id,
                      Dist_to_CAGE = as.numeric(0))
  
  for (i in table$Transcript) {
    
    gr <-
      data_transcript %>% 
      plyranges::filter(transcript_id == i)
    
    if (as.character(strand(gr)) == "-") {
      
      start(gr) <- end(gr)
      GenomicRanges::distance(gr, CAGE)
      table$Dist_to_CAGE[table$Transcript == i] <- min(GenomicRanges::distance(gr, CAGE))
      
    }else {
      
      end(gr) <- start(gr)
      GenomicRanges::distance(gr, CAGE) 
      table$Dist_to_CAGE[table$Transcript == i] <- min(GenomicRanges::distance(gr, CAGE))
    } 
  }
  
  return(table)
}

# Main --------------------------------------------------------------------

# Transcript IDs to filter by
TOI <- bind_rows(!!!Classification) %>% dplyr::select(isoform)

Transcripts <-
  Structures %>% 
  unlist() %>% 
  plyranges::filter(type == "transcript",
                    transcript_id %in% TOI$isoform)


# Distance to nearest CAGE peak

CAGE_dist_table <-
  dist_to_CAGE(data = Transcripts, CAGE = CAGE_seq) %>% 
  left_join(., bind_rows(!!!Classification), by = c("Transcript" = "isoform")) %>% 
  dplyr::select(Transcript,
                Dist_to_CAGE,
                associated_gene,
                Isoform_class)

# Save data ---------------------------------------------------------------

write_csv(CAGE_dist_table, 
          file = here::here("results", "transcripts", "CAGE_dist_table.csv"))
