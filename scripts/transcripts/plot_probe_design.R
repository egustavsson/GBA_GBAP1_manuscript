# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggtranscript)
library(plyranges)

# Load data ---------------------------------------------------------------

probes <-
  read_csv(here::here("data", "probe_design_final.csv"))

ref <-
  rtracklayer::import(here::here("data", "gencode.v38.annotation.gtf.gz"))

# Functions ---------------------------------------------------------------



# Main --------------------------------------------------------------------

GBA_GBAP1_probes <-
  probes %>% 
  dplyr::filter(grepl("GBA", Group.Name)) %>% 
  dplyr::mutate(transcript_id = "Probes",
                gene_name = case_when(Group.Name == "GBA(2629)" ~ "GBA",
                                      Group.Name == "GBAP1(2630)" ~ "GBAP1")) %>% 
  dplyr::select(seqnames,
                start,
                end,
                strand,
                gene_name,
                transcript_id)




GBA_GBAP1_annotation <-
  ref %>% 
  plyranges::filter(gene_name %in% c("GBA", "GBAP1")) %>% 
  data.frame()

GBA_GBAP1_exons <-
  GBA_GBAP1_annotation %>% 
  dplyr::filter(type == "exon")

GBA_GBAP1_introns <- GBA_GBAP1_exons %>% 
  to_intron(group_var = "transcript_id")

GBA_GBAP1_cds <-
  GBA_GBAP1_annotation %>% 
  dplyr::filter(type == "CDS")


# probes are hg19

probe_design_plot <-
  GBA_GBAP1_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range(
    fill = "white",
    height = 0.25
  ) +
  geom_range(
    data = GBA_GBAP1_cds
  ) +
  geom_range(
    data = GBA_GBAP1_probes 
  ) +
  geom_intron(
    data = GBA_GBAP1_introns
  ) +
  labs(
    y = "Transcript name",
    x = "Genomic position"
  ) +
  facet_wrap(vars(gene_name), ncol = 1, scales = "free")



# Save data ---------------------------------------------------------------