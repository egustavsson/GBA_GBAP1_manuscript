library(tidyverse)
library(readr)
library(GenomicRanges)
library(plyranges)
library(diffloop)
library(rtracklayer)

# Arguments ---------------------------------------------------------------

args <-
  list(
    path_to_GENCODEv10 = here::here("raw_data/GENCODE_annotations", "gencode.v10.annotation.gtf.gz"),
    path_to_OMIM = here::here("results", "pseudogenes", "omim_morbid.rda"),
    path_to_GTEx = "/data/GTEx_expression/GTEx_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
    path_to_parent = here::here("raw_data", "gencode.v10.pgene.parents.txt"))

# Load data ---------------------------------------------------------------

## GENCODE v10
GENCODEv10 <-
  rtracklayer::import(args$path_to_GENCODEv10)

## OMIM morbid
load(args$path_to_OMIM)

## GTEx v8 median TPM
GTEx_median_tpm <- 
  read.delim(args$path_to_GTEx,
             header = TRUE,
             skip = 2,
             check.names = FALSE) %>% 
  dplyr::mutate(MeanTPM = rowMeans(.[,3:length(.)]),
                ENSG_to_match = str_extract(Name, "[^.]+")) # So left_join can match regex as ENSG IDs have different numbers after the dot

## Parent genes
parentgene <- 
  read_tsv(args$path_to_parent,
           col_names = TRUE)

# Main --------------------------------------------------------------------

# Set colnames for parent genes

colnames(parentgene) <- c("ID", 
                          "pseudo_name", 
                          "pseudo_chromosome", 
                          "pseudo_strand", 
                          "pseudo_start", 
                          "pseudo_end", 
                          "Parent_gene", 
                          "Parent_transcript", 
                          "Parent_name", 
                          "seqnames", 
                          "Start", 
                          "End", 
                          "Strand")
                                 

# Get OMIM morbid information

parentgene_annotated <- 
  parentgene %>%
  dplyr::mutate(omim_morbid = Parent_name %in% omim_morbid$hgnc_symbol,
                ENSG_to_match = str_extract(Parent_gene, "[^.]+")) %>% # So left_join can match regex as ENSG IDs have different numbers after the dot

  
  # GTEx median TPM
  left_join(., GTEx_median_tpm, by = c("ENSG_to_match" = "ENSG_to_match")) %>% 
  
  # number of tissues the pseudogene is expressed in
  dplyr::mutate(no_tissues_expressed = rowSums(.[,18:71] != 0),
                is_expressed = case_when(MeanTPM != 0 ~ "Yes",
                                         MeanTPM == 0 ~ "No")) 

# Get associated pseudogene annotation

ParentAssociatedPseudo <-
  GENCODEv10 %>% 
  data.frame() %>%
  inner_join(., parentgene, by = c("transcript_id" = "ID")) %>% 
   
  distinct(transcript_id, .keep_all = TRUE) %>% 
  dplyr::mutate(ENSG_to_match = str_extract(Parent_gene, "[^.]+")) %>% 
  dplyr::select(ENSG_to_match,
                associated_pseudo = transcript_type) %>% 
  
  dplyr::mutate(associated_pseudo = case_when(associated_pseudo %in% c("IG_V_pseudogene", 
                                                                       "IG_C_pseudogene", 
                                                                       "IG_J_pseudogene", 
                                                                       "IG_pseudogene") ~ "IG",
                                              associated_pseudo %in% c("TR_J_pseudogene", 
                                                                       "TR_V_pseudogene") ~ "TR",
                                              associated_pseudo %in% c("transcribed_unprocessed_pseudogene", 
                                                                       "unprocessed_pseudogene", 
                                                                       "translated_unprocessed_pseudogene") ~ "Unprocessed",
                                              associated_pseudo %in% c("processed_pseudogene", 
                                                                       "transcribed_processed_pseudogene", 
                                                                       "translated_processed_pseudogene") ~ "Processed",
                                              associated_pseudo %in% c("transcribed_unitary_pseudogene", 
                                                                       "unitary_pseudogene") ~ "Unitary",
                                              associated_pseudo %in% c("polymorphic_pseudogene",
                                                                       "rRNA_pseudogene",
                                                                       "pseudogene") ~ "Other"))

parentgene_annotated <- 
  parentgene_annotated %>% 
  dplyr::left_join(., ParentAssociatedPseudo, by = c("ENSG_to_match" = "ENSG_to_match")) %>% 
  distinct(., Parent_name, .keep_all = TRUE)

# Save data -------------------------------------------------------------------------------------------
save(parentgene_annotated, file = here::here("results", "parentgene_annotated.rda"))
