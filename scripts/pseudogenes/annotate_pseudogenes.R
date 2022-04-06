library(tidyverse)
library(GenomicRanges)
library(plyranges)
library(rtracklayer)
library(here)

# Arguments ---------------------------------------------------------------

args <-
  list(
    path_to_GENCODE = here::here("results", "GENCODE_pseudogenes.rda"),
    path_to_OMIM = here::here("results", "omim_morbid.rda"),
    path_to_GTEx = "/data/GTEx_expression/GTEx_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")

# Load data ---------------------------------------------------------------
load(args$path_to_GENCODE)

GENCODE_v38_pseudo <-
  GENCODE_pseudogenes$version38 %>% 
  plyranges::filter(source == "HAVANA",
                    gene_type != "polymorphic_pseudogene")

load(args$path_to_OMIM)

GTEx_median_tpm <- 
  read.delim(args$path_to_GTEx,
             header = TRUE,
             skip = 2,
             check.names = FALSE) %>% 
  dplyr::mutate(MeanTPM = rowMeans(.[,3:length(.)]),
                ENSG_to_match = str_extract(Name, "[^.]+")) # So left_join can match regex as ENSG IDs have different numbers after the dot

# Main --------------------------------------------------------------------

pseudogene_annotated <-
  
  GENCODE_v38_pseudo %>%
  data.frame() %>% 
  
  # Get OMIM morbid information
  dplyr::mutate(omim_morbid = gene_name %in% omim_morbid$hgnc_symbol,
                ENSG_to_match = str_extract(gene_id, "[^.]+")) %>% # So left_join can match regex as ENSG IDs have different numbers after the dot
  
  # GTEx median TPM
  left_join(., GTEx_median_tpm, by = c("ENSG_to_match" = "ENSG_to_match")) %>% 
  
  # number of tissues the pseudogene is expressed in
  dplyr::mutate(no_tissues_expressed = rowSums(.[,33:86] != 0),
                is_expressed = case_when(MeanTPM != 0 ~ "Yes",
                                         MeanTPM == 0 ~ "No"))
  
# Save data -------------------------------------------------------------------------------------------
save(pseudogene_annotated, file = here::here("results", "pseudogene_annotated.rda"))
