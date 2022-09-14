library(tidyverse)
library(GenomicRanges)
library(plyranges)
library(rtracklayer)
library(here)

# Arguments ---------------------------------------------------------------


# Get GTEx median TPM
GTEx_path <- here::here(tempdir(), "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz")

if(!file.exists(GTEx_path)) {
  
  download.file(
    url = paste0(
      "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/", 
      "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
    ),
    destfile = GTEx_path
  )
  
}

# Args
args <-
  list(
    path_to_ref = here::here("data", "gencode.v38.annotation.gtf.gz"),
    path_to_OMIM = here::here("results", "pseudogenes", "omim_morbid.rda"),
    path_to_GTEx = GTEx_path)

# Load data ---------------------------------------------------------------

GENCODE_v38_pseudo <-
  rtracklayer::import(args$path_to_ref) %>% 
  plyranges::filter(str_detect(gene_type, "pseudogene"),
                    type == "gene",
                    source == "HAVANA",
                    gene_type != "polymorphic_pseudogene") %>% 
  plyranges::mutate(pseudogene_type = case_when(gene_type %in% c("transcribed_unprocessed_pseudogene", 
                                                                 "unprocessed_pseudogene", 
                                                                 "translated_unprocessed_pseudogene") ~ "Unprocessed",
                                                gene_type %in% c("processed_pseudogene", 
                                                                 "transcribed_processed_pseudogene", 
                                                                 "translated_processed_pseudogene") ~ "Processed",
                                                gene_type %in% c("IG_V_pseudogene", 
                                                                 "IG_C_pseudogene", 
                                                                 "IG_J_pseudogene", 
                                                                 "IG_pseudogene",
                                                                 "TR_J_pseudogene", 
                                                                 "TR_V_pseudogene",
                                                                 "transcribed_unitary_pseudogene", 
                                                                 "unitary_pseudogene",
                                                                 "polymorphic_pseudogene",
                                                                 "rRNA_pseudogene",
                                                                 "pseudogene") ~ "Other"))



# OMIM morbid data
load(args$path_to_OMIM)

# GTEx median TPM data
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
save(pseudogene_annotated, file = here::here("results", "pseudogenes", "pseudogene_annotated.rda"))
