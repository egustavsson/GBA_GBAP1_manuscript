library(biomaRt)
library(tidyverse)
library(here)
library(diffloop)

# Load data ---------------------------------------------------------------

# take the mim2gene list 
mim2gene <- readr::read_delim("https://omim.org/static/omim/data/mim2gene.txt", delim = "\t", skip = 4)

# Main --------------------------------------------------------------------

ens_mart <- biomaRt::useEnsembl(biomart = "ensembl", 
                                dataset = "hsapiens_gene_ensembl", 
                                version = 102)

# get all ensembl ids
omim_gene_ens <- mim2gene[["Ensembl Gene ID (Ensembl)"]][!is.na(mim2gene[["Ensembl Gene ID (Ensembl)"]])]

# look for all ensembl for an associated mim morbid phenotype
omim_morbid <- biomaRt::getBM(mart = ens_mart,
                              attributes = c("ensembl_gene_id", 
                                             "hgnc_symbol",
                                             "chromosome_name",
                                             "start_position",
                                             "end_position",
                                             "strand",
                                             "mim_morbid_accession", 
                                             "mim_morbid_description"),
                              filters = "ensembl_gene_id", 
                              values = omim_gene_ens, 
                              uniqueRows = TRUE)

# then filter for only those with a non-NA morbid phenotype
omim_morbid <- omim_morbid %>% 
  as_tibble() %>% 
  dplyr::filter(!is.na(mim_morbid_accession))

# Omim morbid with correct strand information
omim_morbid$strand <- gsub("^1$", "+", gsub("^-1$", "-", omim_morbid$strand))

# Save data ---------------------------------------------------------------
save(omim_morbid, file = here::here("results", "pseudogenes", "omim_morbid.rda"))
