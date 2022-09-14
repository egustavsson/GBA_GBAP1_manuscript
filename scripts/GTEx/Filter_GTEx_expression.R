library(tidyverse)
library(data.table)
library(R.utils)
library(here)

gtex <- fread("/data/GTEx_expression/GTEx_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
samples <- read_delim(file = "/data/GTEx_expression/GTEx_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", delim = "\t")
pheno <- read_delim(file = "/data/GTEx_expression/GTEx_v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", delim = "\t")

# Data filtered for GBA and GBAP1
tissues_to_exclude <- c("Testis", "Brain - Cerebellum", "Brain - Cortex", "Cells - EBV-transformed lymphocytes", "Cells - Cultured fibroblasts",
                        "Bladder", "Cervix - Ectocervix", "Cervix - Endocervix", "Fallopian Tube", "Kidney - Medulla", "Kidney - Cortex")

gtex_filtered <- gtex[gtex$Description %in% c("GBA", "GBAP1"), ]

gtex_filtered <-
  gtex_filtered %>% 
  tidyr::gather(key = "SAMPID", value = "tpm", -Name, -Description) %>% 
  tidyr::separate(SAMPID, into = c("GTEX", "subject"), sep = "-", remove = F) %>% 
  dplyr::mutate(SUBJID = str_c(GTEX, subject, sep = "-")) %>% 
  dplyr::select(-Name) %>% 
  dplyr::inner_join(samples) %>% 
  dplyr::inner_join(pheno) %>% 
  dplyr::filter(!SMTSD %in% tissues_to_exclude) %>% 
  dplyr::mutate(Organ = case_when(str_detect(SMTSD, "Brain") ~ SMTSD,
                                  TRUE ~ SMTS),
                `Brain tissue?` = case_when(str_detect(Organ, "Brain") ~ "TRUE",
                                            TRUE ~ "FALSE"))

saveRDS(gtex_filtered, 
        here::here("results", "GTEx", "gtexv8_GBA_GBAP1_tpm.Rds"))
