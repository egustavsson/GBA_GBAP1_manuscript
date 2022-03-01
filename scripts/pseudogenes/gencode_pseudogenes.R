library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(plyranges)

annotation_list <- 
  list.files("/home/egust/Projects/pseudogenes/raw_data/GENCODE_annotations", full.names = TRUE)


GENCODE_pseudogenes <- lapply(annotation_list, function(x) {
  
  print(paste0("Processing file: ", x))
  
  rtracklayer::import(x) %>% 
    plyranges::filter(str_detect(gene_type, "pseudogene") & type == "gene") %>% 
    plyranges::mutate(gencode_version = as.numeric(gsub('/home/egust/Projects/pseudogenes/raw_data/GENCODE_annotations/gencode.v|*.annotation.gtf.gz', '', x)),
                      pseudogene_type = case_when(gene_type %in% c("IG_V_pseudogene", 
                                                                   "IG_C_pseudogene", 
                                                                   "IG_J_pseudogene", 
                                                                   "IG_pseudogene") ~ "IG",
                                                  gene_type %in% c("TR_J_pseudogene", 
                                                                   "TR_V_pseudogene") ~ "TR",
                                                  gene_type %in% c("transcribed_unprocessed_pseudogene", 
                                                                   "unprocessed_pseudogene", 
                                                                   "translated_unprocessed_pseudogene") ~ "Unprocessed",
                                                  gene_type %in% c("processed_pseudogene", 
                                                                   "transcribed_processed_pseudogene", 
                                                                   "translated_processed_pseudogene") ~ "Processed",
                                                  gene_type %in% c("transcribed_unitary_pseudogene", 
                                                                   "unitary_pseudogene") ~ "Unitary",
                                                  gene_type %in% c("polymorphic_pseudogene",
                                                                   "rRNA_pseudogene",
                                                                   "pseudogene") ~ "Other"))
  
}

) %>% GRangesList(compress = F)

names(GENCODE_pseudogenes) <- sapply(GENCODE_pseudogenes, function(x){
  
  paste0("version", unique(as.numeric(x$gencode_version)))
  
})

# Save data -------------------------------------------------------------------------------------------
save(GENCODE_pseudogenes, file = "/home/egust/Projects/pseudogenes/results/GENCODE_pseudogenes.rda")
