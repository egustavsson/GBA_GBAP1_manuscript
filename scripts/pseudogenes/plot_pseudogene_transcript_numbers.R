# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)

# Arguments ---------------------------------------------------------------

args <-
  list(
    ENCODE_files = list.files(here::here("data", "ENCODE_LR"), full.names = TRUE, recursive = TRUE),
    path_to_pseudo = here::here("results", "pseudogenes", "pseudogene_annotated.rda"),
    path_to_GencodeV38 = here::here("data", "gencode.v38.annotation.gtf.gz")
  )

# Load data ---------------------------------------------------------------

# ENCODE long-read of frontal cortex (BA46)
expression_data = lapply(args$ENCODE_files, function(x) {
  
  read_tsv(x, show_col_types = F) %>%
    dplyr::mutate(sample = gsub(".*/", "",
                                gsub("\\..*","", x)),
                  Tissue = basename(dirname(x)),
                  annot_gene_id = str_extract(annot_gene_id, "[^.]+")) %>% # remove everything after period in ENSG IDs
    dplyr::filter(substr(annot_gene_id, 1, 4) == "ENSG") %>% # Only include known genes
    rename_at(vars(starts_with("rep")), ~paste0("Reads"))
  }
)

# Parent genes
load(args$path_to_pseudo)

pseudo_genes <-
  pseudogene_annotated %>% 
  dplyr::filter(pseudogene_type %in% c("Unprocessed", "Processed")) %>% 
  dplyr::mutate(gene_id = str_extract(gene_id, "[^.]+")) %>% 
  dplyr::select(gene_id, gene_type)

# Functions ---------------------------------------------------------------

TranscriptsExpressedGene <- 
  function(genes, expressionData) {
    
    list_out <- lapply(expressionData, function(x) {
      
      # group by gene and count number of transcripts (rows) per gene
      number_of_transcripts <-
        x %>%
        dplyr::filter(annot_gene_id %in% genes) %>% 
        group_by(annot_gene_name, annot_gene_id, sample, Tissue) %>% 
        summarise(transcripts = n(), .groups = 'keep')
    
    return(number_of_transcripts)
    
    })
  
    return(bind_rows(!!!list_out))
  
  }
    

# Main --------------------------------------------------------------------

pseudogene_transcripts <-
  TranscriptsExpressedGene(genes = pseudo_genes$gene_id, 
                           expressionData = expression_data) %>% 
  # annotate what type of pseudogene
  left_join(., pseudo_genes, by = c("annot_gene_id" = "gene_id")) %>%
  plyranges::mutate(gene_type = case_when(gene_type %in% c("transcribed_unprocessed_pseudogene", 
                                                               "unprocessed_pseudogene", 
                                                               "translated_unprocessed_pseudogene") ~ "Unprocessed",
                                              gene_type %in% c("processed_pseudogene", 
                                                               "transcribed_processed_pseudogene", 
                                                               "translated_processed_pseudogene") ~ "Processed")) %>% 
  aggregate(data = ., transcripts ~ annot_gene_id + gene_type + Tissue, FUN = mean) %>%  # use aggregate to get the mean number of transcripts
  dplyr::mutate(alternative_splicing = case_when(transcripts == 1 ~ "No (single transcript)",
                                                   transcripts >= 1 ~ "Yes (multiple transscripts)"))
number_of_transcripts_pseudogene_plot <-
  pseudogene_transcripts %>% 
  count(alternative_splicing, gene_type, Tissue) %>% 
  group_by(gene_type, Tissue) %>% 
  dplyr::mutate(Percentage = (n/sum(n) *100)) %>% 
  
  
  ggplot(aes(x = gene_type, y = Percentage, fill = alternative_splicing)) +
  geom_bar(stat = "identity", colour="black", width = 0.5) +
  facet_grid(Tissue ~ gene_type, 
             scales="free",
             space = "free") +
  labs(x = "",
       y = "Percentage of pseudogenes alternatively spliced") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_fill_manual("Alternatively spliced pseudogenes",
                    values = c("white", "lightblue"),
                    breaks = c("No (single transcript)", "Yes (multiple transscripts)"),
                    guide = guide_legend(direction = "horizontal",
                                         title.position = "top")) +
  
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.spacing = unit(0,'lines'), 
        legend.position="top",
        legend.title = element_text(size = 16,
                                    face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "bold",
                                   size = 10),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(face = "bold",
                                    size = 12),
        strip.background = element_rect(size = 1, color = "black", fill="gray83"))


# Save data ---------------------------------------------------------------

ggsave(plot = number_of_transcripts_pseudogene_plot, 
       filename = "number_of_transcripts_pseudogene_plot.svg", 
       path = here::here("results", "ENCODE_LR"), 
       width = 6, 
       height = 10, 
       dpi = 600
)
