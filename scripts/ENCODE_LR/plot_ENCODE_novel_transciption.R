# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)

# Arguments ---------------------------------------------------------------

args <-
  list(
    ENCODE_files = list.files(here::here("data", "ENCODE_LR"), full.names = TRUE, recursive = TRUE),
    path_to_parent = here::here("results", "pseudogenes", "parentgene_annotated.rda"),
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
load(args$path_to_parent)

parent_genes <-
  parentgene_annotated %>% 
  dplyr::filter(associated_pseudo %in% c("Unprocessed", "Processed")) %>% 
  dplyr::mutate(gene_id = str_extract(Parent_gene, "[^.]+"),
                gene_type = "parent_gene") %>% 
  dplyr::select(gene_id, gene_type)

# Protein coding genes
prot_coding_genes <-
  rtracklayer::import(args$path_to_GencodeV38) %>% 
  data.frame() %>% 
  dplyr::mutate(gene_id = str_extract(gene_id, "[^.]+")) %>% 
  dplyr::filter(type == "gene",
                gene_type == "protein_coding",
                !gene_id %in% parent_genes$gene_id) %>%
  dplyr::select(gene_id, gene_type)

# Functions ---------------------------------------------------------------

Novel_transcripts_per_sample <- 
  function(list_in) {
    
    list_out <-
      lapply(list_in, function(x) {
        
        # add transcript novelty
        y <-
          x %>%
          # dplyr::mutate(transcript_type = case_when(transcript_novelty == "Known" ~ "known",
          #                                           transcript_novelty != "Known" ~ "novel"))
          
          dplyr::filter(transcript_novelty != "Genomic")
          
          
        # get number of novel and known transcripts per gene
        num_transcripts <-
          y %>%
          group_by(annot_gene_id, transcript_novelty, sample, Tissue) %>% tally()
        
          return(num_transcripts)
      })
    
    return(bind_rows(!!!list_out))
    
  }

# Main --------------------------------------------------------------------

ENCODE_novel_transcripts <-
  Novel_transcripts_per_sample(list_in = expression_data)


# fraction of transcripts per category
ENCODE_transcript_per_category <-
  ENCODE_novel_transcripts %>% 
  pivot_wider(names_from = "transcript_novelty", values_from = "n") %>% 
  ungroup() %>% 
  na.omit() %>% 
  dplyr::mutate(Known = Known / (rowSums(dplyr::select(., Known, NIC, NNC, ISM), na.rm = T)),
                NIC = NIC / (rowSums(dplyr::select(., Known, NIC, NNC, ISM), na.rm = T)),
                NNC = NNC / (rowSums(dplyr::select(., Known, NIC, NNC, ISM), na.rm = T)),
                ISM = ISM / (rowSums(dplyr::select(., Known, NIC, NNC, ISM), na.rm = T))) %>% 
  left_join(., 
            rbind(parent_genes, prot_coding_genes),
            by = c("annot_gene_id" = "gene_id")) %>%  # add gene annotation (parent or protein)
  na.omit() %>% 
  pivot_longer(!c("gene_type", "annot_gene_id", "sample", "Tissue"), names_to = "Transcript_type", values_to = "Fraction") %>% 
  dplyr::filter(Transcript_type == "NNC")

# Plot all
ENCODE_transcript_per_category_plot <- 
  ENCODE_transcript_per_category %>% 
  ggplot(aes(x = gene_type,
             y = Fraction)) +
  geom_violin(aes(fill = gene_type)) +
  geom_boxplot(width = 0.5) +
  ggpubr::stat_compare_means(label.x = 1.2,
                             label.y = 0.65,
                             size = 5) +
  scale_x_discrete(name = "",
                   labels = c("parent_gene" = "Parent genes", "protein_coding" = "Protein coding\ngenes")) + 
  scale_y_continuous(name = "Proportion of trasncripts that are novel\nwith at least one new splice site (%)",
                     labels = function(x) paste0(x * 100, "%"), 
                     limits = c(0, 0.7)) +
  scale_fill_manual(guide = "none",
                    values = c("#5ab4ac", "#ffffff")) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))

# Plot per tissue
ENCODE_transcript_per_category_per_tissue_plot <- 
  ENCODE_transcript_per_category %>% 
  ggplot(aes(x = gene_type,
             y = Fraction)) +
  geom_violin(aes(fill = gene_type)) +
  geom_boxplot(width = 0.5) +
  ggpubr::stat_compare_means(label.x = 1.2,
                             label.y = 0.65,
                             size = 5) +
  scale_x_discrete(name = "",
                   labels = c("parent_gene" = "Parent genes", "protein_coding" = "Protein coding\ngenes")) + 
  scale_y_continuous(name = "Proportion of trasncripts that are novel\nwith at least one new splice site (%)",
                     labels = function(x) paste0(x * 100, "%"), 
                     limits = c(0, 0.7)) +
  scale_fill_manual(guide = "none",
                    values = c("#5ab4ac", "#ffffff")) +
  facet_wrap(vars(Tissue)) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(face = "bold",
                                  size = 12),
        strip.background = element_rect(size = 1, color = "black", fill="gray83"))

# Save data ---------------------------------------------------------------

ggsave(plot = ENCODE_transcript_per_category_plot, 
       filename = "ENCODE_transcript_per_category_plot.svg", 
       path = here::here("results", "ENCODE_LR"), 
       width = 6, 
       height = 5, 
       dpi = 600
)

ggsave(plot = ENCODE_transcript_per_category_per_tissue_plot, 
       filename = "ENCODE_transcript_per_category_per_tissue_plot.png", 
       path = here::here("results", "ENCODE_LR"), 
       width = 12, 
       height = 5, 
       dpi = 600
)
