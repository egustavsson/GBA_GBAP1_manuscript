################################################################################
# Status: finished
# How to use it: just go to line 207 to plot the proportion of genes using a threshold of 5% of individuals.
# ------------------------------------------------------------------------------
# This script contains the code to plot the proportion of parent/protein-coding genes containing
# novel events shared by >= 'threshold'% of individuals.
# The data is displayed by tissue and gene category (i.e. parent_gene 
# and protein-coding categories).
################################################################################

library(tidyverse)
library(DBI)


#####################################
## FUNCTIONS
## ---------------------------------
## This section contains the core functions that generate a data frame with the proportion of 
## parent/protein-coding genes per tissue from the IDB with novel introns shared by 'threshold'% of individuals.
#####################################

get_lncRNA_gene_list <- function() {
  
  ## Stablish connection with AWS S3
  Sys.setenv("AWS_ACCESS_KEY_ID"=Sys.getenv("AWS_ACCESS_KEY_ID"),
             "AWS_SECRET_ACCESS_KEY"=Sys.getenv("AWS_SECRET_ACCESS_KEY"),
             "AWS_DEFAULT_REGION"=Sys.getenv("AWS_DEFAULT_REGION"))
  
  
  ## Load lncRNA gene list
  pseudo_parent_prot_lncRNA <- aws.s3::s3read_using(FUN = read.csv,
                                                    object = "pseudo_parent_prot_lncRNA.csv", 
                                                    bucket = "gba-manuscript") 
  ## Get parent genes
  raw_parent_genes <- pseudo_parent_prot_lncRNA %>% filter(gene_category == "Parent_gene")
  parent_genes <- raw_parent_genes %>% distinct(gene_name, .keep_all = T)
  
  ## Get not parent genes
  raw_not_parent_genes <- pseudo_parent_prot_lncRNA %>% filter(gene_category == "protein_coding")
  not_parent_genes <- raw_not_parent_genes %>% distinct(gene_name, .keep_all = T)
  
  ## Remove overlapping genes between both lists
  common_genes <- intersect(parent_genes %>% distinct(gene_name) %>% pull(),
                            not_parent_genes %>% distinct(gene_name) %>% pull())
  parent_genes <- parent_genes %>%
    filter(!(gene_name %in% common_genes))
  not_parent_genes <- not_parent_genes %>%
    filter(!(gene_name %in% common_genes))
  
  ## Once both lists are tidy, I join them both
  ## The column (gene_category indicates whether they are parent or protein-coding genes)
  parent_not_parent_genes <- rbind(parent_genes, not_parent_genes)
  
  return(parent_not_parent_genes)
  
}


query_IDB <- function(df_lncRNA_genes, 
                      categories = c("Parent_gene", "protein_coding")) {
  
  #########################
  ## This function queries the IDB.
  ## It obtains all novel events from the genes received by parameter (i.e. 'df_lncRNA_genes' parameter)
  ## by each tissue from the IDB.
  ## This query also returns the % of individuals figure sharing each novel event.
  #########################
  
  ## Connect to the IDB
  IDB_path <- "/home/sruiz/PROJECTS/splicing-project-app/intron_db/dependencies/splicing.sqlite"
  con <- dbConnect(RSQLite::SQLite(), IDB_path)
  DBI::dbListTables(con)
  
  ## Get all clusters used (i.e. tissues)
  query <- paste0("SELECT * from 'master'")
  IDB_metadata <- DBI::dbGetQuery(con, query)
  all_tissues_IDB <- IDB_metadata %>% distinct(cluster) %>% pull()
  
  map_df(categories, function(category) {
    
    # category <- c("Parent_gene", "protein_coding")[1]
    gene_list <- df_lncRNA_genes %>%
      dplyr::filter(gene_category == category) %>%
      distinct(gene_name) %>%
      pull()
    
    map_df(all_tissues_IDB, function(tissue_IDB) {
      
      # tissue_IDB <- all_tissues_IDB[1]
      samples_tissue <- IDB_metadata %>% as_tibble() %>% filter(cluster == tissue_IDB) %>% nrow()
      db_IDB <- IDB_metadata %>% filter(cluster == tissue_IDB) %>% distinct(SRA_project) %>% pull()
      
      
      query <- paste0("SELECT distinct(intron.ref_coordinates), 
                      gene.gene_name, 
                      novel.novel_coordinates, 
                      tissue.novel_n_individuals, tissue.novel_mean_counts
                      FROM '", tissue_IDB, "_", db_IDB, "_misspliced' AS tissue
                      INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                      INNER JOIN 'novel' ON novel.novel_junID=tissue.novel_junID
                      INNER JOIN 'gene' ON gene.id=intron.gene_id
                      WHERE gene.gene_name IN ('", paste(gene_list, collapse = "','"), "')")
      
      
      #query <- "SELECT * FROM novel LIMIT 10"
      df_novel_gr <- dbGetQuery(con, query) 
      print(df_novel_gr %>% nrow())
      
      if (df_novel_gr %>% nrow() > 0) {
        
        print(paste0(category, " - ", tissue_IDB))
        
        df_novel_gr %>%
          mutate(tissue = tissue_IDB,
                 samples_tissue = samples_tissue,
                 gene_category = category) %>%
          return()
      } 
      
    })
    
  })
}

filter_IDB_data_by_threshold <- function(df_lncRNA_genes_IDB,
                                         df_lncRNA_genes,
                                         threshold = 5) {
  
  #########################
  ## This function filters the IDB object returned by the function 'query_IDB'.
  ## It obtains all novel events shared by at least % 'threshold' of individuals.
  ## It also calculates the proportion of genes containing introns mis-spliced by at least
  ## 'threshold' of individuals.
  #########################
  
  
  ## Minimum number of % of individuals that the novel junction has to be shared with
  print(paste0(Sys.time(), ": ", threshold, "%..."))
  

  ## Filter parent genes
  parent_genes_processed_tidy <- df_lncRNA_genes_IDB %>%
    filter(novel_n_individuals >= ((threshold * samples_tissue) / 100),
           gene_category == "Parent_gene") %>%
    group_by(tissue) %>%
    distinct(gene_name) %>%
    dplyr::mutate("total_genes_tissue" = n()) %>%
    distinct(tissue, .keep_all = T) %>%
    mutate(proportion = total_genes_tissue/(df_lncRNA_genes %>% 
                                              filter(gene_category == "Parent_gene") %>%
                                              distinct(gene_name) %>%
                                              nrow()),
           gene_category = "parent_gene") %>% select(tissue, proportion, gene_category)
  
  
  ## Filter non-parent genes
  not_parent_genes_processed_tidy <- df_lncRNA_genes_IDB %>%
    filter(novel_n_individuals >= ((threshold * samples_tissue) / 100),
           gene_category == "protein_coding") %>%
    group_by(tissue) %>%
    distinct(gene_name) %>%
    dplyr::mutate("total_genes_tissue" = n()) %>%
    distinct(tissue, .keep_all = T) %>%
    mutate(proportion = total_genes_tissue/(df_lncRNA_genes %>% 
                                              filter(gene_category == "protein_coding") %>%
                                              distinct(gene_name) %>%
                                              nrow()),
           gene_category = "protein_coding") %>% select(tissue, proportion, gene_category)
  
  
  ## Finally, we join both datasets and test them
  return(rbind(x = parent_genes_processed_tidy,
               y = not_parent_genes_processed_tidy))
}

###################################
## PIPELINE
## -------------------------------
## This section contains the code to call the functions required to generate a data frame with the proportion of 
## parent/protein-coding genes per tissue with novel introns shared by 'threshold'% of individuals.
###################################

threshold <- 5

## Load and tidy the lncRNA list of genes
df_lncRNA_genes <- get_lncRNA_gene_list()


## Query the IDB to obtain novel events from lncRNA genes
df_lncRNA_genes_IDB <- query_IDB(df_lncRNA_genes = df_lncRNA_genes)


## Filter by %threshold of individuals and calculate the proportion
df_lncRNA_genes_IDB_tidy <- filter_IDB_data_by_threshold(df_lncRNA_genes_IDB, 
                                                         df_lncRNA_genes, ## We need this parameter to obtain the original number of genes in the lncRNA list
                                                         threshold)


###################################
## PLOT DATA    
## --------------------------------
## This section plots the proportion of parent and protein-coding genes per 
## tissue containing novel introns shared by at least 5% of individuals.
###################################

## Load file
threshold <- 5
df_lncRNA_genes_IDB_tidy <- readRDS(file = "./scripts/intron_db/prop_genes_IDB_threshold_5.rds")

## Order the data-frame by parent_genes desc
df_lncRNA_genes_IDB_tidy <- df_lncRNA_genes_IDB_tidy %>% 
  ungroup() %>%
  mutate(gene_category = str_replace(string = gene_category, pattern = "_", replacement = " ")) %>%
  mutate(gene_category = factor(gene_category, levels = c("parent gene","protein coding"))) %>%
  arrange(gene_category , desc(proportion)) %>%
  mutate(tissue = fct_inorder(tissue))

## Brain tissues will be highlighted in red
colours <- ifelse(str_detect(string = as.factor(df_lncRNA_genes_IDB_tidy$tissue), pattern = "Brain"), "red", "black")

## Plot histogram
ggplot(data = df_lncRNA_genes_IDB_tidy,
       aes(x = tissue, group = gene_category, fill = gene_category)) +
  geom_histogram(aes(y = proportion), position = position_dodge(width = 0.7), stat="identity") +
  ggtitle(paste0("RNA-seq short-read data from GTEx v8\n",
                 "Proportion of genes with novel junctions misspliced in\nat least ", threshold, "% of individuals.")) + 
  ylab("proportion of genes") +
  xlab(NULL) +
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = "12"),
        axis.title = element_text(colour = "black", size = "13"),
        plot.title = element_text(colour = "black", size = "13"),
        strip.text = element_text(colour = "black", size = "13"),
        axis.text.x = element_text(color = colours,
                                   angle = 70, 
                                   vjust = 0.99,
                                   hjust = 0.99),
        legend.text = element_text(size = "13"),
        legend.title = element_text(size = "13"),
        legend.position = "top") +
  guides(fill = guide_legend(title = "Gene type:", ncol = 2, nrow = 1)) 



file_name <- paste0("./results/intron_db/prop_genes_threshold_5_IDB_GTExv8.png")
ggplot2::ggsave(filename = file_name,
                width = 183, height = 183, units = "mm", dpi = 300)

## TOTAL GENES MIS-SPLICED --------------------------------

# joined_genes <- joined_genes %>% 
#   ungroup() %>%
#   #mutate(gene_category = factor(gene_category, levels = c("protein_coding", "parent"))) %>%
#   mutate(gene_category = str_replace(string = gene_category, pattern = "_", replacement = " ")) %>%
#   mutate(gene_category = factor(gene_category, levels = c("parent gene","protein coding"))) %>%
#   arrange(gene_category , desc(total_genes)) %>%
#   mutate(tissue = fct_inorder(tissue))
# colours <- ifelse(str_detect(string = as.factor(joined_genes$tissue), pattern = "Brain"), "red", "black")
# 
# ggplot(data = joined_genes,
#        aes(x = tissue, group = gene_category, fill = gene_category)) +
#   geom_histogram(aes(y = total_genes), 
#                  position = position_dodge(width = 0.7), 
#                  stat="identity") +
#   ylab("total genes mis-spliced") +
#   xlab(NULL) +
#   theme(axis.line = element_line(colour = "black"), 
#         axis.text = element_text(colour = "black", size = "12"),
#         axis.title = element_text(colour = "black", size = "14"),
#         plot.title = element_text(colour = "black", size = "14"),
#         strip.text = element_text(colour = "black", size = "14"),
#         axis.text.x = element_text(color = colours,
#                                    angle = 70, 
#                                    vjust = 0.99,
#                                    hjust = 0.99),
#         legend.text = element_text(size = "14"),
#         legend.title = element_text(size = "14"),
#         legend.position = "top") +
#   guides(fill = guide_legend(title = "Gene type:", ncol = 2, nrow = 1)) 

# ###########################################################
# ############    WILCOXON TEST         #####################
# ###########################################################
# 
# wilcox.test(x = joined_genes %>% filter(gene_category == "parent_gene") %>% pull(proportion),
#             y = joined_genes %>% filter(gene_category == "protein_coding") %>% pull(proportion),
#             alternative = "greater")
# 
# 
# 
# ###########################################################
# ############    LINEAR MODELS         #####################
# ###########################################################
# 
# parent_genes_raw <- readRDS(file = paste0(folder_root, "/raw_parent_genes.rds"))
# parent_genes_processed <- readRDS(file = paste0(folder_root, "/novel_enrichment_parent_genes_40GTExTissues.rds"))
# 
# 
# not_parent_genes_raw <- readRDS(file = paste0(folder_root, "/raw_not_parent_genes.rds"))
# not_parent_genes_processed <- readRDS(file = paste0(folder_root, "/novel_enrichment_not_parent_genes_40GTExTissues.rds"))
# 
# 
# for (tissue in gtex_tissues[11]) {
#   
#   # tissue <- gtex_tissues[11]
#   df_introns <- readRDS(file = paste0("/home/sruiz/PROJECTS/splicing-project/results/pipeline3/missplicing-ratio/", 
#                                       tissue, "/v104/", tissue, "_db_introns.rds")) %>% as_tibble()
#   
#   
#   ## PARENT GENES ---------------------------------------------------------------------------------------------------------------------
#   parent_genes_misspliced <- parent_genes_processed %>%
#     filter(tissue == tissue) %>%
#     distinct(gene_name) %>%
#     mutate(gene_type = "parent",
#            misspliced = T)
#   
#   parent_genes_not_misspliced <- data.frame(gene_name = (setdiff(parent_genes_raw$gene_name,
#                                                                  parent_genes_misspliced$gene_name)),
#                                             gene_type = "parent",
#                                             misspliced = F)
#   
#   all_parent_genes = rbind(parent_genes_processed_grouped,
#                            parent_genes_not_misspliced)
#   
#   
#   ## NOT PARENT GENES ------------------------------------------------------------------------------------------------------------------
#   
#   not_parent_genes_processed_grouped <- not_parent_genes_processed %>%
#     filter(tissue == tissue) %>%
#     group_by(gene_name) %>%
#     dplyr::count() %>%
#     dplyr::rename("num_novel_juncs" = n) %>%
#     mutate(gene_type = "not_parent",
#            num_novel_juncs = 1)
#   
#   not_parent_genes_not_misspliced <- data.frame(gene_name = (setdiff(not_parent_genes_raw$gene_name,
#                                                                      not_parent_genes_processed_grouped$gene_name)),
#                                                 gene_type = "not_parent",
#                                                 num_novel_juncs = 0)
#   
#   all_not_parent_genes = rbind(not_parent_genes_processed_grouped,
#                                not_parent_genes_not_misspliced)
#   
#   
#   ## JOIN PARENT AND NOT PARENT GENES -----------------------------------------------------------------------------------------------------
#   
#   all_genes <- rbind(all_parent_genes, all_not_parent_genes)
#   
#   parent_genes_processed_grouped %>% nrow()
#   parent_genes_not_misspliced %>% nrow()
#   
#   not_parent_genes_processed_grouped %>% nrow()
#   not_parent_genes_not_misspliced %>% nrow()
#   
#   tab <- matrix(c(parent_genes_processed_grouped %>% nrow(),
#                   
#                   not_parent_genes_processed_grouped %>% nrow(),
#                   parent_genes_not_misspliced %>% nrow(),
#                   not_parent_genes_not_misspliced %>% nrow()), ncol=2, byrow=TRUE)
#   
#   colnames(tab) <- c('parent','not-parent')
#   rownames(tab) <- c('misspliced', "not-misspliced")
#   
#   chisq.test(tab)
#   
#   ## ADD GENE FEATURES -----------------------------------------------------------------------------------------------------------------
#   
#   df_introns_tidy <- df_introns %>% 
#     filter(gene_name %in% (all_genes$gene_name %>% unique())) %>%
#     distinct(gene_name, .keep_all = T) %>%
#     select(num_transcripts = transcript_number, 
#            gene_width,
#            tpm,
#            gene_name) %>%
#     unnest(gene_name)
#   
#   df_merged <- merge(x = all_genes,
#                      y = df_introns_tidy,
#                      by="gene_name")
#   
#   df_merged %>% head()
#   
#   ## LINEAR MODELS -----------------------------------------------------------------------------------------------------------------
#   
#   lm(num_novel_juncs ~ gene_type + 
#        num_transcripts + 
#        gene_width + 
#        tpm, 
#      data = df_merged) %>% summary() %>% print()
# }




