################################################################################
# Status: ongoing - NOT FINISHED YET
# ------------------------------------------------------------------------------
# This script contains the code to plot the proportion of 
# novel events per each mis-spliced intron stored within the IDB.
# The data is displayed by tissue and gene category (i.e. parent_gene 
# or protein-coding).
################################################################################



library(tidyverse)
library(DBI)



##################################
## LOAD lncRNA DATA AND TIDY
##################################

## Connect to the IDB
IDB_path <- "/home/sruiz/PROJECTS/splicing-project-app/intron_db/dependencies/splicing.sqlite"
con <- dbConnect(RSQLite::SQLite(), IDB_path)
dbListTables(con)


## Get all clusters used (i.e. tissues)
query <- paste0("SELECT * from 'master'")
metadata <- DBI::dbGetQuery(con, query)
all_tissues_IDB <- metadata %>% distinct(cluster) %>% pull()


## Load gene list
pseudo_parent_prot_lncRNA <- aws.s3::s3read_using(FUN = read.csv,
                                                  object = "pseudo_parent_prot_lncRNA.csv", 
                                                  bucket = "gba-manuscript") 
## Get parent genes
raw_parent_genes <- pseudo_parent_prot_lncRNA %>% filter(gene_category == "Parent_gene")
parent_genes <- QC_pseudo_parent_lncRNA_genes(raw_parent_genes)


## Get not parent genes
raw_not_parent_genes <- pseudo_parent_prot_lncRNA %>% filter(gene_category == "protein_coding")
not_parent_genes <- QC_pseudo_parent_lncRNA_genes(raw_not_parent_genes)

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

###########################################################
## QUERY THE IDB
## Get novel events from parent/not_parent genes
###########################################################

novel_enrichment_genes <- map_df(c("Parent_gene", "protein_coding"), function(category) {
  
  # category <- c("Parent_gene", "protein_coding")[1]
  gene_list <- parent_not_parent_genes %>%
    dplyr::filter(gene_category == category) %>%
    distinct(gene_name) %>%
    pull()
  
  map_df(all_tissues_IDB, function(tissue_IDB) {
    
    # tissue_IDB <- all_tissues_IDB[1]
    samples_tissue <- metadata %>% as_tibble() %>% filter(cluster == tissue_IDB) %>% nrow()
    db_IDB <- metadata %>% filter(cluster == tissue_IDB) %>% distinct(SRA_project) %>% pull()
    
    
    query <- paste0("SELECT distinct(intron.ref_coordinates), gene.gene_name,
                  intron.ref_ss5score, intron.ref_ss3score, intron.clinvar,
                  tissue.ref_type, tissue.ref_n_individuals, tissue.ref_mean_counts, tissue.MSR_D, tissue.MSR_A,
                  novel.novel_ss5score, novel.novel_ss3score, tissue.novel_n_individuals, tissue.novel_mean_counts
                  FROM '", tissue_IDB, "_", db_IDB, "_misspliced' AS tissue
                  INNER JOIN 'intron' ON intron.ref_junID=tissue.ref_junID
                  INNER JOIN 'novel' ON novel.novel_junID=tissue.novel_junID
                  INNER JOIN 'gene' ON gene.id=intron.gene_id
                  WHERE gene.gene_name IN ('", paste(gene_list,collapse = "','"), "')")
    
    
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

# Save data 
#saveRDS(object = novel_enrichment_genes, file = file_name) 



###########################################################
######    FILTER BY N% OF SAMPLES        ##################
###########################################################

## Once the IDB has been queried, the results are stored 
## within the 'novel_enrichment_genes' object.

for (p in c(2.5, 5, 7.5, 10, 20, 40, 60, 80, 90, 95)) {
  
  # p <- 2.5
  print(paste0(Sys.time(), ": ", p, "%..."))
  
  parent_genes_processed_tidy <- novel_enrichment_genes %>%
    filter(novel_n_individuals >= ((p * samples_tissue) / 100),
           gene_category == "Parent_gene") %>%
    #filter(tissue == gtex_tissues[18]) %>%
    group_by(gene_name, tissue) %>%
    dplyr::count() %>%
    dplyr::rename("num_novel_juncs" = n) %>%
    group_by(tissue) %>%
    #mutate(total_novel = sum(num_novel_juncs)) %>%
    mutate(total_genes = n()) %>% 
    distinct(tissue, .keep_all = T) %>%
    select(tissue, total_genes) %>%
    mutate(proportion = total_genes/(raw_parent_genes %>% nrow()),
           gene_category = "parent_gene")
  

  
  
  ## Second, we get the data from the not-parent genes
  
  not_parent_genes_processed_tidy <- novel_enrichment_genes %>%
    filter(novel_n_individuals >= ((p * samples_tissue) / 100),
           gene_category == "protein_coding") %>%
    group_by(gene_name, tissue) %>%
    dplyr::count() %>%
    dplyr::rename("num_novel_juncs" = n) %>%
    group_by(tissue) %>%
    #mutate(total_novel = sum(num_novel_juncs)) %>%
    mutate(total_genes = n()) %>% 
    distinct(tissue, .keep_all = T) %>%
    select(tissue, total_genes) %>%
    mutate(proportion = total_genes/(raw_not_parent_genes %>% nrow()),
           gene_category = "protein_coding") 
  

  
  
  ## Finally, we join both datasets and test them
  joined_genes <- rbind(x = parent_genes_processed_tidy,
                        y = not_parent_genes_processed_tidy) 

}

###########################################################
##############      PLOT DATA         #####################
###########################################################


## PROPORTION NOVEL MIS-SPLICED --------------------------------


joined_genes <- joined_genes %>% 
  ungroup() %>%
  mutate(gene_category = factor(gene_category, levels =c("parent_gene","protein_coding"))) %>%
  arrange(gene_category , desc(proportion)) %>%
  mutate(tissue = fct_inorder(tissue))
colours <- ifelse(str_detect(string = as.factor(joined_genes$tissue), pattern = "Brain"), "red", "black")

ggplot(data = joined_genes,
       aes(x = tissue, group = gene_category, fill = gene_category)) +
  geom_histogram(aes(y = proportion), position = position_dodge(width = 0.7), stat="identity") +
  ggtitle("RNA-seq short-read data from GTEx v8") + 
  ylab("proportion of novel junctions\nacross mis-spliced genes") +
  xlab(NULL) +
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = "12"),
        axis.title = element_text(colour = "black", size = "14"),
        plot.title = element_text(colour = "black", size = "14"),
        strip.text = element_text(colour = "black", size = "14"),
        axis.text.x = element_text(color = colours,
                                   angle = 70, 
                                   vjust = 0.99,
                                   hjust = 0.99),
        legend.text = element_text(size = "14"),
        legend.title = element_text(size = "14"),
        legend.position = "top") +
  guides(fill = guide_legend(title = "Gene type:", ncol = 2, nrow = 1)) 



## TOTAL GENES MIS-SPLICED --------------------------------

joined_genes <- joined_genes %>% 
  ungroup() %>%
  #mutate(gene_category = factor(gene_category, levels = c("protein_coding", "parent"))) %>%
  mutate(gene_category = factor(gene_category, levels = c("parent_gene","protein_coding"))) %>%
  arrange(gene_category , desc(total_genes)) %>%
  mutate(tissue = fct_inorder(tissue))
colours <- ifelse(str_detect(string = as.factor(joined_genes$tissue), pattern = "Brain"), "red", "black")

ggplot(data = joined_genes,
       aes(x = tissue, group = gene_category, fill = gene_category)) +
  geom_histogram(aes(y = total_genes), 
                 position = position_dodge(width = 0.7), 
                 stat="identity") +
  ylab("total genes mis-spliced") +
  xlab(NULL) +
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = "12"),
        axis.title = element_text(colour = "black", size = "14"),
        plot.title = element_text(colour = "black", size = "14"),
        strip.text = element_text(colour = "black", size = "14"),
        axis.text.x = element_text(color = colours,
                                   angle = 70, 
                                   vjust = 0.99,
                                   hjust = 0.99),
        legend.text = element_text(size = "14"),
        legend.title = element_text(size = "14"),
        legend.position = "top") +
  guides(fill = guide_legend(title = "Gene type:", ncol = 2, nrow = 1)) 

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


#####################################
## HELPER FUNCTIONS
#####################################

QC_pseudo_parent_lncRNA_genes <- function(gene_list) {
  
  # Three parent genes present duplicated rows
  gene_list %>% distinct(gene_name) %>% nrow()
  gene_list %>% nrow()
  gene_list %>%
    group_by(gene_name) %>% 
    filter(n()>1) %>%
    as.data.frame()
  
  # Remove duplicated rows
  gene_list <- gene_list %>% distinct(gene_name, .keep_all = T)
  
  
  return(gene_list)
}


