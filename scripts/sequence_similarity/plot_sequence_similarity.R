# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(plyranges)
library(GenomicFeatures)
library(GenomicRanges)
library(biomaRt)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(diffloop)

# Arguments ---------------------------------------------------------------

args <-
  list(
    path_to_parents = here::here("data", "gencode.v10.pgene.parents.txt"),
    path_to_GencodeV10 = here::here("data", "gencode.v10.annotation.gtf.gz"),
    path_to_GencodeV38 = here::here("data", "gencode.v38.annotation.gtf.gz"),
    path_to_ensemblV105 = here::here("data", "Homo_sapiens.GRCh38.105.gtf.gz")
  )

# Load data ---------------------------------------------------------------

# parent_pseudo_homology <- 
#   read_tsv("/home/egust/Projects/pseudogenes/raw_data/similarity.dat", col_names = TRUE)

parentgenes <- 
  read_tsv(args$path_to_parents,
           col_names = TRUE, 
           show_col_types = F)

ref <- 
  list(
    GenV10 = rtracklayer::import(args$path_to_GencodeV10),
    ensV105 = rtracklayer::import(args$path_to_ensemblV105))

# Functions ---------------------------------------------------------------


# will add commandArgs()

split_and_index <-
  function(data, query, subject) {
    
    x <- data[, query]
    x[[1]] <- 
      case_when(
        grepl(pattern = "ENSG000", x[[1]]) ~ stringr::str_remove(x[[1]], "\\..*"), 
        TRUE ~ as.character(x[[1]]))
    x$index <- seq.int(length.out = nrow(x))
    
    y <- data[, subject]
    y[[1]] <- 
      case_when(
        grepl(pattern = "ENSG000", y[[1]]) ~ stringr::str_remove(y[[1]], "\\..*"), 
        TRUE ~ as.character(y[[1]]))
    y$index <- seq.int(length.out = nrow(y))
    
    return(
      list(query = x,
           subject = y))
  }

get_ensembl_canonical <-
  function(gene_list) {
    
    # biomart - get ensembl canonical
    ensembl_105_canonical <-
      getBM(mart = useEnsembl(biomart = "genes", 
                              dataset = "hsapiens_gene_ensembl",
                              version = 105),
            filters = "transcript_is_canonical",
            values = TRUE,
            attributes = c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol"))
    
    for (i in names(gene_list)) {
      
      if (any(grepl(pattern = "ENSG000", gene_list[[i]][[1]]))) {
        
        tmp <- match(gene_list[[i]][[1]], ensembl_105_canonical$ensembl_gene_id)
        gene_list[[i]]$canonical_transcript <- ensembl_105_canonical[tmp, "ensembl_transcript_id"]
        
      }else{
        
        tmp <- match(gene_list[[i]][[1]], ensembl_105_canonical$hgnc_symbol)
        gene_list[[i]]$canonical_transcript <- ensembl_105_canonical[tmp, "ensembl_transcript_id"]
        
      }
    }
    return(gene_list)
    
  }

remove_missing_canonical <- 
  function(data) {
    
    NA_rows <- list()
    
    for (i in names(data)) {
      
      NA_rows[[i]] <- which(is.na(data[[i]]$canonical_transcript))
      
      print(
        paste0(
          "missing ", i, " - ", sum(is.na(data[[i]]$canonical_transcript))))
      
    }
    
    rows_to_remove <-
      unique(
        Reduce(c, NA_rows)
      )
    
    return(
      lapply(data, function(x) {
        
        x[-rows_to_remove, ]
        
      }))
  }

# GRangeslist_canonical <-
#   function(data, reference) {
#     
#     # Edit reference
#     
#     reference <- 
#       ifelse(
#         any(
#           grepl(pattern = "chr", as.character(seqnames(reference)))),
#         print(reference),
#         print(addchr(reference))) %>% 
#         
#         keepStandardChromosomes(pruning.mode = "coarse")
#     
#     lapply(data, function(x) {
#       
#       x <- x %>% 
#         add_column(FiveUTR = NA,
#                    CDS = NA,
#                    ThreeUTR = NA)
#       
#       for (i in 1:nrow(x)) {
#         
#         x$CDS[i] <- unlist(getSeq(Hsapiens, dplyr::filter(ref$ensV105, transcript_id == x$canonical_transcript[i], type == "CDS"))) %>% as.character()
#         
#       }
#       # x$FiveUTR <-
#       #   reference %>%
#       #   plyranges::filter(transcript_id %in% x$canonical_transcript,
#       #                     type == "five_prime_utr")
# 
#       # x$CDS <- unlist(getSeq(Hsapiens, dplyr::filter(reference, transcript_id %in% x$canonical_transcript, type == "CDS")))
# 
# 
#       # x$ThreeUTR <-
#       #   reference %>%
#       #   plyranges::filter(transcript_id %in% x$canonical_transcript,
#       #                     type == "three_prime_utr")
# 
#       return(x)
#     }
#     )
#   }


# Main --------------------------------------------------------------------

# Parentgenes from Gencode v10 and therefore many pseudogene IDs do not ##
# match current annotation and needs to be converted to ensembl IDs ##

GencodeV10_gene_id <- 
  data.frame(gene_name = ref$GenV10$gene_name, 
             gene_id = ref$GenV10$gene_id) %>% 
  unique()

parentgenes$"Pseudo name" <- 
  GencodeV10_gene_id$gene_id[match(parentgenes$Name, GencodeV10_gene_id$gene_name)] 

# Make list of genes to compare
genes <- 
  split_and_index(data = parentgenes,
                  query = "Parent gene",
                  subject = "Pseudo name")

# Make sure list elements of equal size
if(nrow(genes$query) != nrow(genes$subject))
  stop("Number of genes to compare are not matching!")


# Get ensembl canonical transcripts and remove those where a canonical transcript is no found
genes_canonical <- 
  get_ensembl_canonical(gene_list = genes) %>% 
  
  # Remove genes missing canonical transcript
  remove_missing_canonical()

## Since we are interested in how pseudogenes affect multimapping I will use the MANE select ##
## for each parent gene and perform pairwisealignment with the complete pseudogene locus ##

# Clean up ref
# ensembl <- 
#   keepSeqlevels(
#     addchr(ref$ensV105), 
#     pruning.mode = "coarse",
#     value = c("chr1",
#               "chr2",
#               "chr3",
#               "chr4",
#               "chr5",
#               "chr6",
#               "chr7",
#               "chr8",
#               "chr9",
#               "chr10",
#               "chr11",
#               "chr12",
#               "chr13",
#               "chr14",
#               "chr15",
#               "chr16",
#               "chr17",
#               "chr18",
#               "chr19",
#               "chr20",
#               "chr21",
#               "chr22",
#               "chrX",
#               "chrY",
#               "chrMT"))

ensembl <- keepStandardChromosomes(ref$ensV105, pruning.mode = "coarse")


seq_similarity <- function(list, reference, PID) {
  
  # columns for PID output
  list$query[ ,c("PID1", "PID2", "PID3", "PID4")] <- NA 
  
  
  
  
}

for (i in seq_along(genes_canonical$query$index)) {
  
  # Transcript ID for parent MANE select
  TX <- genes_canonical$query$canonical_transcript[genes_canonical$query$index[i]]
  
  # Pseudo gene ID
  Gene <- genes_canonical$subject$`Pseudo name`[genes_canonical$query$index[i]]
  
  # Get sequences from BSgenome.Hsapiens.UCSC.hg38
  unlist(getSeq(Hsapiens, dplyr::filter(ensembl, transcript_id == TX, type == "CDS"))) -> query
  getSeq(Hsapiens, dplyr::filter(ensembl, gene_id == "ENSG00000160766", type == "gene")) -> subject
  
  # Pairwise alignment
  pairwiseAlignment(query, subject) -> align
  
  # Percent Sequence Identity; PID2 = 100 * (identical positions) / (aligned positions)
  genes_canonical$query$similarity[i] <- pid(align, type = "PID4")
  
}

# Save data ---------------------------------------------------------------
