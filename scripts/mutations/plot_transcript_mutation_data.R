## ClinVar data was downloaded and filtered using a bash getClinVarForLoci.sh avilable here: https://github.com/egustavsson/long-read_scripts
## GBA PD data was downloaded from: https://pdgenetics.shinyapps.io/gba1browser/

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(GenomicRanges)
library(GenomicFeatures)
library(ggtranscript)
library(R.utils)
library(rtracklayer)
library(VariantAnnotation)
library(cowplot)

# Arguments ---------------------------------------------------------------

args <-
  list(
    GBA1_browser = here::here("data", "GBA1-PD_browser_data", "GBA1-PD_browser_data-2023-04-03.txt"),
    GTF = here::here("data", "Brain", "Brain_regions_corrected.gtf.cds.gff"),
    clinvar = here::here("data", "ClinVar", "GBA1_clinvar.vcf.gz"))

# Load data ---------------------------------------------------------------

# ClinVar data
clinvar_data <- readVcf(args$clinvar)

# GBA1 browser data https://pdgenetics.shinyapps.io/gba1browser/
PD_browser_data = read.table(args$GBA1_browser, header = T, sep = "\t") %>% 
    dplyr::mutate(strand = "-",
                  seqnames = "chr1",
                  start = Position_hg38,
                  end = Position_hg38,
                  Variant_severity = tolower(Variant.type..Severe..Mild..Risk.Factor.),
                  Variant_severity = case_when(Variant_severity == "risk variant/risk variant" ~ "risk variant",
                                               TRUE ~ Variant_severity)) %>%
  drop_na(Position_hg38) %>% 
  dplyr::filter(Variant_severity != "unknown")

# long-read GTF file with GBA1 transcript structures
lr <- rtracklayer::import(args$GTF)

# also download and load reference annotation 
ref_path <- here::here(tempdir(), "gencode.v38.annotation.gtf.gz")

if(!file.exists(ref_path)) {
  
  download.file(
    url = paste0(
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/", 
      "gencode.v38.annotation.gtf.gz"
    ),
    destfile = ref_path
  )
  
}

R.utils::gunzip(ref_path, remove = TRUE)

ref <- rtracklayer::import(stringr::str_remove(ref_path, "\\.gz"))

# Functions ---------------------------------------------------------------

# Function to convert a collapsed VCF file to a DF and filter for pathogenic variants
PathogenicVarFromVCF <- function(vcf_file, filter_pathogenic = TRUE) {
  
  ranges <- as.data.frame(rowRanges(vcf_file, fixed = F)) %>% rownames_to_column(var = "id")
  
  info <- info(vcf_file)[c("CLNDN", "CLNSIG")] %>% as.data.frame() %>% rownames_to_column(var = "id")
  
  joined_df <- left_join(ranges, info, by = "id")
  
  if (filter_pathogenic == TRUE) {
    pathogenic_strings <- paste(c("Likely_pathogenic",
                                  "Pathogenic",
                                  "Pathogenic/Likely_pathogenic",
                                  "Pathogenic|risk_factor",
                                  "Pathogenic/Likely_pathogenic|risk_factor"),
                                collapse="|")
    
    joined_df <- joined_df %>% dplyr::filter(grepl(pathogenic_strings, CLNSIG))
    
  }
  
  return(joined_df)
}

get_lr_tx_of_interest <- function(lr, pb_ids) {
  
  lr <- lr[lr$transcript_id %in% pb_ids]
  lr_exons <- lr[lr$type == "exon"]
  lr_cds <- lr[lr$type == "CDS"]
  
  # return as list as we need both exons and cds
  lr_exons_cds <- list(
    exons = lr_exons, 
    cds = lr_cds
  )
  
  return(lr_exons_cds)
  
}

get_mane <- function(ref, mane_id) {
  
  # remove any NA transcript ids (i.e. type == "gene")
  mane <- ref[!is.na(ref$transcript_id)] 
  
  # remove to .XX after ENST
  GenomicRanges::mcols(mane)[["transcript_id"]] <- 
    GenomicRanges::mcols(mane)[["transcript_id"]] %>% 
    stringr::str_remove("\\..*")
  
  mane <- mane[mane$transcript_id == mane_id, ]
  mane_exons <- mane[mane$type == "exon"]
  mane_cds <- mane[mane$type == "CDS"]
  
  mane_exons_cds <- list(
    exons = mane_exons, 
    cds = mane_cds
  )
  
  return(mane_exons_cds)
  
}

plot_transcripts_and_mutations <- function(lr_exons_cds, 
                                           mane_exons_cds, 
                                           lr_mane_diffs, 
                                           MANE_canonical = "MANE",
                                           seqnames, 
                                           start, 
                                           end, 
                                           strand,
                                           mutation_data,
                                           clinvar_data
) {
  
  locus_subset <- GRanges(seqnames = seqnames, ranges = IRanges(start = start, end = end), strand = strand)
  
  # merge mane and lr data and convert to data.frame() for plotting
  # convert transcript_id to factor to make sure mane is at top
  transcript_order <- c(
    lr_exons_cds$exons$transcript_id %>% unique(),
    mane_exons_cds$exons$transcript_id %>% unique()
  )
  
  lr_mane_exons_df <- c(lr_exons_cds$exons, mane_exons_cds$exons) %>% 
    as.data.frame() %>% 
    dplyr::mutate(
      transcript_id = transcript_id %>% 
        factor(levels = transcript_order)
    )
  lr_mane_cds_df <- c(lr_exons_cds$cds, mane_exons_cds$cds) %>% 
    as.data.frame() %>% 
    dplyr::mutate(
      transcript_id = transcript_id %>% 
        factor(levels = transcript_order)
    )
  
  # plot diff plot
  diff_plot <- lr_mane_exons_df %>% 
    ggplot(aes(
      xstart = start, 
      xend = end, 
      y = transcript_id
    )) + 
    geom_range(
      height = 0.25, 
      fill = "white"
    ) + 
    geom_range(
      data = lr_mane_cds_df
    ) + 
    geom_intron(
      data = to_intron(lr_mane_exons_df, "transcript_id"), 
      aes(strand = strand), 
      arrow.min.intron.length = 400
    ) + 
    geom_range(
      data = lr_mane_diffs, 
      aes(
        fill = diff_type, 
        colour = diff_type
      ), 
      alpha = 0.2, 
      linetype = 2
    ) + 
    labs(y = "Transcript name",
         x = "") +
    xlim(start(locus_subset), end(locus_subset)) +
    scale_fill_discrete(
      name = paste0("Region in ", MANE_canonical, " transcript:"), 
      labels = c(
        paste0("In ", MANE_canonical), 
        paste0("Not in ", MANE_canonical)
      )
    ) + 
    scale_colour_discrete(
      name = paste0("Region in ", MANE_canonical, " transcript:"), 
      labels = c(
        paste0("In ", MANE_canonical), 
        paste0("Not in ", MANE_canonical)
      )
    ) + 
    theme_bw() + 
    theme(legend.position = "top",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  # plot mutation track
  mutation_plot <-
    mutation_data %>% 
    ggplot(aes(x = start)) +  
    geom_histogram(bins = 200, color = "black", fill = "grey") +
    xlim(start(locus_subset), end(locus_subset)) +
    xlab("Position (hg38)") +
    facet_grid(vars(title, factor(Variant_severity, levels = c("Pathogenic", "severe (null)", "severe", "mild", "risk variant", "unknown")))) +
    theme_bw() +
    theme(panel.spacing.y = unit(0, "cm"), legend.position = "none")
  
  final_plot <-
    plot_grid(diff_plot, 
              mutation_plot,
              ncol = 1, 
              align = "hv", 
              rel_heights = c(0.8, 2), 
              axis = "lr")
  
  return(final_plot)
  
}


# Main --------------------------------------------------------------------

# clinvar
clinvar_df <- PathogenicVarFromVCF(vcf_file = clinvar_data, filter_pathogenic = TRUE)
  
# join clinvar and GBA1 PD browser data
mutation_data <-
  rbind(
    PD_browser_data %>%
      dplyr::mutate(title = "Parkinson's Disease") %>%
      dplyr::select(start, Variant_severity, title),
    clinvar_df %>%
      dplyr::mutate(title = "ClinVar",
                    Variant_severity = "Pathogenic") %>%
      dplyr::select(start, Variant_severity, title)
  )

# get transcripts of interest for ggtranscript plot
gba_lr_exons_cds <- get_lr_tx_of_interest(
  lr = lr, 
  pb_ids = c(
    "PB.845.2627", 
    "PB.845.2629",
    "PB.845.2954"
  ))

# MANE select transcript structure
gba_mane_exons_cds <- get_mane(ref, mane_id = "ENST00000368373")

# obtain differences between MANE and lr exons
gba_lr_mane_diffs <- 
  ggtranscript::to_diff(
    exons = gba_lr_exons_cds$exons %>% as.data.frame(), 
    ref_exons = gba_mane_exons_cds$exons %>% as.data.frame(), 
    group_var = "transcript_id"
  )

# plot transcripts with mutation data annotated
transcript_and_mutation_plot <- 
  plot_transcripts_and_mutations(
    lr_exons_cds = gba_lr_exons_cds, 
    mane_exons_cds = gba_mane_exons_cds, 
    lr_mane_diffs = gba_lr_mane_diffs, 
    seqnames = "chr1", 
    start = 155234000, 
    end = 155242000, 
    strand = "-", 
    mutation_data = mutation_data)

# Save data ---------------------------------------------------------------

ggsave(
  plot = transcript_and_mutation_plot, 
  filename = "transcript_and_mutation_plot.svg", 
  path = here::here("results", "mutations"), 
  width = 10, 
  height = 10, 
  dpi = 600
)
