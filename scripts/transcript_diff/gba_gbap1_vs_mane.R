
# Load libraries ----------------------------------------------------------

# as ggtranscript is in development during the creation of this script
# i have included the code to install the exact version using (0.99.2)
# based on the commit SHA
# devtools::install_github("dzhang32/ggtranscript@88d23b43b42de3e49a67dc8130f5263b6dcf81d1")
library(ggtranscript)

library(tidyverse)
library(here)
library(R.utils)
library(rtracklayer)

# Load data ---------------------------------------------------------------

lr_data <- rtracklayer::import(
  here::here("raw_data", "Brain_regions_corrected.gtf.cds.gff")
  )

# also download and load reference annotation 
gtf_path <- here::here(tempdir(), "gencode.v38.annotation.gtf.gz")

if(!file.exists(gtf_path)) {
  
  download.file(
    url = paste0(
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/", 
      "gencode.v38.annotation.gtf.gz"
    ),
    destfile = gtf_path
  )
  
}
  
R.utils::gunzip(gtf_path, remove = TRUE)

gtf <- rtracklayer::import(stringr::str_remove(gtf_path, "\\.gz"))

# Main --------------------------------------------------------------------

##### GBA #####

# filter for txs of interest
gba_lr <- lr_data[lr_data$gene_id %>% stringr::str_detect("ENSG00000177628")]

gba_txs_of_interest <- c(
  "PB.845.2888",
  "PB.845.2786", 
  "PB.845.2627", 
  "PB.845.2629",
  "PB.845.2848",
  "PB.845.2954"
)

gba_lr <- gba_lr[gba_lr$transcript_id %in% gba_txs_of_interest]
gba_lr_exons <- gba_lr[gba_lr$type == "exon"]
gba_lr_cds <- gba_lr[gba_lr$type == "CDS"]

# obtain MANE select transcript
gba_mane <- gtf[!is.na(gtf$transcript_id)] 
GenomicRanges::mcols(gba_mane)[["transcript_id"]] <- 
  GenomicRanges::mcols(gba_mane)[["transcript_id"]] %>% 
  stringr::str_remove("\\..*")
gba_mane <- gba_mane[gba_mane$transcript_id == "ENST00000368373", ]
gba_mane_exons <- gba_mane[gba_mane$type == "exon"]
gba_mane_cds <- gba_mane[gba_mane$type == "CDS"]

# obtain differences between MANE and lr exons
gba_lr_mane_diffs <- 
  ggtranscript::to_diff(
  exons = gba_lr_exons %>% as.data.frame(), 
  ref_exons = gba_mane_exons %>% as.data.frame(), 
  group_var = "transcript_id"
)

# merge mane and lr data and convert to data.frame() for plotting
gba_lr_mane_exons <- c(gba_lr_exons, gba_mane_exons) %>% as.data.frame()
gba_lr_mane_cds <- c(gba_lr_cds, gba_mane_cds) %>% as.data.frame()

# plot diff plot
gba_lr_mane_diff_plot <- gba_lr_mane_exons %>% 
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
    data = gba_lr_mane_cds
  ) + 
  geom_intron(
    data = to_intron(gba_lr_mane_exons, "transcript_id"), 
    aes(strand = strand), 
    arrow.min.intron.length = 400
  ) + 
  geom_range(
    data = gba_lr_mane_diffs, 
    aes(
      fill = diff_type, 
      colour = diff_type
      ), 
    alpha = 0.2, 
    linetype = 2
  ) + 
  scale_y_discrete(name = "Transcript ID") + 
  scale_x_continuous(name = "Genomic position") + 
  scale_fill_discrete(
    name = "Region in MANE-select transcript", 
    labels = c("In MANE", "Not in MANE")
  ) + 
  scale_colour_discrete(
    name = "Region in MANE-select transcript", 
    labels = c("In MANE", "Not in MANE")
  ) + 
  theme_bw() + 
  theme(legend.position = "top")

# Save data ---------------------------------------------------------------

ggsave(
  plot = gba_lr_mane_diff_plot, 
  filename = "gba_lr_mane_diff_plot.png", 
  path = here::here("results", "transcript_diff"), 
  width = 6, 
  height = 4, 
  dpi = 600
)
