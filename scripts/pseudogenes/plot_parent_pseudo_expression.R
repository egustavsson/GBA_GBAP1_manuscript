# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)

# Load data ---------------------------------------------------------------

load("/home/egust/Projects/pseudogenes/results/pseudogene_annotated.rda")
load("/home/egust/Projects/pseudogenes/results/parentgene_annotated.rda")

# Main --------------------------------------------------------------------

parent_pseudo <-
  list(
    "Parent genes" = parentgene_annotated %>%
      dplyr::filter(associated_pseudo %in% c("Unprocessed", "Processed")) %>%
      dplyr::select(no_tissues_expressed,
                    type = associated_pseudo) %>%
      dplyr::mutate(gene = "Parent genes") %>% 
      na.omit(),
    
    "Pseudogenes" = pseudogene_annotated %>%
      dplyr::filter(pseudogene_type %in% c("Unprocessed", "Processed")) %>%
      dplyr::select(no_tissues_expressed,
                    type = pseudogene_type) %>%
      dplyr::mutate(gene = "Pseudogenes") %>% 
      na.omit()
    )


## KS test of distributions ## 

stat.data <- 
  lapply(parent_pseudo, 
         function(x) ks.test(x$no_tissues_expressed[x$type == "Unprocessed"], 
                             x$no_tissues_expressed[x$type == "Processed"])) %>% 
  
  tibble(
    gene = factor(c("Pseudogenes", "Parent genes")),
    Pvalue = c(ifelse(.$`Parent genes`$p.value < 0.05, "P < .05", paste0("P = ", .$`Parent genes`$p.value)),
             ifelse(.$`Pseudogenes`$p.value < 0.05, "P < .05", paste0("P = ", .$`Parent genes`$p.value))),
    Dvalue = c(paste0("D = ", round(.$`Parent genes`$statistic, digits = 3)),
               paste0("D = ", round(.$`Pseudogenes`$statistic, digits = 3))
               )) %>% 
  dplyr::select(Gene, Pvalue, Dvalue)



# ERs and exons with repeats

  geom_label(data = stat.data, aes(1.5, 0.85, label = stat), size = 2.3)
 


## Plot parent and pseudogene expression ##

parent_pseudo_to_plot <-
  bind_rows(!!!parent_pseudo) %>%
  
  ggplot(aes(x = no_tissues_expressed, 
             fill = type))+
  geom_histogram(aes(y=..density..),
                 position = "dodge", 
                 colour = "black", 
                 binwidth = 2) + 
  geom_text(aes(label = Dvalue),
            data = stat.data,
            vjust = "top",
            hjust = "right",
            inherit.aes = F) +
  labs(x = "Number of tissues expressed", y = "Fraction of genes") +
  scale_fill_manual(values = c("#7570B3", "#E6AB02")) +
  scale_x_continuous(n.breaks = 10) +
  facet_wrap(vars(
    factor(gene, 
           levels = c("Pseudogenes", "Parent genes")))) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(face = "bold",
                                   size = 8),
        axis.text.y = element_text(face = "bold",
                                   size = 10),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(face = "bold",
                                    size = 12),
        strip.background =element_rect(fill="gray80"),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.6))
  

# Save data ---------------------------------------------------------------

ggsave(plot = parent_pseudo_to_plot, 
       filename = "parent_pseudo_expression_plot.png", 
       path = here::here("results", "pseudogenes"), 
       width = 8, 
       height = 4, 
       dpi = 600
)
