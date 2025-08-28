# Make new volcano plot diagrams
library(ggplot2)
library(ggthemes)
library(ggrepel)

# Code from: https://biostatsquid.com/volcano-plots-r-tutorial/

dgrp42=read.csv("rawdata/S42_LC_D0vsS42_HC_D0_deg.csv")
dgrp42$diffexpressed = "NO"
dgrp42$diffexpressed[dgrp42$log2FoldChange>0 & dgrp42$pvalue<0.05]<- "UP"
length(dgrp42$diffexpressed[dgrp42$diffexpressed=="UP"])
dgrp42$diffexpressed[dgrp42$log2FoldChange<0 & dgrp42$pvalue<0.05] <-"DOWN"
length(dgrp42$diffexpressed[dgrp42$diffexpressed=="DOWN"])
length(dgrp42$diffexpressed[dgrp42$diffexpressed=="NO"])

dgrp42$delabel=ifelse(dgrp42$gene_name %in% head(dgrp42[order(dgrp42$pvalue), "gene_name"],30),dgrp42$gene_name,NA)

volcanoplot42 <- ggplot(data = dgrp42, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed,label=delabel)) +
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01301), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) + theme_base() +
  scale_color_manual(values = c("#a63603", "grey", "#fdd0a2")) + # to set the colours of our variable      
  coord_cartesian(ylim = c(0, 22), xlim = c(-10, 14)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = '', #legend_title,
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
  scale_x_continuous(breaks = seq(-10, 14, 10))+scale_y_continuous(breaks = seq(0, 22,10)) + # to customise the breaks in the x axis
  geom_text_repel(max.overlaps = 10) # To show all labels 
volcanoplot42

ggsave("images/DGRP_42_volcano.png",plot=volcanoplot42)


dgrp217=read.csv("rawdata/S217_LC_D0vsS217_HC_D0_deg.csv")
dgrp217$diffexpressed = "NO"
dgrp217$diffexpressed[dgrp217$log2FoldChange>0 & dgrp217$pvalue<0.05]<- "UP"
length(dgrp217$diffexpressed[dgrp217$diffexpressed=="UP"])
dgrp217$diffexpressed[dgrp217$log2FoldChange<0 & dgrp217$pvalue<0.05] <-"DOWN"
length(dgrp217$diffexpressed[dgrp217$diffexpressed=="DOWN"])
length(dgrp217$diffexpressed[dgrp217$diffexpressed=="NO"])

dgrp217$delabel=ifelse(dgrp217$gene_name %in% head(dgrp217[order(dgrp217$pvalue), "gene_name"],30),dgrp217$gene_name,NA)

volcanoplot217 <- ggplot(data = dgrp217, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed,label=delabel)) +
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01301), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) + theme_base() +
  scale_color_manual(values = c("#08519c", "grey", "#c6dbef")) + # to set the colours of our variable      
  coord_cartesian(ylim = c(0, 22), xlim = c(-10, 14)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = '', #legend_title,
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
  scale_x_continuous(breaks = seq(-10, 14, 10))+scale_y_continuous(breaks = seq(0, 22,10)) + # to customise the breaks in the x axis
  geom_text_repel(max.overlaps = 10) # To show all labels 
volcanoplot217

ggsave("images/DGRP_217_volcano.png",plot=volcanoplot217)




### Now to try to re-plot the GO analysis dot-plots
# Used Claude to generate code below:

#data output from Novogene:
go_results=read.csv("rawdata/S42_LC_D0vsS42_HC_D0.all_GOenrich_significant.csv")
kegg=read.csv("rawdata/S42_LC_D0vsS42_HC_D0.all_KEGGenrich_significant.csv")

# Load required libraries
library(ggplot2)
library(dplyr)

# Function to create dotplot for enrichment results
create_enrichment_dotplot <- function(data, title = "Enrichment Analysis", 
                                      max_terms = 20, description_col = "Description") {
  
  # Parse GeneRatio to numeric values
  data$GeneRatio_numeric <- sapply(strsplit(data$GeneRatio, "/"), function(x) {
    as.numeric(x[1]) / as.numeric(x[2])
  })
  
  # Order by p-adjusted value and take top terms
  data_plot <- data %>%
    arrange(padj) %>%
    head(max_terms) %>%
    # Reorder descriptions by GeneRatio for better visualization
    mutate(!!description_col := reorder(get(description_col), GeneRatio_numeric))
  
  # Create the dotplot
  p <- ggplot(data_plot, aes(x = GeneRatio_numeric, y = get(description_col))) +
    geom_point(aes(size = Count, color = padj)) +
    scale_color_gradient(low = "red", high = "blue", 
                         name = "Adjusted\np-value",
                         trans = "log10") +
    scale_size_continuous(name = "Gene\nCount", range = c(2, 10)) +
    labs(x = "Gene Ratio", 
         y = "Description",
         title = title) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9)) +
    # Add guidelines
    theme(panel.grid.major = element_line(color = "grey90", size = 0.5),
          panel.grid.minor = element_line(color = "grey95", size = 0.3))
  
  return(p)
}

# Function specifically for GO results
create_GO_dotplot <- function(go_data, max_terms = 20, category_filter = NULL) {
  
  # Filter by category if specified (e.g., "BP", "MF", "CC")
  if (!is.null(category_filter)) {
    go_data <- go_data[go_data$Category %in% category_filter, ]
  }
  
  # Create the plot
  plot <- create_enrichment_dotplot(go_data, 
                                    title = "GO Enrichment Analysis", 
                                    max_terms = max_terms,
                                    description_col = "Description")
  
  return(plot)
}

# Function specifically for KEGG results
create_KEGG_dotplot <- function(kegg_data, max_terms = 20) {
  
  plot <- create_enrichment_dotplot(kegg_data, 
                                    title = "KEGG Pathway Enrichment Analysis", 
                                    max_terms = max_terms,
                                    description_col = "Description")
  
  return(plot)
}

# Example usage (replace 'go_results' and 'kegg_results' with your actual data frames):

# For GO results:
 go_plot <- create_GO_dotplot(go_results)
 print(go_plot)

# For GO results filtered by category (e.g., only Biological Process):
# go_bp_plot <- create_GO_dotplot(go_results, max_terms = 15, category_filter = "BP")
# print(go_bp_plot)

# For KEGG results:
 kegg_plot <- create_KEGG_dotplot(kegg)
 print(kegg_plot)

# Save plots
 ggsave("images/GO_enrichment_dotplot.png", go_plot, width = 10, height = 8, dpi = 300)
 ggsave("images/KEGG_enrichment_dotplot.png", kegg_plot, width = 10, height = 8, dpi = 300)

# Advanced: Create faceted plot for GO categories
create_GO_faceted_dotplot <- function(go_data, max_terms_per_category = 10) {
  
  # Get top terms per category
  data_plot <- go_data %>%
    group_by(Category) %>%
    arrange(padj) %>%
    slice_head(n = max_terms_per_category) %>%
    ungroup()
  
  # Parse GeneRatio to numeric
  data_plot$GeneRatio_numeric <- sapply(strsplit(data_plot$GeneRatio, "/"), function(x) {
    as.numeric(x[1]) / as.numeric(x[2])
  })
  
  # Reorder within each category
  data_plot <- data_plot %>%
    group_by(Category) %>%
    mutate(Description = reorder_within(Description, GeneRatio_numeric, Category)) %>%
    ungroup()
  
  # Create faceted plot
  p <- ggplot(data_plot, aes(x = GeneRatio_numeric, y = Description)) +
    geom_point(aes(size = Count, color = padj)) +
    scale_color_gradient(low = "red", high = "blue", 
                         name = "Adjusted\np-value",
                         trans = "log10") +
    scale_size_continuous(name = "Gene\nCount", range = c(2, 8)) +
    scale_y_reordered() +
    facet_wrap(~Category, scales = "free_y", ncol = 1) +
    labs(x = "Gene Ratio", 
         y = "GO Term Description",
         title = "GO Enrichment Analysis by Category") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          strip.text = element_text(size = 10, face = "bold"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9))
  
  return(p)
}

# Helper functions for the faceted plot
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}

# Example usage for faceted GO plot:
 go_faceted_plot <- create_GO_faceted_dotplot(go_results)
 print(go_faceted_plot)
 ggsave("images/GO_enrichment_faceted_dotplot.png", go_faceted_plot, width = 12, height = 10, dpi = 300)

