#! usr/bin/env Rscript

library(tidyverse)
library(edgeR)

#' Takes input count data and grouping vector, returning a filtered count data frame.
#' 
#' @param count_data matrix or data frame: A (G x S) counts matrix of gene expression data
#' with S sample columns, G rows of genes, and gene names as row names.
#' 
#' @param group_vector group_vector: A (1 x S) vector by which samples are grouped.
#' 
#' @returns tibble: A (N x 1+S) tibble with N rows of filtered genes, a column of gene names and S sample columns.
filter_count_matrix <- function(count_data, group_vector){
  # Generate DGEList object from counts data.
  unfiltered_dge <- edgeR::DGEList(counts = count_data, group = factor(group_vector))
  
  # Filter DGEList of counts by expression, grouped by vector.
  filtered_matrix <- unfiltered_dge[edgeR::filterByExpr(unfiltered_dge),
                                    keep.lib.sizes=FALSE] %>%
    
    # Normalize library size.
    edgeR::normLibSizes(method = "TMM") %>%
    
    # Normalize gene expression by taking log(counts per million).
    edgeR::cpm(log = TRUE) %>% 
    
    # Output as tibble with genes as first column.
    tibble::as_tibble(rownames = 'gene')
}