#! usr/bin/env Rscript

library(tidyverse)

#' Calculate f-statistic using an ANOVA model.
#' 
#' @param data_matrix matrix (G x S): matrix with S columns as samples columns and 
#' G rows as features.
#' 
#' @param alt matrix (S x V): matrix with rows as samples and columns containing
#' the variables for alternative model in the ANOVA analysis.
#' 
#' @param null matrix (S x V-1): matrix with rows as samples and columns containing
#' the variables for null model in the ANOVA analysis.
#' 
#' @returns numeric (G x 1): vector containing f-statistic value for each feature.
fstats <- function(data_matrix, alt, null) {
  n <- ncol(data_matrix)
  df_alt <- ncol(alt)
  df_null <- ncol(null)
  
  # Compute projection matrices.
  proj_alt <- alt %*% solve(crossprod(alt))
  proj_null <- null %*% solve(crossprod(null))
  
  # Compute residuals.
  residual_alt <- data_matrix %*% (diag(n) - proj_alt %*% t(alt))
  residual_null <- data_matrix %*% (diag(n) - proj_null %*% t(null))
  
  # Compute residual sum of squares.
  rss_alt <- rowSums(residual_alt^2)
  rss_null <- rowSums(residual_null^2)
  
  # Compute F-statistic.
  as.numeric(((rss_null - rss_alt)/(df_alt - df_null))/(rss_alt/(n - df_alt)))
}

#' Calculate the f-statistic for multiple groups.
#' 
#' @param data_matrix matrix (G x S): matrix with S columns as samples columns and 
#' G rows as features.
#' 
#' @param full_model data frame (S x V): data frame with rows as samples and columns containing
#' the variables used in the ANOVA analysis. Column names should end in unique digits.
#' 
#' @returns tibble (G*V x 3): tibble in long format with three columns: gene name, f-statistic,
#' and group tested in ANOVA model for the given gene.
calculate_fstat <- function(data_matrix, full_model){
  # Helper function to remove a column from full model to produce reduced model.
  deselect <- function(full_model, removed){
    full_model %>%
      as.data.frame() %>% 
      dplyr::select(!all_of(removed)) %>% 
      as.matrix()
  }
  
  # Add a column of ones to full model as the intercept.
  full <- full_model %>% dplyr::mutate(Intercept = 1, .before = everything()) %>% as.matrix()
  
  # Iterate over each column of the input data frame used for the full model,
  # removing it from the reduced model to calculate the f-statistic.
  purrr::map(names(full_model), ~ fstats(data_matrix, full, deselect(full, .x))) %>%
    tibble::as_tibble(.name_repair = 'unique_quiet') %>%
    
    # Set the column names of the output to the columns used for the full model.
    dplyr::rename_with(~ colnames(full_model), .cols = everything()) %>%
    
    # Add row names of data_matrix as a column to output.
    dplyr::mutate(gene = rownames(data_matrix), .before = everything()) %>% 
    
    # Pivot from wide to long, setting the cluster column to integer.
    tidyr::pivot_longer(!gene, names_to = 'cluster', values_to = 'fstat',) %>%
    #names_transform = list(cluster = as.integer)) %>%
    
    # Reorder output.
    dplyr::select(gene, fstat, cluster) %>%
    dplyr::arrange(cluster)
}

#' Calculate the beta coefficients for linear model of gene expression.
#' 
#' @param expression_data matrix (G x S): matrix with S columns as samples columns and 
#' G rows as features.
#' 
#' @param variables data frame (S x V): data frame with rows as samples and columns containing
#' the variables used in the ANOVA analysis. Column names should end in unique digits.
#' 
#' @returns tibble (G*V x 4): tibble in long format with three columns: gene name, beta, sign
#' indicating whether beta > or < 0, and variable of beta coefficient.
calculate_beta_coefficient <- function(expression_data, variables) {
  # Helper function to produce beta coefficients from a linear model.
  beta_coefficient <- function(a, b){
    linearmod <- lm(a ~ ., data = b)
    
    # Drop intercept value from output.
    lm(data.frame(scale(linearmod$model)))$coefficients[-1]
  }
  
  # Iterate over expression data rows, calculating beta coefficient of gene expression
  # as a linear model of model variables.
  apply(expression_data, MARGIN = 1, function(x) beta_coefficient(x, variables)) %>%
    tibble::as_tibble(rownames = 'cluster') %>%
    tidyr::pivot_longer(!cluster, names_to = 'gene', values_to = 'beta') %>% 
    
    # Add sign column indicating whether beta coefficient is positive or negative.
    dplyr::mutate(sign = ifelse(beta > 0, 1, -1)) %>%
    
    # Reorder output.
    dplyr::select(gene, beta, sign, cluster) %>% 
    dplyr::arrange(cluster)
}

#' Combine f-statistic values with sign of beta coefficient to rank genes
#' according to correlation.
signed_fstat <- function(expression_data, variables) {
  dplyr::left_join(x = calculate_fstat(expression_data, variables),
                   y = calculate_beta_coefficient(expression_data, variables),
                   by = c('gene', 'cluster')) %>%
    dplyr::mutate(sign_fstat = fstat*sign) %>%
    dplyr::select(gene, sign_fstat, cluster)
}

#' Optional function to add a small amount of noise to break ties to input vector.
#' Noise is proportional to sizes of values in input, to minimize likelihood of altering
#' ranking position with neighboring values.
add_noise <- function(input_vector, seed = NULL) {
  set.seed(seed)
  noise <- runif(length(input_vector), min=0, max=0.5) %>%
    round(digits = 3)  
  input_vector + (noise*input_vector/100000000)
}
