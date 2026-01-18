#' =============================================================================
#' Helper Functions for Niche Breadth Analysis
#' =============================================================================
#'
#' Utility functions used across simulation and metric calculation scripts.
#' =============================================================================

#' Set Max and Min Values
#'
#' Scales a vector to a new range defined by new_max and new_min.
#'
#' @param x Numeric vector to scale
#' @param new_max New maximum value
#' @param new_min New minimum value
#' @return Scaled numeric vector
set_max_min <- function(x, new_max, new_min) {
  current_max <- max(x)
  current_min <- min(x)
  current_range <- current_max - current_min
  new_range <- new_max - new_min

  x_scaled <- (x - current_min) * new_range / current_range + new_min
  x_scaled[x == current_max] <- new_max
  x_scaled[x == current_min] <- new_min

  return(x_scaled)
}


#' Weighted Variance
#'
#' Calculates the weighted variance for each column of a matrix.
#'
#' @param x Matrix or data.frame of values
#' @param w Vector of weights (same length as number of rows in x)
#' @return Vector of weighted variances (one per column)
weighted.variance <- function(x, w) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    x <- as.matrix(x)
  }

  w.var <- numeric(ncol(x))

  for (i in 1:ncol(x)) {
    w.mean <- sum(w * x[, i]) / sum(w)
    w.ss <- sum(w * (x[, i] - w.mean)^2)
    w.var[i] <- w.ss / sum(w)
  }
  return(w.var)
}


#' Eliminate Columns by Sum Range
#'
#' Filters columns of a data frame based on column sum thresholds.
#' Useful for removing species that are too rare or too common.
#'
#' @param df Data frame (species as columns)
#' @param min_sum Minimum column sum to keep
#' @param max_sum Maximum column sum to keep
#' @return Data frame with filtered columns
eliminate_cols <- function(df, min_sum, max_sum) {
  col_sums <- colSums(df)
  keep_cols <- which(col_sums >= min_sum & col_sums <= max_sum)
  return(df[, keep_cols, drop = FALSE])
}


#' Melt Data Frame for Niche Analysis
#'
#' Converts wide-format species data to long format for hypervolume calculations.
#'
#' @param data Species presence-absence matrix (sites x species)
#' @param env_vars Environmental variables matrix (sites x n_env)
#' @param species_cols Character vector of species column names
#' @return Long-format data frame with species, environmental variables
melt_species_data <- function(data, env_vars, species_cols = NULL) {
  if (is.null(species_cols)) {
    species_cols <- colnames(data)
  }

  df_combined <- cbind(as.data.frame(data), as.data.frame(env_vars))
  env_cols <- colnames(env_vars)

  df_melted <- tidyr::gather(df_combined, key = "species", value = "PA",
                              all_of(species_cols))

  # Rearrange columns and filter
  df_melted <- df_melted[, c("species", "PA", env_cols)]
  df_melted <- df_melted[df_melted$PA != 0, ]
  df_melted <- df_melted[, !(colnames(df_melted) %in% "PA")]

  return(df_melted)
}


#' Aggregate Niche Breadth Results
#'
#' Aggregates metric values across iterations by computing means per species.
#'
#' @param df Data frame with iteration, sci.name, and metric columns
#' @param metric_cols Character vector of metric column names
#' @return Aggregated data frame (one row per species)
aggregate_niche_breadth <- function(df, metric_cols = NULL) {
  if (is.null(metric_cols)) {
    # Auto-detect metric columns (exclude iteration and sci.name)
    metric_cols <- setdiff(colnames(df), c("iteration", "sci.name", "niche_breadth"))
  }

  result <- df %>%
    dplyr::group_by(sci.name) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(metric_cols),
                                    ~mean(.x, na.rm = TRUE)))

  return(result)
}


#' Scale Data Frame Columns
#'
#' Scales specified columns of a data frame to have mean 0 and sd 1.
#'
#' @param df Data frame to scale
#' @param cols Columns to scale (default: all numeric except first two)
#' @return Data frame with scaled columns
scale_metrics <- function(df, cols = NULL) {
  if (is.null(cols)) {
    # Default: scale all numeric columns except first two (iteration, sci.name)
    cols <- 3:ncol(df)
  }

  df[, cols] <- scale(df[, cols])
  return(df)
}


#' Calculate Spearman Correlation with Oracle
#'
#' Computes Spearman correlation between each metric and the oracle niche breadth.
#'
#' @param df Data frame with niche_breadth and metric columns
#' @param oracle_col Name of the oracle niche breadth column
#' @param metric_cols Names of metric columns to correlate
#' @return Data frame with metric names and correlations
calculate_metric_correlations <- function(df, oracle_col = "niche_breadth",
                                           metric_cols = NULL) {
  if (is.null(metric_cols)) {
    metric_cols <- setdiff(colnames(df), c("iteration", "sci.name", oracle_col))
  }

  correlations <- sapply(metric_cols, function(col) {
    cor(df[[oracle_col]], df[[col]], method = "spearman", use = "complete.obs")
  })

  result <- data.frame(
    metric = metric_cols,
    correlation = correlations,
    stringsAsFactors = FALSE
  )

  return(result)
}


#' Flatten Correlation Matrix
#'
#' Extracts the lower triangle of a correlation matrix as a vector.
#'
#' @param mat Correlation matrix
#' @return Vector of lower-triangle values
flatten_corr <- function(mat) {
  mat[lower.tri(mat, diag = FALSE)]
}


#' Calculate Coefficient of Variation
#'
#' Computes CV for each metric in a data frame.
#'
#' @param df Data frame with metric columns
#' @param metric_cols Names of metric columns
#' @return Named vector of CVs
calculate_cv <- function(df, metric_cols) {
  sapply(metric_cols, function(m) {
    x <- df[[m]]
    mval <- mean(x, na.rm = TRUE)
    if (is.na(mval) || mval == 0) NA else sd(x, na.rm = TRUE) / mval
  })
}


#' Safe Print
#'
#' Prints a message with optional verbosity control.
#'
#' @param msg Message to print
#' @param verbose Logical; if FALSE, suppresses output
safe_print <- function(msg, verbose = TRUE) {
  if (verbose) {
    message(msg)
  }
}
