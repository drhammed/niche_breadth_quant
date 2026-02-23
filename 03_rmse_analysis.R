#' =============================================================================
#' 03_rmse_analysis.R
#' Supplementary performance evaluation: Rank-based RMSE
#' =============================================================================
#'
#' Computes rank-based RMSE between oracle niche breadth and each metric
#' across all 12 simulation scenarios. Complements the Spearman correlation
#' analysis in 02_metrics_agreement.Rmd.
#'
#' Rank RMSE captures the magnitude of rank displacement: a metric could have
#' moderate Spearman rho but low rank RMSE if errors are small and evenly
#' distributed, or high rank RMSE if a few species are grossly misranked.
#' =============================================================================

# --- Packages ----------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, knitr)

# --- Configuration -----------------------------------------------------------
source("config.R")
source("R/helper_functions.R")

METRICS <- c("SimpSSI", "beta.a", "beta.w", "om_tol", "nr_hv",
             "hv_blond", "nb_Gam", "nb_latent", "nb_dist")

# --- Load results ------------------------------------------------------------
results_dir <- "results/simulations"
rds_files   <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)

all_results <- list()
for (file in rds_files) {
  obj_name <- tools::file_path_sans_ext(basename(file))
  all_results[[obj_name]] <- readRDS(file)
}

cat("Loaded", length(all_results), "scenario files\n")

# --- Rank-based RMSE function ------------------------------------------------
#' For a single iteration (n species), rank both oracle and metric values,
#' then compute RMSE between the two rank vectors.
#'
#' Metrics where higher values indicate narrower niches (inverse relationship)
#' are handled by using absolute rank differences --- the Spearman rho sign
#' already captures direction; rank RMSE captures displacement magnitude.

rank_rmse <- function(oracle, estimated) {
  # Remove NAs pairwise
  valid <- complete.cases(oracle, estimated)
  if (sum(valid) < 3) return(NA_real_)

  oracle_ranks   <- rank(oracle[valid])
  metric_ranks   <- rank(estimated[valid])
  # Also check the inverse ranking (for metrics inversely related to breadth)
  rmse_forward   <- sqrt(mean((oracle_ranks - metric_ranks)^2))
  rmse_inverse   <- sqrt(mean((oracle_ranks - (max(metric_ranks) + 1 - metric_ranks))^2))
  # Return the lower of the two (best-case alignment)
  min(rmse_forward, rmse_inverse)
}

# --- Compute rank RMSE per iteration per scenario ----------------------------
rmse_results <- list()

for (scenario_name in names(all_results)) {
  df <- all_results[[scenario_name]]
  iterations <- unique(df$iteration)

  iter_rmse <- list()
  for (iter in iterations) {
    iter_df <- df[df$iteration == iter, ]

    row <- data.frame(
      scenario  = scenario_name,
      iteration = iter,
      stringsAsFactors = FALSE
    )

    for (metric in METRICS) {
      if (metric %in% names(iter_df)) {
        row[[metric]] <- rank_rmse(iter_df$niche_breadth, iter_df[[metric]])
      } else {
        row[[metric]] <- NA_real_
      }
    }

    iter_rmse[[length(iter_rmse) + 1]] <- row
  }

  rmse_results[[scenario_name]] <- bind_rows(iter_rmse)
}

rmse_all <- bind_rows(rmse_results)

# --- Summarise: mean rank RMSE per scenario ----------------------------------
rmse_summary <- rmse_all %>%
  group_by(scenario) %>%
  summarise(across(all_of(METRICS), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  arrange(scenario)

cat("\n=== Mean Rank RMSE per Scenario ===\n\n")
print(as.data.frame(rmse_summary), digits = 2)

# --- Summarise: overall mean rank RMSE per metric ----------------------------
overall_rmse <- rmse_all %>%
  summarise(across(all_of(METRICS), ~ mean(.x, na.rm = TRUE)))

cat("\n=== Overall Mean Rank RMSE (across all scenarios) ===\n\n")
print(as.data.frame(overall_rmse), digits = 2)

# --- Save outputs ------------------------------------------------------------
write.csv(rmse_summary, "results/rank_rmse_summary.csv", row.names = FALSE)
cat("\nSaved: results/rank_rmse_summary.csv\n")

# --- Comparison table: Spearman rho alongside rank RMSE ----------------------
# Load existing Spearman correlations if available
spearman_file <- "results/oracle_correlations_summary.csv"
if (file.exists(spearman_file)) {
  spearman_df <- read.csv(spearman_file)

  cat("\n=== Side-by-side: Mean |Spearman rho| vs Mean Rank RMSE ===\n\n")

  # Reshape for comparison
  rmse_long <- rmse_summary %>%
    pivot_longer(cols = all_of(METRICS), names_to = "metric",
                 values_to = "rank_rmse")

  spearman_long <- spearman_df %>%
    pivot_longer(cols = any_of(METRICS), names_to = "metric",
                 values_to = "spearman_rho")

  # Overall means
  rmse_overall <- rmse_long %>%
    group_by(metric) %>%
    summarise(mean_rank_rmse = mean(rank_rmse, na.rm = TRUE), .groups = "drop")

  spearman_overall <- spearman_long %>%
    group_by(metric) %>%
    summarise(mean_abs_rho = mean(abs(spearman_rho), na.rm = TRUE), .groups = "drop")

  comparison <- left_join(spearman_overall, rmse_overall, by = "metric") %>%
    arrange(mean_rank_rmse)

  print(as.data.frame(comparison), digits = 3)
  write.csv(comparison, "results/spearman_vs_rmse_comparison.csv", row.names = FALSE)
  cat("\nSaved: results/spearman_vs_rmse_comparison.csv\n")
}
