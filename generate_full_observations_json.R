# Full 900-observation JSON for all simulation scenarios
# Each scenario has 30 iterations × 30 species = 900 rows
# Columns: iteration, sci.name, niche_breadth (oracle), + 9 niche breadth metrics

library(jsonlite)

sim_dir <- "results/simulations"

# List all simulation RDS files
rds_files <- list.files(sim_dir, pattern = "\\.rds$", full.names = TRUE)

# Read each scenario into a named list
full_observations <- list()

for (f in rds_files) {
  scenario_name <- tools::file_path_sans_ext(basename(f))
  df <- readRDS(f)

  # Keep relevant columns: iteration, sci.name, niche_breadth, and the 9 metrics
  metric_cols <- c("SimpSSI", "beta.a", "beta.w", "om_tol",
                    "nr_hv", "hv_blond", "nb_Gam", "nb_latent", "nb_dist")
  keep_cols <- intersect(c("iteration", "sci.name", "niche_breadth", metric_cols), colnames(df))
  df <- df[, keep_cols, drop = FALSE]

  # Round numeric columns for readability
  num_cols <- sapply(df, is.numeric)
  df[num_cols] <- lapply(df[num_cols], round, digits = 4)

  full_observations[[scenario_name]] <- as.list(df)
}

# Sort by scenario name for consistency
full_observations <- full_observations[sort(names(full_observations))]

# Write to JSON
output_path <- "results/full_simulation_mat.json"
write_json(full_observations, output_path, pretty = TRUE, auto_unbox = TRUE)

cat("Written", length(full_observations), "scenarios to", output_path, "\n")
cat("Each scenario contains", nrow(readRDS(rds_files[1])), "observations\n")
