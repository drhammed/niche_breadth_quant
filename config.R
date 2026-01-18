#' =============================================================================
#' Configuration file for Chapter 1- Niche Breadth Quantification 
#' =============================================================================
#'
#' This file contains all configurable parameters for the simulation study.
#' Modify these values to change simulation settings.
#' =============================================================================

# ------------------------------------------------------------------------------
# Simulation Parameters
# ------------------------------------------------------------------------------

N_SPECIES <- 30        # Number of species per simulation
N_SITES   <- 500       # Number of sites (communities)
N_REPS    <- 30        # Number of simulation iterations per scenario
RANDOM_SEED <- 9999    # Random seed for reproducibility

# ------------------------------------------------------------------------------
# Parallel Processing
# ------------------------------------------------------------------------------

NUM_CORES <- 6    # Number of CPU cores for parallel processing
                  # Our lab can use up to 32 I believe, but I'm using 6 here based on my computer 
                  # Adjust based on your machine's resources.

# ------------------------------------------------------------------------------
# Co-occurrence Metric Parameters
# ------------------------------------------------------------------------------

CO_OCCUR_REPS    <- 100  # Number of random subsamples for co-occurrence metrics
CO_OCCUR_PSAMPLE <- 4    # Size of random subsample (plots)
CO_OCCUR_PSAMPLE2 <- 2   # Size of secondary subsample for Simpson index

# ------------------------------------------------------------------------------
# Hypervolume Parameters
# ------------------------------------------------------------------------------

NICHE_ROVER_NSAMPLES <- 1000  # Number of samples for nicheROVER
HV_SAMPLES_PER_POINT <- 10   # Samples per point for Blonder hypervolume

# ------------------------------------------------------------------------------
# Latent Variable Model Parameters
# ------------------------------------------------------------------------------

NLV <- 5  # Number of latent variables for ecoCopula

# ------------------------------------------------------------------------------
# Breadth Distribution Parameters
# ------------------------------------------------------------------------------

# Uniform distribution bounds
BREADTH_UNIFORM_MIN <- 0.1
BREADTH_UNIFORM_MAX <- 1.2

# Normal distribution parameters
BREADTH_NORMAL_MEAN_1ENV <- 0.02
BREADTH_NORMAL_MEAN_2ENV <- 0.01
BREADTH_NORMAL_SD <- 0.8

# Gamma distribution parameters
BREADTH_GAMMA_SHAPE <- 2
BREADTH_GAMMA_RATE  <- 2

# Gamma floor values (pmax threshold)

BREADTH_GAMMA_FLOOR_1ENV <- 0.08
BREADTH_GAMMA_FLOOR_2ENV <- 0.05


BREADTH_NORMAL_FLOOR <- 0.08

# ------------------------------------------------------------------------------
# Scenario Selection
# ------------------------------------------------------------------------------
#
# Specify which scenarios to run by setting filters below.
# Set to NULL to include all options for that parameter.
#
# Examples:
#   RUN_RESPONSE_TYPES <- c("symmetric")           # Only symmetric
#   RUN_RESPONSE_TYPES <- c("symmetric", "asymmetric")  # Both (default)
#   RUN_BREADTH_DISTRIBUTIONS <- c("normal")       # Only normal
#   RUN_N_ENVS <- c(2)                             # Only 2-environment
#   RUN_N_ENVS <- NULL                             # All environments (default)
# ------------------------------------------------------------------------------

RUN_RESPONSE_TYPES <- NULL   # NULL = all; or c("symmetric"), c("asymmetric"), c("symmetric", "asymmetric")
RUN_BREADTH_DISTRIBUTIONS <- NULL # NULL = all; or c("uniform"), c("normal"), c("gamma"), or combinations
RUN_N_ENVS <- NULL                # NULL = all; or c(1), c(2), c(1, 2)

# ------------------------------------------------------------------------------
# Scenario Definitions
# ------------------------------------------------------------------------------

# Define all 12 possible scenarios
ALL_SCENARIOS <- data.frame(
  scenario_id = 1:12,
  n_env = rep(c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2), 1),
  response_type = rep(c("symmetric", "symmetric", "symmetric",
                        "asymmetric", "asymmetric", "asymmetric"), 2),
  breadth_distribution = rep(c("uniform", "normal", "gamma"), 4),
  stringsAsFactors = FALSE
)

# Add descriptive names
ALL_SCENARIOS$name <- paste0(
  ifelse(ALL_SCENARIOS$response_type == "symmetric", "Sym_", "Asy_"),
  ALL_SCENARIOS$breadth_distribution, "_env",
  ALL_SCENARIOS$n_env
)

# Apply filters to select scenarios
SCENARIOS <- ALL_SCENARIOS
if (!is.null(RUN_RESPONSE_TYPES)) {
  SCENARIOS <- SCENARIOS[SCENARIOS$response_type %in% RUN_RESPONSE_TYPES, ]
}
if (!is.null(RUN_BREADTH_DISTRIBUTIONS)) {
  SCENARIOS <- SCENARIOS[SCENARIOS$breadth_distribution %in% RUN_BREADTH_DISTRIBUTIONS, ]
}
if (!is.null(RUN_N_ENVS)) {
  SCENARIOS <- SCENARIOS[SCENARIOS$n_env %in% RUN_N_ENVS, ]
}

# Reset row names
rownames(SCENARIOS) <- NULL

# ------------------------------------------------------------------------------
# Output File Names
# ------------------------------------------------------------------------------

# Function to generate output filename for each scenario
get_output_filename <- function(scenario_name, type = "rds") {
  paste0("result_df_", scenario_name, ".", type)
}

# ------------------------------------------------------------------------------
# Empirical Data Settings
# ------------------------------------------------------------------------------

# BBS species data file
SPECIES_DATA_FILE <- "data/spp_mod.csv"

# WorldClim bioclim data directory
BIOCLIM_DIR <- "data/bioclim"

# Minimum species occurrences for empirical analysis
MIN_SPECIES_OCCURRENCES <- 500

# Number of species to select for empirical validation
N_EMPIRICAL_SPECIES <- 30

# Year range for BBS data subsetting
BBS_YEAR_MIN <- 2015

# PCA components to retain for climate variables
N_PCA_COMPONENTS <- 4

# Species list for empirical analysis
# These are the 30 species used in this paper
# If set to NULL, it'll instead randomly sample species
EMPIRICAL_SPECIES <- c(
  "Contopus virens",
  "Dolichonyx oryzivorus",
  "Parkesia noveboracensis",
  "Seiurus aurocapilla",
  "Aphelocoma woodhouseii",
  "Patagioenas fasciata",
  "Setophaga striata",
  "Cyanocitta cristata",
  "Piranga rubra",
  "Coragyps atratus",
  "Passerina caerulea",
  "Tachycineta bicolor",
  "Empidonax oberholseri",
  "Corvus brachyrhynchos / ossifragus",
  "Piranga olivacea",
  "Vireo olivaceus",
  "Melanerpes carolinus",
  "Setophaga petechia",
  "Nucifraga columbiana",
  "Pheucticus melanocephalus",
  "Rhynchophanes mccownii",
  "Loxia curvirostra",
  "Carduelis flammea / hornemanni",
  "Sayornis saya",
  "Vermivora cyanoptera",
  "Corvus ossifragus",
  "Loxia leucoptera",
  "Woodpecker sp.",
  "Spinus tristis",
  "Tympanuchus phasianellus"
)

# ------------------------------------------------------------------------------
# Print Configuration Summary
# ------------------------------------------------------------------------------

print_config <- function() {
  cat("\n==================== SIMULATION CONFIG ====================\n")
  cat(sprintf("Species: %d | Sites: %d | Iterations: %d\n",
              N_SPECIES, N_SITES, N_REPS))
  cat(sprintf("Parallel cores: %d | Random seed: %d\n",
              NUM_CORES, RANDOM_SEED))
  cat(sprintf("\nScenarios to run: %d of %d\n", nrow(SCENARIOS), nrow(ALL_SCENARIOS)))
  cat("------------------------------------------------------------------\n")
  for (i in 1:nrow(SCENARIOS)) {
    cat(sprintf("  [%2d] %s\n", SCENARIOS$scenario_id[i], SCENARIOS$name[i]))
  }
  cat("==================================================================\n\n")
}
