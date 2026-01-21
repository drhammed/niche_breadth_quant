#' =============================================================================
#' Unified Simulation Functions for Species Distribution
#' =============================================================================
#'
#' All Gaussian response functions for simulating species distributions.
#'
#' =============================================================================

# Source helper functions
# source("R/helper_functions.R")

#' Generate Abundances - Symmetric Gaussian (Single Species)
#'
#' Simulates abundances for a single species across multiple environmental
#' gradients using a symmetric Gaussian response function.
#'
#' @param x Matrix of environmental values (sites x n_env)
#' @param optima Vector of species optima (one per environment)
#' @param breadth Vector of niche breadths (one per environment)
#' @param h Vector of peak height parameters (one per environment)
#' @return List with Abundance.oracle, Abundance.sim, and niche.breadth.Oracle
generate_abundances_symmetric <- function(x, optima, breadth, h) {
  n.sites <- nrow(x)
  n.x <- ncol(x)

  Abundance.oracle <- matrix(0, n.sites, n.x)
  niche.breadth.Oracle <- numeric(n.x)

  for (i in 1:n.x) {
    # Symmetric Gaussian response
    Abundance.oracle[, i] <- rpois(n.sites,
      30 * h[i] * sqrt(exp(-(x[, i] - optima[i])^2 / (2 * breadth[i]^2))))

    # Calculate weighted variance for oracle niche breadth
    if (sum(Abundance.oracle[, i]) == 0) {
      niche.breadth.Oracle[i] <- 0
    } else {
      w.mean <- sum(Abundance.oracle[, i] * x[, i]) / sum(Abundance.oracle[, i])
      w.ss <- sum(Abundance.oracle[, i] * (x[, i] - w.mean)^2)
      w.var <- w.ss / sum(Abundance.oracle[, i])
      niche.breadth.Oracle[i] <- sqrt(w.var)
    }
  }

  # Average across environments
  Abundance.oracle.avg <- apply(Abundance.oracle, 1, sum) / n.x
  Abundance.sim <- rpois(length(Abundance.oracle.avg), Abundance.oracle.avg)

  return(list(
    Abundance.oracle = Abundance.oracle.avg,
    Abundance.sim = Abundance.sim,
    niche.breadth.Oracle = niche.breadth.Oracle
  ))
}


#' Generate Abundances - Asymmetric Gaussian (Single Species)
#'
#' Simulates abundances for a single species across multiple environmental
#' gradients using an asymmetric Gaussian response function.
#'
#' @param x Matrix of environmental values (sites x n_env)
#' @param optima Vector of species optima (one per environment)
#' @param peak Scalar peak abundance value
#' @param left_breadth Vector of left-side breadths (one per environment)
#' @param right_breadth Vector of right-side breadths (one per environment)
#' @param strength.x Vector of predictor strengths (one per environment)
#' @return List with Abundance.oracle, Abundance.sim, and niche.breadth.Oracle
generate_abundances_asymmetric <- function(x, optima, peak, left_breadth,
                                            right_breadth, strength.x) {
  n.sites <- NROW(x)
  n.x <- NCOL(x)

  left_matrix <- matrix(0, n.sites, n.x)
  right_matrix <- matrix(0, n.sites, n.x)
  Abundance.oracle <- matrix(0, n.sites, n.x)
  niche.breadth.Oracle <- 0

  for (i in 1:n.x) {
    # Asymmetric Gaussian - left side
    left_matrix[, i] <- ifelse(x[, i] < optima[i],
      exp((log(peak) - (1/strength.x[i]) * ((x[, i] - optima[i]) / left_breadth[i])^2)), 0)

    # Asymmetric Gaussian - right side
    right_matrix[, i] <- ifelse(x[, i] >= optima[i],
      exp((log(peak) - (1/strength.x[i]) * ((x[, i] - optima[i]) / right_breadth[i])^2)), 0)

    Abundance.oracle[, i] <- left_matrix[, i] + right_matrix[, i]
    niche.breadth.Oracle <- niche.breadth.Oracle +
      sqrt(weighted.variance(as.matrix(x[, i]), Abundance.oracle[, i]))
  }

  Abundance.oracle.avg <- apply(Abundance.oracle, 1, sum) / n.x
  Abundance.sim <- rpois(length(Abundance.oracle.avg), Abundance.oracle.avg)
  breadth.Oracle <- niche.breadth.Oracle / n.x

  return(list(
    Abundance.oracle = Abundance.oracle.avg,
    Abundance.sim = Abundance.sim,
    niche.breadth.Oracle = niche.breadth.Oracle
  ))
}


#' Generate Breadth Parameters
#'
#' Generates niche breadth parameter values from specified distribution.
#' Uses scenario-specific parameters.
#'
#' @param n Total number of values to generate
#' @param distribution One of "uniform", "normal", or "gamma"
#' @param n_env Number of environmental variables (1 or 2)
#' @param response_type Either "symmetric" or "asymmetric"
#' @param params List of distribution parameters (optional, overrides defaults)
#' @return Vector of breadth values
generate_breadth_params <- function(n, distribution = "uniform", n_env = 1,
                                     response_type = "symmetric", params = NULL) {

  # Parameters for each scenario

  if (distribution == "uniform") {
    # Uniform: same for all scenarios [0.1, 1.2]
    breadth <- runif(n, min = BREADTH_UNIFORM_MIN, max = BREADTH_UNIFORM_MAX)

  } else if (distribution == "normal") {
    # Normal distribution - parameters vary by scenario
    if (response_type == "symmetric" && n_env == 1) {
      # Symmetric 1-env: mean=0.02, sd=0.8, floor=0.08
      b <- rnorm(n, mean = BREADTH_NORMAL_MEAN_SYM_1ENV, sd = BREADTH_NORMAL_SD)
      breadth <- pmax(b, BREADTH_NORMAL_FLOOR_SYM_1ENV)
    } else if (response_type == "symmetric" && n_env == 2) {
      # Symmetric 2-env: mean=0.01, sd=0.8, NO floor
      breadth <- rnorm(n, mean = BREADTH_NORMAL_MEAN_2ENV, sd = BREADTH_NORMAL_SD)
    } else if (response_type == "asymmetric") {
      # Asymmetric (1-env and 2-env): mean=0.01, sd=0.8, NO floor
      breadth <- rnorm(n, mean = BREADTH_NORMAL_MEAN_ASY_1ENV, sd = BREADTH_NORMAL_SD)
    }

  } else if (distribution == "gamma") {
    # Gamma distribution - parameters vary by scenario
    b <- rgamma(n, shape = BREADTH_GAMMA_SHAPE, rate = BREADTH_GAMMA_RATE)

    if (response_type == "symmetric" && n_env == 1) {
      # Symmetric 1-env: pmax floor=0.08
      breadth <- pmax(b, BREADTH_GAMMA_FLOOR_SYM_1ENV)
    } else if (response_type == "symmetric" && n_env == 2) {
      # Symmetric 2-env: pmax floor=0.05
      breadth <- pmax(b, BREADTH_GAMMA_FLOOR_SYM_2ENV)
    } else if (response_type == "asymmetric" && n_env == 1) {
      # Asymmetric 1-env: use set_max_min transformation [0.05, 5]
      breadth <- set_max_min(b, new_max = BREADTH_GAMMA_ASY_1ENV_MAX,
                             new_min = BREADTH_GAMMA_ASY_1ENV_MIN)
    } else if (response_type == "asymmetric" && n_env == 2) {
      # Asymmetric 2-env: pmax floor=0.05
      breadth <- pmax(b, BREADTH_GAMMA_FLOOR_ASY_2ENV)
    }

  } else {
    stop("Unknown distribution: ", distribution)
  }

  return(breadth)
}


#' Generate Community - Unified Function
#'
#' Main function to generate simulated community data. Handles all combinations
#' of symmetric/asymmetric response and uniform/normal/gamma breadth distributions.
#'
#' @param n_species Number of species to simulate
#' @param n_sites Number of sites (communities)
#' @param n_env Number of environmental variables (1 or 2)
#' @param response_type Either "symmetric" or "asymmetric"
#' @param breadth_distribution Either "uniform", "normal", or "gamma"
#' @param breadth_params Optional list of distribution parameters
#' @return List containing:
#'   - community_abundance_matrix: sites x species abundance matrix
#'   - breadth.Oracle: true niche breadths per species
#'   - x: environmental gradient matrix
#'   - Abundance.oracle: expected abundances
generate_community <- function(n_species, n_sites, n_env,
                                response_type = "symmetric",
                                breadth_distribution = "uniform",
                                breadth_params = NULL) {

  community_abundance_matrix <- matrix(0, n_sites, n_species)
  Abundance.oracle <- matrix(0, n_sites, n_species)
  x <- NULL

  # Repeat until valid community generated (all species/sites have occurrences)
  repeat {
    # Generate environmental gradients
    x <- matrix(rnorm(n_env * n_sites), n_sites, n_env)
    colnames(x) <- paste0("env.", 1:n_env)

    # Generate species optima
    optima <- matrix(runif(n_env * n_species, -2, 2), n_env, n_species)

    # Generate response based on type
    if (response_type == "symmetric") {
      # Symmetric Gaussian parameters
      h <- matrix(runif(n_env * n_species, min = 0.3, max = 1), n_env, n_species)
      breadth <- matrix(
        generate_breadth_params(n_env * n_species, breadth_distribution, n_env,
                                response_type = "symmetric", params = breadth_params),
        n_env, n_species
      )

      breadth.Oracle <- numeric(n_species)
      for (species in 1:n_species) {
        sim <- generate_abundances_symmetric(x, optima[, species],
                                              breadth[, species], h[, species])
        community_abundance_matrix[, species] <- sim$Abundance.sim
        breadth.Oracle[species] <- mean(sim$niche.breadth.Oracle)
        Abundance.oracle[, species] <- sim$Abundance.oracle
      }

    } else if (response_type == "asymmetric") {
      # Asymmetric Gaussian parameters
      peak <- round(runif(n_species, 10, 500))
      strength.x <- matrix(1, n_env, n_species)  # Fixed at 1

      left_breadth <- matrix(
        generate_breadth_params(n_env * n_species, breadth_distribution, n_env,
                                response_type = "asymmetric", params = breadth_params),
        n_env, n_species
      )
      right_breadth <- matrix(
        generate_breadth_params(n_env * n_species, breadth_distribution, n_env,
                                response_type = "asymmetric", params = breadth_params),
        n_env, n_species
      )

      breadth.Oracle <- numeric(n_species)
      for (species in 1:n_species) {
        sim <- generate_abundances_asymmetric(x, optima[, species], peak[species],
                                               left_breadth[, species],
                                               right_breadth[, species],
                                               strength.x[, species])
        community_abundance_matrix[, species] <- sim$Abundance.sim
        breadth.Oracle[species] <- sim$niche.breadth.Oracle
        Abundance.oracle[, species] <- sim$Abundance.oracle
      }

    } else {
      stop("Unknown response_type: ", response_type)
    }

    # Check validity: all species and sites have non-zero occurrences
    n_species_valid <- sum(colSums(community_abundance_matrix) != 0)
    n_sites_valid <- sum(rowSums(community_abundance_matrix) != 0)

    if ((n_species_valid == n_species) && (n_sites_valid == n_sites)) {
      break
    }
  }

  # Name columns
  colnames(community_abundance_matrix) <- paste0("sp", 1:n_species)
  names(breadth.Oracle) <- paste0("sp", 1:n_species)

  return(list(
    community_abundance_matrix = community_abundance_matrix,
    breadth.Oracle = breadth.Oracle,
    x = x,
    Abundance.oracle = Abundance.oracle
  ))
}


#' Run Single Scenario Simulation
#'
#' Wrapper function to run a complete simulation for one scenario.
#'
#' @param scenario_row Row from SCENARIOS data frame
#' @param n_species Number of species
#' @param n_sites Number of sites
#' @param seed Random seed for this scenario
#' @return List with community data and parameters
run_scenario_simulation <- function(scenario_row, n_species, n_sites, seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  result <- generate_community(
    n_species = n_species,
    n_sites = n_sites,
    n_env = scenario_row$n_env,
    response_type = scenario_row$response_type,
    breadth_distribution = scenario_row$breadth_distribution
  )

  result$scenario <- scenario_row
  return(result)
}
