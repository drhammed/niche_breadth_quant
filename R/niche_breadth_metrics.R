#' =============================================================================
#' Niche Breadth Metrics
#' =============================================================================
#'
#' This file contains all 9 niche breadth estimation methods:
#'   1. SimpSSI (Multi-site Simpson)
#'   2. Beta.a (Additive beta partitioning)
#'   3. Beta.w (Whittaker's beta)
#'   4. om_tol (OMI tolerance)
#'   5. nr_hv (nicheROVER hypervolume)
#'   6. hv_blond (Blonder hypervolume)
#'   7. nb_Gam (GAM-based niche breadth)
#'   8. nb_latent (Latent variable model)
#'   9. nb_dist (Average Hellinger distance)
#'
#' Functions by Hammed Akande and Pedro Peres-Neto.
#' =============================================================================

# =============================================================================
# 1-3. CO-OCCURRENCE BASED METRICS (SimpSSI, Beta.a, Beta.w)
# =============================================================================

#' Calculate Co-occurrence Based Niche Breadth Metrics
#'
#' Computes multi-site Simpson (SimpSSI), additive beta (Beta.a), and
#' Whittaker's beta (Beta.w) for each species.
#'
#' @param sim.com Presence-absence matrix (sites x species)
#' @param n.species Number of species
#' @param reps Number of Monte Carlo repetitions
#' @param psample Number of plots for subset selection
#' @param psample2 Number of plots for Simpson calculation
#' @param species_names Optional vector of species names. If NULL (default),
#'   uses generic names "sp1", "sp2", etc. For empirical data, pass
#'   colnames(sim.com) to use real species names.
#' @return Data frame with sci.name, multi.sim, Beta.a, Beta.w
co_occur <- function(sim.com, n.species = NULL, reps = 100,
                      psample = 4, psample2 = 2, species_names = NULL) {

  if (is.null(n.species)) n.species <- ncol(sim.com)

  # Use provided species names or generate generic ones
  if (is.null(species_names)) {
    species_names <- paste0("sp", 1:n.species)
  }

  # Internal helper functions
  minbibj <- function(matrix) {
    nr <- dim(matrix)[1]
    sumbibj <- 0
    for (i in 1:(nr - 1)) {
      for (j in (i + 1):nr) {
        bi <- sum(matrix[i, ] & (!matrix[j, ]))
        bj <- sum(matrix[j, ] & (!matrix[i, ]))
        sumbibj <- sumbibj + min(bi, bj)
      }
    }
    sumbibj
  }

  zeroColumn <- function(matrix) {
    sum <- 0
    nc <- dim(matrix)[2]
    for (i in 1:nc) if (!sum(matrix[, i])) sum <- sum + 1
    sum
  }

  Simpson.multi <- function(x) {
    matrix <- as.matrix(x)
    sumSi <- sum(matrix)
    St <- ncol(matrix) - zeroColumn(matrix)
    a <- sumSi - St
    index <- a / (minbibj(matrix) + a)
    index
  }

  # Convert to 2-column format for generalist-specialist algorithm
  spp.vec <- NULL
  plot.vec <- NULL
  n_sites_sel <- nrow(sim.com)

  for (i in 1:n_sites_sel) {
    vec.true <- as.logical(sim.com[i, ])
    plot.vec <- c(plot.vec, rep(i, length = sum(sim.com[i, ])))
    spp.vec <- c(spp.vec, c(1:n.species)[vec.true])
  }

  out.simul <- data.frame(plot.vec, spp.vec)
  GOSmat <- out.simul

  # Species matrix - use provided species names or generic ones
  SppMat <- data.frame(
    sort(unique(GOSmat[, 2])),
    species_names[sort(unique(GOSmat[, 2]))]
  )

  plotID <- factor(GOSmat[, 1])
  SppID <- GOSmat[, 2]
  Nplots <- length(levels(plotID))
  richness <- tapply(SppID, plotID, length)
  max.rich <- max(richness)
  metacom <- table(plotID, SppID)

  # Select species occurring in >= psample plots
  plots.per.spp <- tapply(plotID, SppID, length)
  Species <- sort(unique(GOSmat[, 2]))[plots.per.spp >= psample]
  Nspp <- length(Species)

  # Initialize outputs
  Theta.out <- data.frame(
    sci.name = character(Nspp),
    multi.sim = numeric(Nspp),
    Beta.a = numeric(Nspp),
    Beta.w = numeric(Nspp),
    stringsAsFactors = FALSE
  )

  for (sp in 1:Nspp) {
    lab <- as.numeric(labels(metacom)[2][[1]])
    xlab <- c(1:dim(metacom)[2])
    metacol <- xlab[lab == Species[sp]]
    sp.plots <- as.logical(metacom[, metacol])
    sp.metacom <- metacom[sp.plots, ]
    Np <- dim(sp.metacom)[1]
    wide <- length(xlab)

    # Simpson calculation
    multi.sim.rep <- rep(0, reps)
    for (k in 1:reps) {
      rowselect <- sample(row.names(sp.metacom), psample2)
      trueRow <- is.element(row.names(sp.metacom), rowselect)
      rand.mat <- sp.metacom[trueRow, ]
      multi.sim.rep[k] <- 1 - (Simpson.multi(rand.mat))
    }
    multi.sim <- mean(multi.sim.rep)

    # Monte Carlo for beta metrics
    rpmat <- matrix(c(1:Np), reps, Np, byrow = TRUE)
    rpmat <- t(apply(rpmat, 1, function(x) sample(x, psample)))
    mc.mat <- array(0, dim = c(psample, wide, reps))
    for (i in 1:reps) {
      mc.mat[, , i] <- sp.metacom[rpmat[i, ], ]
    }

    colsum <- apply(mc.mat, c(2, 3), sum)
    colsum[colsum > 0] <- 1
    rich.vec <- colSums(colsum) - 1
    rich.vec2 <- colSums(colsum)
    mc.mat[mc.mat > 0] <- 1
    rowsum <- apply(mc.mat, c(1, 3), sum)
    Walpha.vec <- colMeans(rowsum)

    Beta.a.vec <- rich.vec - Walpha.vec
    Beta.w.vec <- rich.vec2 / Walpha.vec

    # Store results
    Theta.out$sci.name[sp] <- as.character(SppMat[, 2][SppMat[, 1] == Species[sp]])
    Theta.out$multi.sim[sp] <- multi.sim
    Theta.out$Beta.a[sp] <- mean(Beta.a.vec)
    Theta.out$Beta.w[sp] <- mean(Beta.w.vec)
  }

  return(Theta.out)
}


# =============================================================================
# 4. OMI TOLERANCE (om_tol)
# =============================================================================

#' Calculate OMI-based Niche Breadth
#'
#' Uses Outlying Mean Index (OMI) from ade4 package to estimate niche tolerance.
#'
#' @param env_vars Environmental variables (data frame or matrix)
#' @param PA Presence-absence matrix (sites x species)
#' @param nf Number of PCA axes (default: 2)
#' @param niche.breadth Optional oracle niche breadth for comparison
#' @return Data frame with sci.name and om_tol
omi_params_fun <- function(env_vars, PA, nf = 2, niche.breadth = NULL) {
  # Perform PCA on environmental variables
  pca.Env.virt <- ade4::dudi.pca(env_vars, center = TRUE, scale = TRUE,
                                  scannf = FALSE, nf = nf)

  # Calculate niche parameters using niche() and niche.param()
  omi <- ade4::niche(pca.Env.virt, as.data.frame(PA), scannf = FALSE)
  omi.param <- ade4::niche.param(omi)

  # Extract tolerance values
  om_tol <- as.data.frame(omi.param)$Tol

  # Build result data frame
  if (!is.null(niche.breadth)) {
    nb <- as.data.frame(niche.breadth)
    nb$sci.name <- row.names(nb)
    om_df <- as.data.frame(cbind(nb, om_tol))
    colnames(om_df) <- c("niche_b", "sci.name", "om_tol")
  } else {
    om_df <- data.frame(
      sci.name = rownames(omi.param),
      om_tol = om_tol,
      stringsAsFactors = FALSE
    )
  }

  return(om_df)
}


# =============================================================================
# 5. NICHE ROVER HYPERVOLUME (nr_hv)
# =============================================================================

#' Calculate nicheROVER Hypervolume
#'
#' Uses Bayesian approach (nicheROVER package) to estimate niche hypervolume.
#'
#' @param data Presence-absence matrix (sites x species)
#' @param env_vars Environmental variables (sites x n_env)
#' @param nsamples Number of posterior samples (default: 1000)
#' @param niche.breadth Optional oracle niche breadth for comparison
#' @return Data frame with sci.name and nr_hypervolume
nr_hypervolume_fun <- function(data, env_vars, nsamples = 1000,
                                niche.breadth = NULL) {

  # Melt data to long format
  df_melted <- cbind(as.data.frame(data), as.data.frame(env_vars))
  species_cols <- colnames(data)
  env_cols <- colnames(env_vars)

  df_melted <- tidyr::gather(df_melted, key = "species", value = "PA",
                              tidyr::all_of(species_cols))
  df_melted <- df_melted[, c("species", "PA", env_cols)]
  df_melted <- df_melted[df_melted$PA != 0, ]
  df_melted <- df_melted[, !(colnames(df_melted) %in% "PA")]

  # Custom niw.post function with regularization
  set.seed(999)
  niw.post.custom <- function(nsamples, X, Psi) {
    par <- nicheROVER::niw.coeffs(X, lambda = rnorm(ncol(as.matrix(X))),
                                  kappa = 20, Psi = Psi, nu = 5)
    # Add regularization for numerical stability
    par$Psi <- par$Psi + diag(ncol(par$Psi)) * 1e-2
    nicheROVER::rniw(nsamples, par$lambda, par$kappa, par$Psi, par$nu)
  }

  # Define dimensionality and initialize
  d <- length(env_cols)
  psi <- crossprod(matrix(rnorm(d^2), d, d))

  # Calculate niche parameters for each species
  sim.par <- tapply(1:nrow(df_melted), df_melted$species, function(ii) {
    X <- df_melted[ii, env_cols, drop = FALSE]
    if (ncol(X) == 1) {
      X <- as.matrix(X)
    }
    niw.post.custom(nsamples = nsamples, X = X, Psi = psi)
  })

  # Calculate posterior niche sizes
  sim.size <- sapply(sim.par, function(spec) {
    apply(spec$Sigma, 3, nicheROVER::niche.size, alpha = 0.95)
  })

  # Point estimates
  rover_nz <- rbind(est = colMeans(sim.size), se = apply(sim.size, 2, sd))

  # Build result
  nicheRover_hyp <- data.frame(
    sci.name = colnames(rover_nz),
    nr_hypervolume = rover_nz[1, ],
    stringsAsFactors = FALSE
  )

  # Add oracle
  if (!is.null(niche.breadth)) {
    nb <- as.data.frame(niche.breadth)
    nb$sci.name <- row.names(nb)
    nb_sim <- nb[nb$sci.name %in% nicheRover_hyp$sci.name, ]
    nicheRover_hyp <- merge(nb_sim, nicheRover_hyp, by = "sci.name", sort = FALSE)
    colnames(nicheRover_hyp)[2] <- "niche_b"
  }

  return(nicheRover_hyp)
}


# =============================================================================
# 6. BLONDER HYPERVOLUME (hv_blond)
# =============================================================================

#' Calculate Blonder Hypervolume
#'
#' Uses hypervolume package to estimate Gaussian kernel density hypervolumes.
#'
#' @param data Presence-absence matrix (sites x species)
#' @param env_vars Environmental variables (sites x n_env)
#' @param niche.breadth Optional oracle niche breadth for comparison
#' @return Data frame with sci.name and hypervolume
hypervolume_blond_fun <- function(data, env_vars, niche.breadth = NULL) {

  # Melt data
  df_melted <- cbind(as.data.frame(data), as.data.frame(env_vars))
  species_cols <- colnames(data)
  env_cols <- colnames(env_vars)

  df_melted <- tidyr::gather(df_melted, key = "species", value = "PA",
                              tidyr::all_of(species_cols))
  df_melted <- df_melted[, c("species", "PA", env_cols)]
  df_melted <- df_melted[df_melted$PA != 0, ]
  df_melted <- df_melted[, !(colnames(df_melted) %in% "PA")]

  # Calculate hypervolume for each species
  species_values <- unique(df_melted$species)
  hv_results <- data.frame(
    sci.name = character(),
    hypervolume = numeric(),
    stringsAsFactors = FALSE
  )

  for (sp in species_values) {
    sim_data <- df_melted[df_melted$species == sp, env_cols, drop = FALSE]
    tryCatch({
      hv <- hypervolume::hypervolume_gaussian(sim_data, name = sp,
                                               samples.per.point = 10)
      sp_hv <- hypervolume::get_volume(hv)
      hv_results <- rbind(hv_results,
                          data.frame(sci.name = sp, hypervolume = sp_hv))
    }, error = function(e) {
      warning(paste("Error for species", sp, ":", e$message))
      hv_results <- rbind(hv_results,
                          data.frame(sci.name = sp, hypervolume = NA))
    })
  }

  # Add oracle
  if (!is.null(niche.breadth)) {
    nb <- as.data.frame(niche.breadth)
    nb$sci.name <- row.names(nb)
    hv_results <- merge(nb, hv_results, by.x = "sci.name", by.y = "sci.name")
    colnames(hv_results)[2] <- "niche_b"
  }

  return(hv_results)
}


# =============================================================================
# 7. GAM-BASED NICHE BREADTH (nb_Gam)
# =============================================================================

#' Fit GAM Model
#'
#' Fit a GAM with error handling.
#'
#' @param data.df Data frame with response and predictors
#' @param predictors Predictor column names
#' @return List with gam_model and warning message
fit_gam <- function(data.df, predictors) {
  formula_str <- paste("response ~",
                       paste(paste("s(", predictors, ")", sep = ""),
                             collapse = " + "))
  formula_obj <- as.formula(formula_str)

  result <- NULL
  warningMessage <- NULL

  tryCatch({
    gam_model <- withCallingHandlers({
      mgcv::gam(formula_obj, data = data.df, family = binomial(link = "logit"))
    }, warning = function(w) {
      warningMessage <<- w$message
      invokeRestart("muffleWarning")
    })

    if (inherits(gam_model, "gam") && any(gam_model$converged != 0)) {
      warningMessage <- "Convergence Issues"
    }

    result <- list(gam_model = gam_model, warning = warningMessage)
  }, error = function(e) {
    result <- list(gam_model = NULL,
                   warning = paste("Error:", conditionMessage(e)))
  })

  return(result)
}


#' Calculate GAM-based Niche Breadth
#'
#' Uses Generalized Additive Models to estimate niche breadth based on
#' average distance among predicted values.
#'
#' @param distribution_data Presence-absence matrix (sites x species)
#' @param env_variables Environmental variables (sites x n_env)
#' @param verbose Logical; print progress?
#' @return Data frame with Niche_Breadth and Gam_Warning per species
estimate_nicheBreadth_Gam <- function(distribution_data, env_variables,
                                       verbose = TRUE) {
  Env_variables <- data.frame(env_variables)
  colnames(Env_variables) <- paste0("Env", 1:NCOL(env_variables))

  n.species <- ncol(distribution_data)
  n.sites <- nrow(distribution_data)
  gam_Warning <- character(n.species)
  predicted_values <- matrix(0, n.sites, n.species)

  for (i in 1:n.species) {
    data.df <- data.frame(response = distribution_data[, i], Env_variables)
    gam_result <- fit_gam(data.df, predictors = colnames(Env_variables))

    gam_Warning[i] <- if (!is.null(gam_result$warning)) {
      gam_result$warning
    } else {
      "noWarning"
    }

    if (!is.null(gam_result$gam_model)) {
      if (verbose) {
        print(c(i, summary(gam_result$gam_model)$r.sq))
      }
      predicted_values[, i] <- predict(gam_result$gam_model, type = "response")
    }
  }

  niche.breadth <- apply(as.matrix(dist(t(predicted_values))), 2, mean)
  names(niche.breadth) <- colnames(distribution_data)

  results <- data.frame(
    Niche_Breadth = niche.breadth,
    Gam_Warning = gam_Warning,
    stringsAsFactors = FALSE
  )

  return(results)
}


# =============================================================================
# 8. LATENT VARIABLE MODEL (nb_latent)
# =============================================================================

#' Calculate Latent Variable-based Niche Breadth
#'
#' Uses ecoCopula package to fit latent variable models and estimate niche
#' breadth from predicted probabilities.
#'
#' @param distribution_data Presence-absence matrix (sites x species)
#' @param env_variables Environmental variables (sites x n_env)
#' @param nlv Number of latent variables (default: 5)
#' @param verbose Logical; print progress?
#' @return Named vector of niche breadths
estimate_nicheBreadth_Latents <- function(distribution_data, env_variables,
                                           nlv = 5, verbose = TRUE) {
  n.species <- ncol(distribution_data)

  # Filter columns to avoid fitting issues
  distribution_data.smaller <- eliminate_cols(
    distribution_data,
    min(colSums(distribution_data)),
    max(colSums(distribution_data))
  )

  # Fit stacked SDM
  eco_PA <- ecoCopula::stackedsdm(distribution_data.smaller, ~ 1,
                                   data = distribution_data, family = "binomial")
  eco_lvs <- ecoCopula::cord(eco_PA, nlv = nlv)
  eco_lvs <- matrix(eco_lvs$scores, nrow(distribution_data.smaller), nlv)

  predicted_values <- matrix(0, nrow = nrow(distribution_data.smaller),
                             ncol = n.species)

  for (i in 1:n.species) {
    logistic_model <- mgcv::gam(distribution_data[, i] ~ eco_lvs,
                                family = binomial(link = "logit"))
    if (verbose) {
      print(c(i, summary(logistic_model)$r.sq))
    }
    predicted_values[, i] <- predict(logistic_model, type = "response")
  }

  niche.breadth <- apply(as.matrix(dist(t(predicted_values))), 2, mean)
  names(niche.breadth) <- colnames(distribution_data)

  return(niche.breadth)
}


# =============================================================================
# 9. AVERAGE HELLINGER DISTANCE (nb_dist)
# =============================================================================

#' Calculate Average Hellinger Distance Niche Breadth
#'
#' Simple metric based on average Hellinger distance among species.
#'
#' @param distribution_data Presence-absence or abundance matrix
#' @param env_variables Environmental variables (unused but kept for consistency)
#' @return Named vector of niche breadths
estimate_nicheBreadth_avg.Dist <- function(distribution_data, env_variables = NULL) {
  spe.Hellinger <- vegan::decostand(distribution_data, method = "hellinger")
  niche.breadth <- apply(as.matrix(dist(t(spe.Hellinger))), 2, mean)
  names(niche.breadth) <- colnames(distribution_data)

  return(niche.breadth)
}


# =============================================================================
# WRAPPER FUNCTION: CALCULATE ALL METRICS
# =============================================================================

#' Calculate All Niche Breadth Metrics
#'
#' Wrapper function that calculates all 9 niche breadth metrics for a given
#' community dataset.
#'
#' @param sim.com Presence-absence matrix (sites x species)
#' @param env_vars Environmental variables matrix (sites x n_env)
#' @param niche.breadth Optional oracle niche breadth (named vector)
#' @param co_occur_params List of parameters for co_occur function
#' @param nr_params List of parameters for nicheROVER function
#' @param nlv Number of latent variables for ecoCopula
#' @param verbose Logical; print progress?
#' @return Data frame with all metrics per species
calculate_all_metrics <- function(sim.com, env_vars, niche.breadth = NULL,
                                   co_occur_params = list(reps = 100,
                                                          psample = 4,
                                                          psample2 = 2),
                                   nr_params = list(nsamples = 1000),
                                   nlv = 5, verbose = TRUE) {

  n.species <- ncol(sim.com)

  # Ensure column names exist
  if (is.null(colnames(sim.com))) {
    colnames(sim.com) <- paste0("sp", 1:n.species)
  }

  # Name environmental columns
  env_vars <- as.data.frame(env_vars)
  if (is.null(colnames(env_vars))) {
    colnames(env_vars) <- paste0("env.", 1:ncol(env_vars))
  }

  if (verbose) message("Calculating co-occurrence metrics...")
  co_occ_result <- co_occur(sim.com, n.species,
                            reps = co_occur_params$reps,
                            psample = co_occur_params$psample,
                            psample2 = co_occur_params$psample2)

  if (verbose) message("Calculating OMI tolerance...")
  omi_result <- omi_params_fun(env_vars, sim.com, niche.breadth = niche.breadth)

  if (verbose) message("Calculating nicheROVER hypervolume...")
  nr_result <- nr_hypervolume_fun(sim.com, env_vars,
                                   nsamples = nr_params$nsamples,
                                   niche.breadth = niche.breadth)

  if (verbose) message("Calculating Blondel hypervolume...")
  hv_result <- hypervolume_blond_fun(sim.com, env_vars,
                                      niche.breadth = niche.breadth)

  if (verbose) message("Calculating GAM-based niche breadth...")
  gam_result <- estimate_nicheBreadth_Gam(sim.com, env_vars, verbose = verbose)

  if (verbose) message("Calculating latent variable niche breadth...")
  latent_result <- estimate_nicheBreadth_Latents(sim.com, env_vars,
                                                  nlv = nlv, verbose = verbose)

  if (verbose) message("Calculating average distance niche breadth...")
  dist_result <- estimate_nicheBreadth_avg.Dist(sim.com, env_vars)

  # Combine all results
  result_df <- data.frame(
    sci.name = co_occ_result$sci.name,
    SimpSSI = co_occ_result$multi.sim,
    beta.a = co_occ_result$Beta.a,
    beta.w = co_occ_result$Beta.w,
    om_tol = omi_result$om_tol[match(co_occ_result$sci.name, omi_result$sci.name)],
    nr_hv = nr_result$nr_hypervolume[match(co_occ_result$sci.name, nr_result$sci.name)],
    hv_blond = hv_result$hypervolume[match(co_occ_result$sci.name, hv_result$sci.name)],
    nb_Gam = gam_result$Niche_Breadth[match(co_occ_result$sci.name, rownames(gam_result))],
    nb_latent = latent_result[match(co_occ_result$sci.name, names(latent_result))],
    nb_dist = dist_result[match(co_occ_result$sci.name, names(dist_result))],
    stringsAsFactors = FALSE
  )

  return(result_df)
}
