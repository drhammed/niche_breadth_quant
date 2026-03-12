#' =============================================================================
#' Niche Breadth Metrics
#' =============================================================================
#'
#' This file contains all 11 niche breadth estimation methods:
#'   1. SimpSSI (Multi-site Simpson)
#'   2. Beta.a (Additive beta partitioning)
#'   3. Beta.w (Whittaker's beta)
#'   4. om_tol (OMI tolerance)
#'   5. nr_hv (nicheROVER hypervolume)
#'   6. hv_blond (Blonder hypervolume)
#'   7. nb_Gam (GAM-based niche breadth, general version)
#'   8. nb_latent (Latent variable model, covariance-based)
#'   9. nb_dist (Weighted average environmental distance)
#'  10. LPCD_mean (Alpha diversity contribution mean)
#'  11. LPCD_sd (Alpha diversity contribution SD)
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
#' @param nsamples Number of posterior samples
#' @param niche.breadth Optional oracle niche breadth for comparison
#' @param empirical Logical; if TRUE, uses empirical data settings
#' @return Data frame with sci.name and nr_hypervolume
nr_hypervolume_fun <- function(data, env_vars, nsamples = 10,
                                niche.breadth = NULL, empirical = FALSE) {

  # Melt data to long format
  df_melted <- cbind(as.data.frame(data), as.data.frame(env_vars))
  species_cols <- colnames(data)
  env_cols <- colnames(env_vars)
  n_env <- length(env_cols)

  df_melted <- tidyr::gather(df_melted, key = "species", value = "PA",
                              tidyr::all_of(species_cols))
  df_melted <- df_melted[, c("species", "PA", env_cols)]
  df_melted <- df_melted[df_melted$PA != 0, ]
  df_melted <- df_melted[, !(colnames(df_melted) %in% "PA")]

  # For empirical data: scale environmental variables
  if (empirical) {
    df_melted[env_cols] <- scale(df_melted[env_cols])
  }

  # Define dimensionality and initialize psi
  # set.seed(999)
  d <- n_env
  psi <- crossprod(matrix(rnorm(d^2), d, d))

  # Generate posterior parameters based on data type
  if (empirical) {
    # Empirical data: Custom niw.post with specific params
    custom_niw_post_empirical <- function(nsamples, X, Psi) {
      par <- nicheROVER::niw.coeffs(X, lambda = rnorm(ncol(X)), kappa = 20,
                                     Psi = Psi, nu = 5)
      # Add regularization
      par$Psi <- par$Psi + diag(ncol(par$Psi)) * 1e-2
      nicheROVER::rniw(nsamples, par$lambda, par$kappa, par$Psi, par$nu)
    }

    sim.par <- tapply(1:nrow(df_melted), df_melted$species, function(ii) {
      X <- df_melted[ii, env_cols, drop = FALSE]
      custom_niw_post_empirical(nsamples = nsamples, X = X, Psi = psi)
    })
  } else if (n_env == 1) {
    # Simulation single environment: Custom niw.post with regularization
    custom_niw_post <- function(nsamples, X, Psi) {
      par <- nicheROVER::niw.coeffs(X)
      # Add regularization to avoid singularity issues
      par$Psi <- par$Psi + diag(ncol(par$Psi)) * 1e-2
      nicheROVER::rniw(nsamples, par$lambda, par$kappa, par$Psi, par$nu)
    }

    sim.par <- tapply(1:nrow(df_melted), df_melted$species, function(ii) {
      X <- df_melted[ii, env_cols, drop = FALSE]
      custom_niw_post(nsamples = nsamples, X = X, Psi = psi)
    })
  } else {
    # Simulation two or more environments: Default nicheROVER::niw.post
    sim.par <- tapply(1:nrow(df_melted), df_melted$species, function(ii) {
      X <- df_melted[ii, env_cols, drop = FALSE]
      nicheROVER::niw.post(nsamples = nsamples, X = X)
    })
  }

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
#' Implementation varies by data type:
#' - Simulation: samples.per.point = 10, default bandwidth
#' - Empirical: Scaled data, Silverman bandwidth, dynamic samples.per.point
#'
#' @param data Presence-absence matrix (sites x species)
#' @param env_vars Environmental variables (sites x n_env)
#' @param niche.breadth Optional oracle niche breadth for comparison
#' @param empirical Logical; if TRUE, uses empirical data settings
#' @return Data frame with sci.name and hypervolume
hypervolume_blond_fun <- function(data, env_vars, niche.breadth = NULL,
                                   empirical = FALSE) {

  # Melt data
  df_melted <- cbind(as.data.frame(data), as.data.frame(env_vars))
  species_cols <- colnames(data)
  env_cols <- colnames(env_vars)

  df_melted <- tidyr::gather(df_melted, key = "species", value = "PA",
                              tidyr::all_of(species_cols))
  df_melted <- df_melted[, c("species", "PA", env_cols)]
  df_melted <- df_melted[df_melted$PA != 0, ]
  df_melted <- df_melted[, !(colnames(df_melted) %in% "PA")]

  # For empirical data: scale environmental variables
  if (empirical) {
    df_melted[env_cols] <- scale(df_melted[env_cols])
  }

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
      if (empirical) {
        # Empirical mode: Silverman bandwidth, dynamic samples.per.point, weights
        bandwidth <- hypervolume::estimate_bandwidth(sim_data, method = "silverman")
        bandwidth[is.na(bandwidth) | bandwidth == 0] <- 0.01
        samples_per_point <- ceiling((10^(3 + sqrt(ncol(sim_data)))) / nrow(sim_data))

        hv <- hypervolume::hypervolume_gaussian(
          sim_data,
          name = sp,
          weight = rep(1 / nrow(sim_data), nrow(sim_data)),
          samples.per.point = samples_per_point,
          kde.bandwidth = bandwidth
        )
      } else {
        # Simulation mode: fixed samples.per.point = 10
        hv <- hypervolume::hypervolume_gaussian(sim_data, name = sp,
                                                 samples.per.point = 10)
      }
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
# 7. GAM-BASED NICHE BREADTH (nb_Gam) - General version
# =============================================================================

#' Weighted Covariance Matrix
#'
#' Computes weighted covariance matrix for calculating niche breadth.
#'
#' @param X Matrix of environmental variables
#' @param w Vector of weights (e.g., predicted probabilities)
#' @return Weighted covariance matrix
weighted_cov <- function(X, w) {
  X <- as.matrix(X)
  w <- as.numeric(w)

  ok <- complete.cases(X) & is.finite(w)
  X <- X[ok, , drop = FALSE]
  w <- w[ok]

  if (nrow(X) == 0 || sum(w) <= 0) {
    return(matrix(NA, ncol = ncol(X), nrow = ncol(X)))
  }

  w <- w / sum(w)
  mu <- colSums(X * w)
  X_centered <- sweep(X, 2, mu, "-")
  cov_w <- t(X_centered) %*% (X_centered * w)

  return(cov_w)
}

#' Summarize Weighted Covariance into Scalar Breadth
#'
#' @param cov_mat Weighted covariance matrix
#' @param summary_method One of "mean_sd", "rms_sd", "trace", "determinant"
#' @return Scalar niche breadth value
summarize_breadth_from_cov <- function(cov_mat,
                                       summary_method = c("mean_sd",
                                                          "rms_sd",
                                                          "trace",
                                                          "determinant")) {
  summary_method <- match.arg(summary_method)

  if (any(is.na(cov_mat))) return(NA_real_)

  eig <- try(eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values,
             silent = TRUE)
  if (inherits(eig, "try-error")) return(NA_real_)

  eig[eig < 0 & abs(eig) < 1e-10] <- 0
  if (any(eig < 0)) return(NA_real_)

  sds <- sqrt(diag(cov_mat))

  out <- switch(summary_method,
                mean_sd     = mean(sds),
                rms_sd      = sqrt(mean(diag(cov_mat))),
                trace       = sum(diag(cov_mat)),
                determinant = {
                  det_val <- det(cov_mat)
                  if (!is.finite(det_val) || det_val < 0) NA_real_ else det_val^(1 / (2 * ncol(cov_mat)))
                }
  )

  return(out)
}

#' Build GAM Formula
#'
#' Constructs a GAM formula supporting independent, joint, or hybrid models.
#'
#' @param response_name Name of the response variable
#' @param predictors Character vector of predictor names
#' @param model_type One of "independent", "joint", "hybrid"
#' @param joint_engine One of "te", "ti"
#' @param pairwise_interactions Optional list of variable pairs for hybrid model
#' @param k Optional basis dimension
#' @return A formula object
build_gam_formula <- function(response_name = "response",
                              predictors,
                              model_type = c("independent", "joint", "hybrid"),
                              joint_engine = c("te", "ti"),
                              pairwise_interactions = NULL,
                              k = NULL) {
  model_type <- match.arg(model_type)
  joint_engine <- match.arg(joint_engine)

  if (length(predictors) == 0) stop("No predictors supplied.")

  smooth_1d <- function(v) {
    if (is.null(k)) {
      paste0("s(", v, ")")
    } else {
      paste0("s(", v, ", k = ", k, ")")
    }
  }

  smooth_joint <- function(vars) {
    engine <- joint_engine
    vars_txt <- paste(vars, collapse = ", ")
    if (is.null(k)) {
      paste0(engine, "(", vars_txt, ")")
    } else if (length(k) == 1) {
      paste0(engine, "(", vars_txt, ", k = ", k, ")")
    } else {
      k_txt <- paste(k, collapse = ", ")
      paste0(engine, "(", vars_txt, ", k = c(", k_txt, "))")
    }
  }

  n_env <- length(predictors)

  if (n_env == 1) {
    rhs <- smooth_1d(predictors)
  } else if (model_type == "independent") {
    rhs <- paste(vapply(predictors, smooth_1d, character(1)), collapse = " + ")
  } else if (model_type == "joint") {
    rhs <- smooth_joint(predictors)
  } else if (model_type == "hybrid") {
    main_terms <- vapply(predictors, smooth_1d, character(1))

    if (is.null(pairwise_interactions)) {
      pairwise_interactions <- utils::combn(predictors, 2, simplify = FALSE)
    } else {
      pairwise_interactions <- lapply(pairwise_interactions, as.character)
    }

    pair_terms <- vapply(pairwise_interactions, smooth_joint, character(1))
    rhs <- paste(c(main_terms, pair_terms), collapse = " + ")
  }

  stats::as.formula(paste(response_name, "~", rhs))
}

#' Fit GAM Model (General)
#'
#' Fit a GAM with error handling, supporting multiple model types.
#'
#' @param data.df Data frame with response and predictors
#' @param predictors Predictor column names
#' @param model_type One of "independent", "joint", "hybrid"
#' @param joint_engine One of "te", "ti"
#' @param pairwise_interactions Optional list of variable pairs
#' @param family GLM family (default: binomial logit)
#' @param method Fitting method (default: "REML")
#' @param k Optional basis dimension
#' @return List with gam_model, warning, and formula
fit_gam_general <- function(data.df,
                            predictors,
                            model_type = c("independent", "joint", "hybrid"),
                            joint_engine = c("te", "ti"),
                            pairwise_interactions = NULL,
                            family = stats::binomial(link = "logit"),
                            method = "REML",
                            k = NULL) {
  model_type <- match.arg(model_type)
  joint_engine <- match.arg(joint_engine)

  formula_obj <- build_gam_formula(
    response_name = "response",
    predictors = predictors,
    model_type = model_type,
    joint_engine = joint_engine,
    pairwise_interactions = pairwise_interactions,
    k = k
  )

  result <- NULL
  warningMessage <- NULL

  tryCatch({
    gam_model <- withCallingHandlers({
      mgcv::gam(formula_obj, data = data.df, family = family, method = method)
    }, warning = function(w) {
      warningMessage <<- w$message
      invokeRestart("muffleWarning")
    })

    if (inherits(gam_model, "gam") && !isTRUE(gam_model$converged)) {
      warningMessage <- "Convergence issues"
    }

    result <- list(
      gam_model = gam_model,
      warning = warningMessage,
      formula = formula_obj
    )
  }, error = function(e) {
    result <- list(
      gam_model = NULL,
      warning = paste("Error:", conditionMessage(e)),
      formula = formula_obj
    )
  })

  return(result)
}

#' Calculate GAM-based Niche Breadth (General)
#'
#' Uses Generalized Additive Models to estimate niche breadth based on
#' weighted covariance of environmental variables.
#'
#' @param distribution_data Presence-absence matrix (sites x species)
#' @param env_variables Environmental variables (sites x n_env)
#' @param summary_method One of "mean_sd", "rms_sd", "trace", "determinant"
#' @param model_type One of "independent", "joint", "hybrid"
#' @param joint_engine One of "te", "ti"
#' @param pairwise_interactions Optional list of variable pairs
#' @param family GLM family (default: binomial logit)
#' @param method Fitting method (default: "REML")
#' @param k Optional basis dimension
#' @param verbose Logical; print progress?
#' @return Data frame with Niche_Breadth, Gam_Warning, and Formula per species
estimate_nicheBreadth_Gam_general <- function(distribution_data,
                                              env_variables,
                                              summary_method = c("mean_sd",
                                                                 "rms_sd",
                                                                 "trace",
                                                                 "determinant"),
                                              model_type = c("independent",
                                                             "joint",
                                                             "hybrid"),
                                              joint_engine = c("te", "ti"),
                                              pairwise_interactions = NULL,
                                              family = stats::binomial(link = "logit"),
                                              method = "REML",
                                              k = NULL,
                                              verbose = TRUE) {
  summary_method <- match.arg(summary_method)
  model_type <- match.arg(model_type)
  joint_engine <- match.arg(joint_engine)

  Env_variables <- as.data.frame(env_variables)
  colnames(Env_variables) <- paste0("Env", seq_len(ncol(Env_variables)))

  distribution_data <- as.data.frame(distribution_data)

  n.species <- ncol(distribution_data)
  n.sites <- nrow(distribution_data)

  if (nrow(Env_variables) != n.sites) {
    stop("distribution_data and env_variables must have the same number of rows.")
  }

  gam_Warning <- character(n.species)
  niche.breadth <- rep(NA_real_, n.species)
  gam_formula <- character(n.species)
  fitted_values <- matrix(NA_real_, nrow = n.sites, ncol = n.species)
  colnames(fitted_values) <- colnames(distribution_data)

  for (i in seq_len(n.species)) {
    data.df <- data.frame(
      response = distribution_data[[i]],
      Env_variables,
      check.names = FALSE
    )

    gam_result <- fit_gam_general(
      data.df = data.df,
      predictors = colnames(Env_variables),
      model_type = model_type,
      joint_engine = joint_engine,
      pairwise_interactions = pairwise_interactions,
      family = family,
      method = method,
      k = k
    )

    gam_formula[i] <- deparse(gam_result$formula)

    gam_Warning[i] <- if (!is.null(gam_result$warning)) {
      gam_result$warning
    } else {
      "noWarning"
    }

    if (!is.null(gam_result$gam_model)) {
      pred <- try(stats::predict(gam_result$gam_model, type = "response"),
                  silent = TRUE)

      if (!inherits(pred, "try-error")) {
        fitted_values[, i] <- pred
        cov_w <- weighted_cov(Env_variables, pred)
        niche.breadth[i] <- summarize_breadth_from_cov(
          cov_mat = cov_w,
          summary_method = summary_method
        )

        if (verbose) {
          rsq <- try(summary(gam_result$gam_model)$r.sq, silent = TRUE)
          if (!inherits(rsq, "try-error")) {
            print(c(species = i, r_sq = rsq, breadth = niche.breadth[i]))
          } else {
            print(c(species = i, breadth = niche.breadth[i]))
          }
        }
      }
    }
  }

  results <- data.frame(
    Niche_Breadth = niche.breadth,
    Gam_Warning = gam_Warning,
    Formula = gam_formula,
    stringsAsFactors = FALSE,
    row.names = colnames(distribution_data)
  )

  attr(results, "fitted_values") <- fitted_values

  return(results)
}


# =============================================================================
# 8. LATENT VARIABLE MODEL (nb_latent) - Covariance-based
# =============================================================================

#' Calculate Latent Variable-based Niche Breadth
#'
#' Uses ecoCopula package to fit latent variable models and estimate niche
#' breadth from the covariance of latent scores at occupied sites.
#'
#' @param distribution_data Presence-absence matrix (sites x species)
#' @param nlv Number of latent variables (default: 3)
#' @param verbose Logical; print progress?
#' @return Named vector of niche breadths
estimate_nicheBreadth_Latents <- function(distribution_data, nlv = 3, verbose = TRUE) {
  n.species <- ncol(distribution_data)

  # Ensure species names exist
  spp_names <- colnames(distribution_data)
  if (is.null(spp_names)) {
    spp_names <- paste0("sp", seq_len(n.species))
    colnames(distribution_data) <- spp_names
  }

  # Filter columns to avoid fitting issues in ecoCopula
  distribution_data.smaller <- eliminate_cols(
    distribution_data,
    min(colSums(distribution_data)),
    max(colSums(distribution_data))
  )

  # Fit stacked SDM and extract latent variables
  eco_PA <- ecoCopula::stackedsdm(
    distribution_data.smaller,
    ~ 1,
    data = distribution_data,
    family = "binomial"
  )

  eco_lvs <- ecoCopula::cord(eco_PA, nlv = nlv)
  eco_lvs <- as.data.frame(matrix(eco_lvs$scores, nrow(distribution_data), nlv))
  colnames(eco_lvs) <- paste0("LV", seq_len(nlv))

  niche.breadth <- rep(NA_real_, n.species)
  names(niche.breadth) <- spp_names

  for (i in seq_len(n.species)) {
    spp_name <- spp_names[i]
    spp_occ <- distribution_data[, i] > 0

    LV_i <- eco_lvs[spp_occ, , drop = FALSE]

    if (nrow(LV_i) < 2) {
      niche.breadth[i] <- NA
      if (verbose) message(spp_name, ": fewer than 2 occupied sites; returning NA")
      next
    }

    if (ncol(LV_i) == 1) {
      niche.breadth[i] <- sqrt(stats::var(LV_i[, 1], na.rm = TRUE))
    } else {
      S <- stats::cov(LV_i, use = "complete.obs")
      niche.breadth[i] <- sqrt(det(S + diag(1e-8, ncol(S))))
    }

    if (verbose) {
      message(spp_name, ": breadth = ", round(niche.breadth[i], 5))
    }
  }

  return(niche.breadth)
}


# =============================================================================
# 9. WEIGHTED AVERAGE ENVIRONMENTAL DISTANCE (nb_dist)
# =============================================================================

#' Calculate Weighted Average Environmental Distance Niche Breadth
#'
#' Computes niche breadth as the weighted average pairwise environmental
#' distance among occupied sites, where weights are abundance products.
#'
#' @param distribution_data Presence-absence or abundance matrix (sites x species)
#' @param env_variables Environmental variables (sites x n_env)
#' @return Named vector of niche breadths
estimate_nicheBreadth_avgDist_weighted <- function(distribution_data, env_variables) {
  if (is.null(env_variables)) {
    stop("env_variables must be provided.")
  }

  distribution_data <- as.matrix(distribution_data)
  env_variables <- as.matrix(env_variables)

  if (nrow(distribution_data) != nrow(env_variables)) {
    stop("distribution_data and env_variables must have the same number of rows.")
  }

  n_species <- ncol(distribution_data)
  niche_breadth <- rep(NA_real_, n_species)
  names(niche_breadth) <- colnames(distribution_data)

  env_dist <- as.matrix(dist(env_variables))

  for (i in seq_len(n_species)) {
    w <- distribution_data[, i]
    occ <- w > 0

    if (sum(occ) < 2) {
      niche_breadth[i] <- NA_real_
      next
    }

    w_sub <- w[occ]
    d_sub <- env_dist[occ, occ, drop = FALSE]

    W <- outer(w_sub, w_sub)
    niche_breadth[i] <- sum(d_sub[lower.tri(d_sub)] * W[lower.tri(W)]) / sum(W[lower.tri(W)])
  }

  niche_breadth
}


# =============================================================================
# 10. ALPHA DIVERSITY CONTRIBUTION MEAN (LPCD_mean)
# =============================================================================

#' Calculate Alpha Diversity Contribution Mean Niche Breadth
#'
#' Estimates niche breadth based on species' mean contributions to local
#' alpha diversity across sites. Based on occurrence data only.
#'
#' @param abund_mat Abundance matrix (sites x species)
#' @param na.as.zero Logical; treat NAs as zeros? (default: TRUE)
#' @return Named vector of niche breadths
estimate_nicheBreadth_from_alphaMean <- function(abund_mat, na.as.zero = TRUE) {

  n.sites <- nrow(abund_mat)
  n.species <- ncol(abund_mat)

  Sum.sites <- rowSums(abund_mat)
  total_sum <- sum(abund_mat)

  LPCD.alpha <- matrix(0, n.sites, n.species)

  for (i in 1:n.sites) {
    for (j in 1:n.species) {
      n_ij <- abund_mat[i, j]
      n_i  <- Sum.sites[i]

      if (n_ij > 0) {
        p_ij <- n_ij / n_i
        LPCD.alpha[i, j] <- (n_ij / total_sum) * log(1 + p_ij)
      } else {
        LPCD.alpha[i, j] <- if (na.as.zero) 0 else NA_real_
      }
    }
  }

  niche.breadth <- apply(LPCD.alpha, 2, mean, na.rm = TRUE)
  names(niche.breadth) <- colnames(abund_mat)

  return(niche.breadth)
}


# =============================================================================
# 11. ALPHA DIVERSITY CONTRIBUTION SD (LPCD_sd)
# =============================================================================

#' Calculate Alpha Diversity Contribution SD Niche Breadth
#'
#' Estimates niche breadth as the inverse of the standard deviation of
#' species' contributions to local alpha diversity. Based on occurrence only.
#'
#' @param abund_mat Abundance matrix (sites x species)
#' @param na.as.zero Logical; treat NAs as zeros? (default: TRUE)
#' @return Named vector of niche breadths
estimate_nicheBreadth_from_alphaSD <- function(abund_mat, na.as.zero = TRUE) {

  n.sites <- nrow(abund_mat)
  n.species <- ncol(abund_mat)

  Sum.sites <- rowSums(abund_mat)
  total_sum <- sum(abund_mat)

  LPCD.alpha <- matrix(0, n.sites, n.species)

  for (i in 1:n.sites) {
    for (j in 1:n.species) {
      n_ij <- abund_mat[i, j]
      n_i  <- Sum.sites[i]

      if (n_ij > 0) {
        p_ij <- n_ij / n_i
        LPCD.alpha[i, j] <- (n_ij / total_sum) * log(1 + p_ij)
      } else {
        LPCD.alpha[i, j] <- if (na.as.zero) 0 else NA_real_
      }
    }
  }

  alpha_sd <- apply(LPCD.alpha, 2, sd, na.rm = TRUE)
  niche.breadth <- 1 / (alpha_sd + 1e-8)
  names(niche.breadth) <- colnames(abund_mat)

  return(niche.breadth)
}


# =============================================================================
# WRAPPER FUNCTION: CALCULATE ALL METRICS
# =============================================================================

#' Calculate All Niche Breadth Metrics
#'
#' Wrapper function that calculates all 11 niche breadth metrics for a given
#' community dataset.
#'
#' @param sim.com Presence-absence matrix (sites x species)
#' @param env_vars Environmental variables matrix (sites x n_env)
#' @param niche.breadth Optional oracle niche breadth (named vector)
#' @param co_occur_params List of parameters for co_occur function
#' @param nr_params List of parameters for nicheROVER function
#' @param gam_params List of parameters for GAM function
#' @param nlv Number of latent variables for ecoCopula
#' @param verbose Logical; print progress?
#' @param empirical Logical; if TRUE, uses empirical data settings for hypervolume
#'   functions (scaling, different bandwidth/niw params). Default FALSE for simulations.
#' @return Data frame with all metrics per species
calculate_all_metrics <- function(sim.com, env_vars, niche.breadth = NULL,
                                   co_occur_params = list(reps = 100,
                                                          psample = 4,
                                                          psample2 = 2),
                                   nr_params = list(nsamples = 10),
                                   gam_params = list(summary_method = "determinant",
                                                     model_type = "hybrid",
                                                     joint_engine = "ti"),
                                   nlv = 5, verbose = TRUE,
                                   empirical = FALSE) {

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
                                   niche.breadth = niche.breadth,
                                   empirical = empirical)

  if (verbose) message("Calculating Blonder hypervolume...")
  hv_result <- hypervolume_blond_fun(sim.com, env_vars,
                                      niche.breadth = niche.breadth,
                                      empirical = empirical)

  if (verbose) message("Calculating GAM-based niche breadth...")
  gam_result <- estimate_nicheBreadth_Gam_general(
    sim.com, env_vars,
    summary_method = gam_params$summary_method,
    model_type = gam_params$model_type,
    joint_engine = gam_params$joint_engine,
    verbose = verbose
  )

  if (verbose) message("Calculating latent variable niche breadth...")
  latent_result <- estimate_nicheBreadth_Latents(sim.com,
                                                  nlv = nlv, verbose = verbose)

  if (verbose) message("Calculating weighted env distance niche breadth...")
  dist_result <- estimate_nicheBreadth_avgDist_weighted(sim.com, env_vars)

  if (verbose) message("Calculating LPCD mean niche breadth...")
  lpcd_mean_result <- estimate_nicheBreadth_from_alphaMean(sim.com)

  if (verbose) message("Calculating LPCD SD niche breadth...")
  lpcd_sd_result <- estimate_nicheBreadth_from_alphaSD(sim.com)

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
    LPCD_mean = lpcd_mean_result[match(co_occ_result$sci.name, names(lpcd_mean_result))],
    LPCD_sd = lpcd_sd_result[match(co_occ_result$sci.name, names(lpcd_sd_result))],
    stringsAsFactors = FALSE
  )

  return(result_df)
}
