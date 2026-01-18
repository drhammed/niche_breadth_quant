#' =============================================================================
#' Co-occurrence Based Theta Functions
#' =============================================================================
#'
#' Set of functions for estimating species niche breadth based on compositional
#' data using co-occurrence based theta metric introduced by Fridley et al. (2007).
#'
#' Author: David Zeleny (zeleny.david@gmail.com)
#' Based on codes by Jason Fridley (Fridley et al. 2007) and David Zeleny (Zeleny 2009)
#'
#' References:
#' - Fridley et al. (2007) J. Ecology 95: 707-722
#' - Zeleny (2009) J. Ecology 97: 10-17
#' - Botta-Dukat (2012) J. Vegetation Science 23: 201-207
#' =============================================================================

#' Calculate Theta Metric
#'
#' Main function for calculating co-occurrence based theta metric of species
#' habitat specialization.
#'
#' @param input.matrix Community data matrix (samples x species). If not
#'   presence-absence, will be automatically transformed.
#' @param thresh Minimal frequency of species (default: 5)
#' @param psample Size of random subsample (default: 5)
#' @param reps Number of random subsamples (default: 10)
#' @param method Beta-diversity algorithm: "additive", "multiplicative",
#'   "pairwise.jaccard", "pairwise.sorensen", "pairwise.simpson",
#'   "multi.sorensen", "multi.simpson", "beals", or "beta.div"
#' @param q Generalization parameter for multiplicative beta (default: 0)
#' @param rarefaction Logical; use rarefaction for multiplicative method?
#' @param parallel Logical; use parallel computation?
#' @param no.cores Number of cores for parallel (default: 2)
#' @param remove.out Logical; remove outliers?
#' @param verbal Logical; show progress?
#' @return Data frame with sci.name, local.avgS, occur.freq, meanco, theta, theta.sd
#' @export
calculate.theta <- function(input.matrix, species.data = NULL, thresh = 5,
                             psample = 5, reps = 10, method = "multiplicative",
                             q = 0, rarefaction = TRUE,
                             beta.div.method = 'hellinger',
                             beta.div.sqrt.D = FALSE, beta.div.samp = TRUE,
                             beals.file = NULL, pa.transform = FALSE,
                             force.subsample = FALSE, parallel = FALSE,
                             no.cores = 2, remove.out = FALSE,
                             out.metric = 'sorensen', verbal = FALSE,
                             juicer = FALSE, tcltk = FALSE) {

  METHODS <- c('additive', 'multiplicative', 'pairwise.jaccard',
               'pairwise.sorensen', 'pairwise.simpson', 'multi.sorensen',
               'multi.simpson', 'rao', 'beals', 'beta.div')
  method.n <- pmatch(method, METHODS)
  if (is.na(method.n)) stop('invalid method')
  if (method.n == -1) stop('ambiguous method')
  method <- METHODS[method.n]

  win.pb <- NULL

  if (is.na(reps) || reps < 2)
    stop("Number of random subsamples must be integer >= 2")
  if (is.na(thresh) || thresh < 2)
    stop("Minimum frequency of species must be integer >= 2")
  if (thresh < psample)
    stop("Minimum frequency of species must be >= size of random subsamples")

  if (!is.matrix(input.matrix)) input.matrix <- as.matrix(input.matrix)
  if (is.null(row.names(input.matrix)))
    row.names(input.matrix) <- seq(1, nrow(input.matrix))

  if (method %in% c('additive', 'multiplicative', 'pairwise.jaccard',
                    'pairwise.sorensen', 'pairwise.simpson',
                    'multi.sorensen', 'multi.simpson', 'beals'))
    pa.transform <- TRUE

  if (pa.transform) input.matrix <- ifelse(input.matrix > 0, 1, 0)

  # Species selection
  Nplots <- nrow(input.matrix)
  plots.per.spp <- colSums(input.matrix > 0)
  select.spp <- plots.per.spp[plots.per.spp >= thresh]
  Nspp <- length(select.spp)

  # For beals method - transform to beals smoothed form
  if (method == "beals") {
    if (is.null(beals.file)) {
      beals.matrix <- beals.2(input.matrix, include = TRUE, verbal = verbal)
      for (co in seq(1, ncol(input.matrix))) {
        if (sum(input.matrix[, co]) > 0) {
          beals.temp <- beals.matrix[, co][as.logical(input.matrix[, co])]
          stats.temp <- fivenum(beals.temp)
          iqr <- diff(stats.temp[c(2, 4)])
          beals.thresh <- min(beals.temp[!(beals.temp < stats.temp[2] - 1.5 * iqr)])
          beals.matrix[, co] <- as.numeric(beals.matrix[, co] >= beals.thresh)
        } else {
          beals.matrix[, co] <- 0
        }
      }
    } else {
      beals.matrix <- as.matrix(read.delim(file = beals.file, row.names = 1,
                                           check.names = FALSE))
    }
  }

  if (!parallel) {
    temp.res <- lapply(1:Nspp, FUN = function(sp) {
      if (method == 'beals') {
        temp.matrix <- beals.matrix[input.matrix[, colnames(input.matrix) ==
                                                   names(select.spp[sp])] > 0, ]
      } else {
        temp.matrix <- input.matrix[input.matrix[, colnames(input.matrix) ==
                                                   names(select.spp[sp])] > 0, ]
      }
      temp.matrix <- temp.matrix[, colSums(temp.matrix) > 0]
      sci.name <- labels(select.spp[sp])
      calculate.theta.0(temp.matrix = temp.matrix, sci.name = sci.name, sp = sp,
                        remove.out = remove.out, out.metric = out.metric,
                        thresh = thresh, psample = psample, reps = reps,
                        method = method, q = q, rarefaction = rarefaction,
                        beta.div.method = beta.div.method,
                        beta.div.sqrt.D = beta.div.sqrt.D,
                        beta.div.samp = beta.div.samp,
                        force.subsample = force.subsample, parallel = parallel,
                        win.pb = win.pb, verbal = verbal, juicer = juicer)
    })
  }

  if (parallel) {
    workers <- parallel::makeCluster(no.cores)
    parallel::clusterExport(workers, c('calculate.theta.0', 'input.matrix',
                                       'select.spp', 'remove.out', 'thresh',
                                       'psample', 'reps', 'method', 'parallel'),
                           envir = environment())
    temp.res <- parallel::parLapply(workers, 1:Nspp, fun = function(sp) {
      if (method == 'beals') {
        temp.matrix <- beals.matrix[input.matrix[, colnames(input.matrix) ==
                                                   names(select.spp[sp])] > 0, ]
      } else {
        temp.matrix <- input.matrix[input.matrix[, colnames(input.matrix) ==
                                                   names(select.spp[sp])] > 0, ]
      }
      temp.matrix <- temp.matrix[, colSums(temp.matrix) > 0]
      sci.name <- labels(select.spp[sp])
      calculate.theta.0(temp.matrix = temp.matrix, sci.name = sci.name, sp = sp,
                        remove.out = remove.out, out.metric = out.metric,
                        thresh = thresh, psample = psample, reps = reps,
                        method = method, q = q, rarefaction = rarefaction,
                        beta.div.method = beta.div.method,
                        beta.div.sqrt.D = beta.div.sqrt.D,
                        beta.div.samp = beta.div.samp,
                        force.subsample = force.subsample, parallel = parallel,
                        win.pb = NULL, verbal = verbal, juicer = juicer)
    })
    parallel::stopCluster(workers)
  }

  theta.out <- do.call(rbind.data.frame, temp.res)
  rownames(theta.out) <- NULL
  names(theta.out) <- c('sci.name', 'local.avgS', 'occur.freq', 'meanco',
                        'theta', 'theta.sd')
  theta.out$sci.name <- as.character(theta.out$sci.name)

  if (!is.null(species.data)) {
    theta.out <- as.data.frame(cbind(sci.name = theta.out[, 1],
                                     species.data[as.character(theta.out[, 'sci.name']), 1:2],
                                     theta.out[, -1]))
  }

  return(theta.out)
}


#' Internal Theta Calculation Function
#'
#' @keywords internal
calculate.theta.0 <- function(temp.matrix, sci.name, sp, remove.out, out.metric,
                               thresh, psample, reps, method, rarefaction, q,
                               beta.div.method, beta.div.sqrt.D, beta.div.samp,
                               force.subsample, parallel, win.pb, verbal, juicer) {

  # Outlier removal
  if (remove.out) {
    if (out.metric == 'sorensen')
      veg.dist <- as.matrix(vegan::vegdist(temp.matrix > 0))
    if (out.metric == 'binary.euclidean')
      veg.dist <- as.matrix(dist(temp.matrix > 0))
    if (out.metric == 'euclidean')
      veg.dist <- as.matrix(dist(temp.matrix))
    diag(veg.dist) <- NA
    distances <- rowMeans(veg.dist, na.rm = TRUE)
    outliers <- distances > (mean(distances) + 2 * sd(distances))
    temp.matrix <- temp.matrix[!outliers, ]
    temp.matrix <- temp.matrix[, colSums(temp.matrix) > 0]
  }

  # Subsampling methods
  if (method %in% c('additive', 'multiplicative', 'multi.sorensen',
                    'multi.simpson', 'beals', 'rao') &
      !(method %in% c('multiplicative', 'beals') & rarefaction) ||
      (method %in% c('pairwise.jaccard', 'pairwise.sorensen',
                     'pairwise.simpson', 'beta.div') & force.subsample)) {

    if (!nrow(temp.matrix) < thresh) {
      rn.temp.matrix <- matrix(rownames(temp.matrix), ncol = reps,
                               nrow = nrow(temp.matrix), byrow = FALSE)
      sample.temp.matrix <- apply(rn.temp.matrix, 2,
                                  FUN = function(x) sample(x, psample))

      mc.mat <- array(0, dim = c(psample, ncol(temp.matrix), reps))
      for (i in 1:reps) mc.mat[, , i] <- temp.matrix[sample.temp.matrix[, i], ]
      total.rich <- colSums(apply(mc.mat, c(2, 3), sum) > 0)
      mean.alpha <- colMeans(apply(mc.mat > 0, c(1, 3), sum))

      if (method == "multiplicative") {
        if (q == 0) {
          Wbeta.vec <- total.rich / mean.alpha
        } else {
          Wbeta.vec <- unlist(lapply(1:reps, FUN = function(i)
            vegetarian::d(mc.mat[, , i], lev = 'beta', q = q)))
        }
      }
      if (method == "beals") Wbeta.vec <- total.rich / mean.alpha
      if (method == "additive") Wbeta.vec <- total.rich - mean.alpha
      if (method == "pairwise.jaccard")
        Wbeta.vec <- unlist(lapply(1:reps, FUN = function(i)
          mean(betapart::beta.pair(mc.mat[, , i], index = 'jaccard')$beta.jac)))
      if (method == "pairwise.sorensen")
        Wbeta.vec <- unlist(lapply(1:reps, FUN = function(i)
          mean(betapart::beta.pair(mc.mat[, , i], index = 'sorensen')$beta.sor)))
      if (method == "pairwise.simpson")
        Wbeta.vec <- unlist(lapply(1:reps, FUN = function(i)
          mean(betapart::beta.pair(mc.mat[, , i], index = 'sorensen')$beta.sim)))
      if (method == "multi.sorensen")
        Wbeta.vec <- unlist(lapply(1:reps, FUN = function(i)
          betapart::beta.multi(mc.mat[, , i], index = 'sorensen')$beta.SOR))
      if (method == "multi.simpson")
        Wbeta.vec <- unlist(lapply(1:reps, FUN = function(i)
          betapart::beta.multi(mc.mat[, , i], index = 'sorensen')$beta.SIM))
      if (method == "rao")
        Wbeta.vec <- unlist(lapply(1:reps, FUN = function(i)
          cati::RaoRel(t(mc.mat[, , i]), dfunc = NULL, dphyl = NULL,
                       Jost = TRUE)$TD$Beta_prop))
      if (method == "beta.div")
        Wbeta.vec <- unlist(lapply(1:reps, FUN = function(i)
          adespatial::beta.div(mc.mat[, , i], method = beta.div.method,
                               sqrt.D = beta.div.sqrt.D, nperm = 0)$SStotal_BDtotal[2]))

      theta <- mean(Wbeta.vec)
      theta.sd <- sd(Wbeta.vec)
      meanco <- mean(total.rich)
      local.avgS <- mean(mean.alpha)
      occur.freq <- nrow(temp.matrix)

      result <- list(sci.name, local.avgS, occur.freq, meanco, theta, theta.sd)
      return(result)
    }
  }

  # Non-subsampling methods
  if (method %in% c('pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson',
                    'beta.div') & !force.subsample) {
    if (!nrow(temp.matrix) < thresh) {
      total.rich <- sum(colSums(temp.matrix) > 0)
      mean.alpha <- mean(rowSums(temp.matrix > 0))

      if (method == "pairwise.jaccard")
        theta <- mean(betapart::beta.pair(temp.matrix, index = 'jaccard')$beta.jac)
      if (method == "pairwise.sorensen")
        theta <- mean(betapart::beta.pair(temp.matrix, index = 'sorensen')$beta.sor)
      if (method == "pairwise.simpson")
        theta <- mean(betapart::beta.pair(temp.matrix, index = 'sorensen')$beta.sim)
      if (method == "beta.div")
        theta <- adespatial::beta.div(temp.matrix, method = beta.div.method,
                                      sqrt.D = beta.div.sqrt.D,
                                      nperm = 0)$SStotal_BDtotal[2]

      meanco <- total.rich
      local.avgS <- mean.alpha
      occur.freq <- nrow(temp.matrix)

      result <- list(sci.name, local.avgS, occur.freq, meanco, theta, theta.sd = 0)
      return(result)
    }
  }

  # Rarefaction method for multiplicative/beals
  if (method %in% c('multiplicative', 'beals') & rarefaction) {
    if (!nrow(temp.matrix) < thresh) {
      theta <- beta.raref(comm = temp.matrix, sites = psample)
      total.rich <- sum(colSums(temp.matrix) > 0)
      meanco <- total.rich
      local.avgS <- theta$alpha
      occur.freq <- nrow(temp.matrix)
      result <- list(sci.name, local.avgS, occur.freq, meanco, theta$beta, theta$beta.sd)
      return(result)
    }
  }
}


#' Beals Smoothing Function
#'
#' Modified from vegan package for calculating species pool using Beals method.
#'
#' @param x Input compositional matrix
#' @param include Logical; include target species?
#' @param verbal Logical; show progress?
#' @return Beals smoothed matrix
#' @keywords internal
beals.2 <- function(x, include = TRUE, verbal = FALSE) {
  x <- as.matrix(x)
  x[x > 0] <- 1
  refX <- x
  incSp <- include
  refX <- as.matrix(refX)
  M <- crossprod(refX, refX)
  C <- diag(M)
  M <- sweep(M, 2, replace(C, C == 0, 1), "/")
  if (!incSp) for (i in 1:ncol(refX)) M[i, i] <- 0
  S <- rowSums(x)
  b <- x
  for (i in 1:nrow(x)) {
    b[i, ] <- rowSums(sweep(M, 2, x[i, ], "*"))
  }
  SM <- rep(S, ncol(x))
  if (!incSp) SM <- SM - x
  b <- b / replace(SM, SM == 0, 1)
  b
}


#' Beta Rarefaction
#'
#' Calculates rarefied beta diversity.
#'
#' @param comm Community matrix
#' @param sites Number of sites for rarefaction
#' @return List with alpha, beta, and beta.sd
#' @keywords internal
beta.raref <- function(comm, sites) {
  # Simple implementation - can be expanded
  n <- nrow(comm)
  if (n < sites) return(list(alpha = NA, beta = NA, beta.sd = NA))

  reps <- 100
  betas <- numeric(reps)

  for (i in 1:reps) {
    idx <- sample(1:n, sites)
    sub <- comm[idx, ]
    sub <- sub[, colSums(sub) > 0]
    gamma <- ncol(sub)
    alpha <- mean(rowSums(sub > 0))
    betas[i] <- gamma / alpha
  }

  list(
    alpha = mean(rowSums(comm > 0)),
    beta = mean(betas, na.rm = TRUE),
    beta.sd = sd(betas, na.rm = TRUE)
  )
}
