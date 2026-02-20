#' Generate Fisher's log-series abundances with a dominant species
#'
#' @description
#' Creates a rank-abundance vector for \emph{S} species following Fisher's
#' log-series, while allocating a fixed fraction of individuals to a single
#' dominant species ("A"). The remaining individuals are distributed across
#' species B... in decreasing expectation \eqn{\propto \alpha x^{r} / r},
#' where \eqn{r = 2, \ldots, S}.
#'
#' @details
#' The algorithm:
#' \enumerate{
#'   \item Reserves \code{round(n_individuals * dominant_fraction)} individuals for
#'         species \code{"A"}.
#'   \item Computes relative expectations for ranks \code{2:n_species} using
#'         \eqn{\alpha x^{r} / r} and scales them to the remaining individuals.
#'   \item Rounds to integers and adjusts the second species (if present) so that
#'         the final sum equals \code{n_individuals}.
#'   \item Names the vector with \code{LETTERS[1:n_species]} and drops any
#'         species whose rounded abundance is \code{0}.
#' }
#'
#' @param n_species Integer (>= 1). Total number of species.
#' @param n_individuals Integer (>= 1). Total number of individuals to allocate.
#' @param dominant_fraction Numeric in \eqn{[0,1]}. Fraction of individuals assigned to
#'   species \code{"A"}; the remainder are split among species \code{"B"} ... according to
#'   the log-series.
#' @param alpha Numeric (> 0). Fisher's \eqn{\alpha} parameter controlling tail weight.
#' @param x Numeric, typically close to 1. The log-series scaling parameter.
#'
#' @return
#' A named numeric vector of integer abundances whose names are \code{"A"}, \code{"B"},
#' \code{"C"}, ... up to \code{n_species}. The vector sums to \code{n_individuals}.
#'
#' @note
#' Rounding may produce zeros for some tail species; these are omitted from the
#' return value. If you need a fixed-length vector including zeros, reindex
#' against \code{LETTERS[1:n_species]} after the call.
#'
#' @examples
#' abund <- generate_fisher_log_series(
#'   n_species = 10, n_individuals = 1000,
#'   dominant_fraction = 0.3, alpha = 3, x = 0.95
#' )
#' sum(abund) # 1000
#'
#' @seealso \code{\link{generate_heterogeneous_distribution}}
#' @export
generate_fisher_log_series <- function(
  n_species,
  n_individuals,
  dominant_fraction,
  alpha,
  x
) {
  stopifnot(n_species >= 1L, n_individuals >= 1L)

  if (n_species == 1L) {
    # With a single species, all individuals are assigned to A.
    out <- c(as.integer(n_individuals))
    names(out) <- "A"
    return(out)
  }

  n_dominant <- round(n_individuals * dominant_fraction)
  n_remaining <- n_individuals - n_dominant

  ranks <- 2:n_species
  rel <- alpha * (x^ranks) / ranks

  # a) Distribute n_remaining among non-dominants using largest remainder method
  abund_others <- integer(length(rel))
  if (n_remaining > 0) {
    if (sum(rel) > 0) {
      scaled <- rel / sum(rel) * n_remaining
      base <- floor(scaled)
      rem <- n_remaining - sum(base)

      if (rem > 0L) {
        frac <- scaled - base
        o <- order(frac, decreasing = TRUE)
        base[o[seq_len(rem)]] <- base[o[seq_len(rem)]] + 1L
      }
      abund_others <- as.integer(base)
    } else {
      # if rel is all zero, distribute uniformly
      abund_others <- rep(0L, length(rel))
      if (length(abund_others) > 0) {
        distrib <- as.integer(stats::rmultinom(1, size = n_remaining, prob = rep(1, length(rel)))[, 1])
        abund_others <- distrib
      }
    }
  }

  # b) Combine and ensure total sum is correct
  all_abund <- c(n_dominant, abund_others)

  # The sum should be correct now because `abund_others` sums to `n_remaining`
  # and `n_dominant + n_remaining` is `n_individuals`.
  # A final guard for floating point issues is still good practice.
  current_sum <- sum(all_abund)
  if (current_sum != n_individuals) {
    diff <- n_individuals - current_sum
    # Add to dominant species, as its count was the only one simply rounded.
    all_abund[1L] <- all_abund[1L] + diff
  }

  names(all_abund) <- LETTERS[seq_len(n_species)]
  all_abund[all_abund > 0]
}


#' Generate a heterogeneous community (point process + environment + interactions)
#'
#' @description
#' Generates individual locations and species identities over an arbitrary polygon
#' domain by combining:
#' \enumerate{
#'   \item a spatial point process for baseline locations (clustered, inhibited, or Poisson),
#'   \item environmental filtering for gradient-responsive species, and
#'   \item local interspecific interactions within a neighbourhood radius.
#' }
#' The dominant species "A" and the pool of non-dominant species can use
#' different point-process models. Environmental suitability is applied as a
#' Gaussian preference around a per-species optimum, and local interactions are
#' incorporated as a multiplicative modifier computed from nearby already-assigned
#' individuals.
#'
#' @details
#' \strong{Abundances.} Total individuals per species are generated with
#' \code{\link{generate_sad}()} using \code{P$SAD_MODEL} (default: \code{"fisher"}).
#' The built-in Fisher option allocates a fixed fraction to species "A" and
#' distributes the remainder by a log-series across \code{B, C, ...}; other SAD
#' models generate a full \code{A..} vector directly.
#'
#' \strong{Baseline locations (point processes).} Locations are simulated using
#' the process names in \code{P}. Supported values (case-insensitive):
#' \describe{
#'   \item{\code{"poisson"}}{Homogeneous Poisson process (Complete Spatial Randomness).}
#'   \item{\code{"thomas"}}{Thomas (Neyman-Scott) cluster process. A fast C++ implementation
#'   is used if available.}
#'   \item{\code{"strauss"}}{Strauss process for inhibition. A fast C++ MCMC implementation
#'   is used if available.}
#'   \item{\code{"geyer"}}{Geyer saturation process. A fast C++ MCMC implementation
#'   is used if available.}
#' }
#'
#' \strong{Environmental filtering.} For species listed in \code{P$GRADIENT},
#' the assignment probability for each species is proportional to
#' \eqn{\exp\{-(x - \mu)^2 / (2 \sigma^2)\}}, where \eqn{x} is the normalized
#' environmental value at the point, \eqn{\mu} the species optimum, and
#' \eqn{\sigma} the tolerance. Environmental values are attached by nearest
#' neighbour join to the grid returned by \code{\link{create_environmental_gradients}()}.
#'
#' \strong{Local interactions.} For each candidate point and focal species,
#' the abundance-independent interaction modifier is the geometric mean of the
#' corresponding coefficients in \code{P$INTERACTION_MATRIX} for neighbours found
#' within \code{P$INTERACTION_RADIUS} (using up to 5 nearest already-assigned
#' individuals). Coefficients \code{> 1} favour co-occurrence; \code{< 1} penalize it.
#'
#' \strong{Tie-breaking and robustness.} If all weights for a step are
#' non-finite or non-positive, a uniform assignment is used for that step to
#' avoid dead ends. Only points with a non-empty \code{species} are returned.
#'
#' @param domain An \pkg{sf} polygon (or multipolygon) defining the sampling
#'   domain; must have a valid CRS (projected coordinates recommended).
#' @param P A fully materialized parameter list, typically from
#'   \code{\link{load_config}()}, containing at least:
#'   \itemize{
#'     \item \emph{Community size:} \code{N_SPECIES}, \code{N_INDIVIDUALS},
#'       \code{DOMINANT_FRACTION}, \code{FISHER_ALPHA}, \code{FISHER_X}.
#'     \item \emph{Model family:} \code{MODEL_FAMILY} in \code{"manual"},
#'       \code{"niche_filtering"}, \code{"neutral_csr"},
#'       \code{"neutral_hubbell_like"}, or \code{"hybrid"}.
#'     \item \emph{Environment:} \code{SAMPLING_RESOLUTION}, \code{ENVIRONMENTAL_NOISE},
#'       and \code{GRADIENT} (tibble with \code{species}, \code{gradient} in
#'       \{temperature, elevation, rainfall\}, \code{optimum} in \eqn{[0,1]},
#'       and \code{tol > 0}). If absent, species are treated as neutral.
#'     \item \emph{Neutral/hybrid controls (if family is neutral/hybrid):}
#'       \code{NEUTRAL_M}, \code{NEUTRAL_NU}, \code{NEUTRAL_META_MODEL},
#'       \code{DISPERSAL_KERNEL}, \code{DISPERSAL_SCALE}, \code{DISPERSAL_ALPHA},
#'       and \code{HYBRID_ENV_WEIGHT}.
#'     \item \emph{Point-process selection (strings):}
#'       \code{SPATIAL_PROCESS_A} and \code{SPATIAL_PROCESS_OTHERS} in
#'       \code{"poisson"}, \code{"thomas"}, \code{"strauss"}, or \code{"geyer"}.
#'     \item \emph{Thomas (A) params (if used):}
#'       \code{A_PARENT_INTENSITY} (parents per area; optional),
#'       \code{A_MEAN_OFFSPRING} (mean children per parent),
#'       \code{A_CLUSTER_SCALE} (Gaussian sd of offspring displacement; map units).
#'     \item \emph{Strauss/Geyer (others) params (if used):}
#'       \itemize{
#'         \item \strong{Strauss (inhibition surrogate):} \code{OTHERS_R} (interaction radius) and
#'           \code{OTHERS_S} (inhibition strength in \eqn{(0,1]}; smaller values yield stronger inhibition).
#'         \item \strong{Geyer (saturation):} \code{OTHERS_R} (interaction radius), \code{OTHERS_GAMMA}
#'           (interaction parameter; \code{< 1} inhibition, \code{> 1} clustering), and \code{OTHERS_S}
#'           (saturation count; positive integer).
#'         \item \code{OTHERS_BETA} (baseline intensity/multiplier; optional; used by some engines).
#'       }
#'
#' \strong{OTHERS_* quick reference.} These parameters are shared across several
#' point-process code paths, so the same name can feed different underlying
#' samplers depending on \code{SPATIAL_PROCESS_OTHERS}:
#'
#' \tabular{lllll}{
#' \strong{Parameter} \tab \strong{Meaning} \tab \strong{Used when} \tab \strong{Default} \tab \strong{Constraints}\cr
#' \code{OTHERS_R} \tab interaction radius \tab Strauss, Geyer \tab \code{1} \tab \code{> 0} (map units)\cr
#' \code{OTHERS_S} \tab Strauss: inhibition strength; Geyer: saturation count \tab Strauss, Geyer \tab \code{2} \tab Strauss: \code{(0,1]}; Geyer: integer \code{>= 0}\cr
#' \code{OTHERS_GAMMA} \tab Geyer interaction parameter \tab Geyer \tab \code{NA} (falls back to \code{OTHERS_S} in some dispatchers) \tab \code{> 0} (\code{<1} inhibition, \code{>1} clustering)\cr
#' \code{OTHERS_BETA} \tab baseline intensity / multiplier \tab some engines / wrappers \tab \code{NA} \tab \code{>= 0} (units depend on engine)\cr
#' }
#'
#' \emph{Note:} in internal dispatchers, missing \code{OTHERS_GAMMA} may fall back
#' to \code{OTHERS_S} for backwards compatibility in some Strauss/Geyer paths.
#' Prefer setting \code{OTHERS_GAMMA} explicitly for Geyer.
#'     \item \emph{Local interactions:} \code{INTERACTION_RADIUS} (map units) and
#'       \code{INTERACTION_MATRIX} (S x S numeric, dimnames = species letters).
#'   }
#'
#' @return An \pkg{sf} POINT layer with a character column \code{species} and
#'   appended environmental columns (e.g. \code{temperature_C}, \code{elevation_m},
#'   \code{rainfall_mm}). Rows correspond to simulated individuals retained after
#'   assignment.
#'
#' @section Notes:
#' \itemize{
#'   \item Use a projected CRS (e.g. metres) so that process radii and cluster
#'   scales are in meaningful linear units.
#'   \item When \code{SPATIAL_PROCESS_OTHERS} is inhibitory and \code{n} is large
#'   relative to \code{OTHERS_R} and domain area, you may hit feasibility limits;
#'   adjust \code{n}, \code{r}, or choose a different process.
#'   \item Setting \code{INTERACTION_RADIUS = 0} or an all-ones matrix disables
#'   local interaction effects.
#' }
#'
#' @seealso
#' \code{\link{create_environmental_gradients}()},
#' \code{\link{generate_sad}()},
#' \code{\link{generate_fisher_log_series}()},
#' \code{\link{load_config}()},
#' \code{\link{load_interactions}()}
#'
#' @examples
#' \dontrun{
#' P <- load_config("simul_init.txt")
#' domain <- create_sampling_domain()
#'
#' # Example: A clustered (Thomas), others mildly inhibited (Strauss)
#' P$SPATIAL_PROCESS_A <- "thomas"
#' P$A_PARENT_INTENSITY <- NA
#' P$A_MEAN_OFFSPRING <- 10
#' P$A_CLUSTER_SCALE <- 0.8
#'
#' P$SPATIAL_PROCESS_OTHERS <- "strauss"
#' P$OTHERS_R <- 0.6
#' P$OTHERS_S <- 0.5
#'
#' pts <- generate_heterogeneous_distribution(domain, P)
#' plot(sf::st_geometry(pts), pch = 16, cex = 0.4)
#' }
#'
#' @export
generate_heterogeneous_distribution <- function(domain, P) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  family <- tolower(as.character(P$MODEL_FAMILY %||% "manual"))
  family <- gsub("-", "_", family)
  if (family %in% c("neutral_hubbell_like", "hybrid")) {
    return(spesim_simulate_neutral_recruitment(domain, P))
  }

  # --- 1) Target abundances --------------------------------------------
  abund <- generate_sad(
    n_species = P$N_SPECIES,
    n_individuals = P$N_INDIVIDUALS,
    model = P$SAD_MODEL %||% "fisher",
    # fisher
    dominant_fraction = P$DOMINANT_FRACTION,
    alpha = P$FISHER_ALPHA,
    x = P$FISHER_X,
    # geometric / zipf
    k = P$GEOMETRIC_K,
    exponent = P$ZIPF_EXPONENT,
    q = P$ZIPF_Q,
    # lognormal families
    meanlog = P$LOGNORMAL_MEANLOG,
    sdlog = P$LOGNORMAL_SDLOG,
    # poisson-gamma
    shape = P$POIGAMMA_SHAPE,
    rate = P$POIGAMMA_RATE,
    # neutral
    theta = P$ZSM_THETA,
    m = P$ZSM_M,
    # custom
    sad = P$SAD_VECTOR
  )
  n_A <- unname(abund["A"] %||% 0)
  other_species <- rep(
    names(abund)[names(abund) != "A"],
    times = abund[names(abund) != "A"]
  )

  # --- 2) Environmental grid (for later joins/weights) -----------------
  env_grid <- create_environmental_gradients(
    domain,
    P$SAMPLING_RESOLUTION,
    P$ENVIRONMENTAL_NOISE,
    covariates = P$ENV_COVARIATES %||% NULL,
    drivers = unique(c(as.character(P$ENV_DRIVERS %||% character(0)),
                       as.character(P$GRADIENT$gradient %||% character(0))))
  )
  env_sf <- sf::st_as_sf(
    env_grid,
    coords = c("x", "y"),
    crs = sf::st_crs(domain)
  )

  # --- 3) Simulate A / others -----------------------------------------
  crs_dom <- sf::st_crs(domain)

  # Dominant A
  proc_A <- tolower(P$SPATIAL_PROCESS_A %||% "poisson")
  args_A <- list(
    A_PARENT_INTENSITY = P$A_PARENT_INTENSITY,
    A_MEAN_OFFSPRING = P$A_MEAN_OFFSPRING,
    A_CLUSTER_SCALE = P$A_CLUSTER_SCALE
  )

  # Base points for species A
  pts_A_sf <- simulate_points_dispatch(proc_A, domain, n_A, args = args_A)

  # Optional environmental filtering for A: resample A locations from an
  # oversampled candidate set with weights proportional to Gaussian suitability.
  if (!is.null(P$GRADIENT) && NROW(P$GRADIENT) > 0 && ("A" %in% P$GRADIENT$species) && n_A > 0) {
    os <- as.numeric(P$ENV_FILTER_OVERSAMPLE %||% 6)
    if (!is.finite(os) || os < 1) os <- 6
    os <- as.integer(ceiling(os))

    rowA <- P$GRADIENT[P$GRADIENT$species == "A", , drop = FALSE][1, ]
    env_col <- .resolve_env_col(rowA$gradient, names(env_grid))

    if (!is.na(env_col)) {
      n_cand <- as.integer(max(n_A, n_A * os))

      # Generate a larger pool of points to sample from
      cand_sf <- simulate_points_dispatch(proc_A, domain, n_cand, args = args_A)

      if (nrow(cand_sf) < n_cand) {
        # top up if dispatcher returned fewer than requested
        cand_sf <- rbind(cand_sf, sf::st_sample(domain, size = n_cand - nrow(cand_sf), type = "random"))
      }

      cand_sf <- sf::st_join(cand_sf, env_sf, join = sf::st_nearest_feature)
      vals <- cand_sf[[env_col]]

      w <- exp(-((vals - rowA$optimum)^2) / (2 * rowA$tol^2))
      if (!all(is.finite(w)) || all(w <= 0)) w <- rep(1, length(w))

      sel <- sample.int(nrow(cand_sf), size = n_A, replace = FALSE, prob = w)

      # Overwrite pts_A_sf with the filtered sample
      pts_A_sf <- cand_sf[sel, ]
    }
  }

  # Others as one pool of points; species assigned later by counts
  n_all_others <- length(other_species)
  proc_O <- tolower(P$SPATIAL_PROCESS_OTHERS %||% "poisson")
  args_O <- list(
      OTHERS_MU = P$OTHERS_MU,
      OTHERS_SIGMA = P$OTHERS_SIGMA,
      OTHERS_BETA = P$OTHERS_BETA,
      OTHERS_R = P$OTHERS_R,
      OTHERS_S = P$OTHERS_S,
      OTHERS_GAMMA = P$OTHERS_GAMMA
  )
  pts_O_sf <- simulate_points_dispatch(proc_O, domain, n_all_others, args = args_O)

  # Extract sfc for the next step, which expects it.
  sfc_A <- sf::st_geometry(pts_A_sf)
  sfc_O <- sf::st_geometry(pts_O_sf)

  pts_A <- if (length(sfc_A) > 0) {
    sf::st_sf(species = rep("A", length(sfc_A)), geometry = sfc_A)
  } else {
    sf::st_sf(
      species = character(0),
      geometry = sf::st_sfc(crs = crs_dom)
    )
  }

  # For non-dominants, we initially store placeholder species labels and assign
  # species identities *after* point generation.
  pts_O <- if (length(sfc_O) > 0) {
    sf::st_sf(
      species = rep(NA_character_, length(sfc_O)),
      geometry = sfc_O
    )
  } else {
    sf::st_sf(
      species = character(0),
      geometry = sf::st_sfc(crs = crs_dom)
    )
  }

  pts <- rbind(pts_A, pts_O)
  if (nrow(pts) == 0) {
    return(sf::st_sf(
      species = character(0),
      geometry = sf::st_sfc(crs = crs_dom)
    ))
  }

  # --- 5) Species assignment: environment + optional local interactions -----
  #
  # The simulator generates locations first, then assigns species identities for
  # the non-dominant pool to match target abundances. For gradient-responsive
  # species, assignment is weighted by Gaussian environmental suitability.
  # If INTERACTION_RADIUS > 0 and INTERACTION_MATRIX is non-neutral, assignment
  # is also weighted by a local-interaction modifier computed from up to 5
  # nearest already-assigned neighbours within INTERACTION_RADIUS.

  # Attach environment once (used for weighted assignment; joined again at end)
  pts_env <- sf::st_join(pts, env_sf, join = sf::st_nearest_feature)

  # Identify which points are in the non-dominant pool
  idx_A <- which(pts_env$species == "A")
  idx_O <- which(is.na(pts_env$species))

  # Target non-dominant counts (exact)
  target_counts <- table(other_species)
  target_species <- names(target_counts)

  has_grad <- !is.null(P$GRADIENT) && NROW(P$GRADIENT) > 0
  use_interactions <- isTRUE(as.numeric(P$INTERACTION_RADIUS %||% 0) > 0) &&
    !is.null(P$INTERACTION_MATRIX) &&
    !isTRUE(all(as.numeric(P$INTERACTION_MATRIX) == 1))

  # If we have nothing to weight by, do a simple random permutation assignment.
  if (length(idx_O) > 0) {
    if (!has_grad && !use_interactions) {
      pts_env$species[idx_O] <- sample(other_species, size = length(idx_O), replace = FALSE)
    } else {
      # Precompute 5-nearest neighbours for all points (for interaction modifier)
      xy_all <- sf::st_coordinates(pts_env)
      k_nn <- min(5L, max(0L, nrow(xy_all) - 1L))
      nn_idx <- NULL
      nn_dist <- NULL
      if (k_nn >= 1L) {
        knn <- FNN::get.knn(xy_all, k = k_nn)
        nn_idx <- knn$nn.index
        nn_dist <- knn$nn.dist
      }

      # Remaining quota per species
      remain <- as.integer(target_counts)
      names(remain) <- target_species

      # Assigned labels (A fixed; others NA until assigned)
      assigned <- rep(NA_character_, nrow(pts_env))
      assigned[idx_A] <- "A"

      # Interaction matrix (ensure dimnames)
      IM <- P$INTERACTION_MATRIX
      if (is.null(IM)) {
        IM <- matrix(1, P$N_SPECIES, P$N_SPECIES,
          dimnames = list(LETTERS[1:P$N_SPECIES], LETTERS[1:P$N_SPECIES])
        )
      }
      if (is.null(dimnames(IM))) {
        dimnames(IM) <- list(LETTERS[1:nrow(IM)], LETTERS[1:ncol(IM)])
      }

      # Helper: environmental suitability for a focal species at point i
      env_weight <- function(sp, i) {
        if (!has_grad) return(1)
        row <- P$GRADIENT[P$GRADIENT$species == sp, , drop = FALSE]
        if (NROW(row) == 0) return(1)
        row <- row[1, ]
        env_col <- .resolve_env_col(row$gradient, names(env_grid))
        if (is.na(env_col) || !(env_col %in% names(pts_env))) return(1)
        val <- as.numeric(pts_env[[env_col]][i])
        if (!is.finite(val)) return(1)
        tol <- as.numeric(row$tol)
        if (!is.finite(tol) || tol <= 0) return(1)
        opt <- as.numeric(row$optimum)
        w <- exp(-((val - opt)^2) / (2 * tol^2))
        if (!is.finite(w) || w <= 0) return(1e-12)
        w
      }

      # Helper: local-interaction modifier using up to 5 nearest assigned neighbours
      interaction_weight <- function(sp, i) {
        if (!use_interactions || is.null(nn_idx)) return(1)
        # neighbours for i (precomputed KNN), filter to those already assigned and within radius
        neigh <- nn_idx[i, ]
        dist <- nn_dist[i, ]
        keep <- which(is.finite(dist) & dist <= as.numeric(P$INTERACTION_RADIUS))
        if (!length(keep)) return(1)
        neigh <- neigh[keep]
        neigh_sp <- assigned[neigh]
        neigh_sp <- neigh_sp[!is.na(neigh_sp)]
        if (!length(neigh_sp)) return(1)

        # geometric mean of coefficients
        vals <- as.numeric(IM[sp, neigh_sp, drop = TRUE])
        vals[!is.finite(vals) | vals <= 0] <- 1e-12
        exp(mean(log(vals)))
      }

      # Sequential assignment in random order ("already-assigned" neighbours)
      order_O <- sample(idx_O, size = length(idx_O), replace = FALSE)
      for (i in order_O) {
        avail <- names(remain)[remain > 0]
        if (!length(avail)) break
        if (length(avail) == 1L) {
          assigned[i] <- avail[1]
          remain[avail[1]] <- remain[avail[1]] - 1L
          next
        }

        w <- vapply(avail, function(sp) env_weight(sp, i) * interaction_weight(sp, i), numeric(1))
        w[!is.finite(w) | w < 0] <- 0
        if (sum(w) <= 0) {
          assigned[i] <- sample(avail, 1L)
        } else {
          assigned[i] <- sample(avail, 1L, prob = w)
        }
        remain[assigned[i]] <- remain[assigned[i]] - 1L
      }

      # Fill any remaining (should be rare) deterministically
      if (any(is.na(assigned[idx_O]))) {
        still <- idx_O[is.na(assigned[idx_O])]
        avail <- rep(names(remain), pmax(0L, remain))
        if (length(avail) == length(still)) {
          assigned[still] <- sample(avail, size = length(still), replace = FALSE)
        } else if (length(avail) > 0) {
          assigned[still] <- sample(avail, size = length(still), replace = TRUE)
        } else {
          assigned[still] <- sample(target_species, size = length(still), replace = TRUE)
        }
      }

      pts_env$species <- assigned
    }
  }

  # Return to a consistent schema
  pts$species <- as.character(pts_env$species)
# --- 6) Final env join and return -----------------------------------
  lin_st <- .linear_domain_state(domain, P)
  if (isTRUE(lin_st$enabled)) {
    xy_all <- sf::st_coordinates(pts)
    xy_fix <- .constrain_xy_to_linear(xy_all, lin_st, domain)
    sf::st_geometry(pts) <- .as_sfc_points(xy_fix, crs_dom)
    pts$linear_pos <- .xy_to_linear_pos(xy_fix, lin_st)
  }
  pts <- sf::st_join(pts, env_sf, join = sf::st_nearest_feature)
  pts$species <- as.character(pts$species)
  pts
}
