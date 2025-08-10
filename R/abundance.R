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
generate_fisher_log_series <- function(n_species, n_individuals, dominant_fraction, alpha, x) {
  stopifnot(n_species >= 1L, n_individuals >= 1L)
  n_dominant <- round(n_individuals * dominant_fraction)
  n_remaining <- n_individuals - n_dominant

  if (n_species == 1L) {
    out <- c(n_dominant)
    names(out) <- "A"
    return(out)
  }

  ranks <- 2:n_species
  rel <- alpha * (x^ranks) / ranks
  abund <- if (sum(rel) > 0) round(rel / sum(rel) * n_remaining) else rep(0L, length(rel))

  all_abund <- c(n_dominant, abund)
  adj <- n_individuals - sum(all_abund)
  if (length(all_abund) >= 2L) all_abund[2L] <- all_abund[2L] + adj else all_abund[1L] <- all_abund[1L] + adj

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
#' \strong{Abundances.} Total individuals per species are drawn with
#' \code{\link{generate_fisher_log_series}()}, allocating a fixed fraction to
#' species "A" and distributing the remainder by a log-series across
#' \code{B, C, ...}.
#'
#' \strong{Baseline locations (point processes).} Locations are simulated using
#' the process names in \code{P}. Supported values (case-insensitive):
#' \describe{
#'   \item{\code{"poisson"}}{Homogeneous Poisson inside \code{domain}.}
#'   \item{\code{"thomas"}}{Thomas (Neyman-Scott) cluster process; uses parent intensity
#'   (derived if missing), mean offspring, and Gaussian cluster scale. Implemented internally.}
#'   \item{\code{"strauss"}}{Inhibition via a sequential-inhibition surrogate using
#'   interaction radius \code{OTHERS_R} and inhibition strength \code{OTHERS_S} in \eqn{(0,1]}.
#'   Lower \code{OTHERS_S} increases inhibition.}
#'   \item{\code{"geyer"}}{Geyer saturation surrogate: acceptance probability proportional to
#'   \eqn{\gamma^{\min(m, s)}} where \eqn{m} neighbours fall within \code{OTHERS_R} and
#'   saturation count \code{OTHERS_S}. Implemented internally.}
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
#'     \item \emph{Environment:} \code{SAMPLING_RESOLUTION}, \code{ENVIRONMENTAL_NOISE},
#'       and \code{GRADIENT} (tibble with \code{species}, \code{gradient} in
#'       \{temperature, elevation, rainfall\}, \code{optimum} in \eqn{[0,1]},
#'       and \code{tol > 0}). If absent, species are treated as neutral.
#'     \item \emph{Point-process selection (strings):}
#'       \code{SPATIAL_PROCESS_A} and \code{SPATIAL_PROCESS_OTHERS} in
#'       \code{"poisson"}, \code{"thomas"}, \code{"strauss"}, or \code{"geyer"}.
#'     \item \emph{Thomas (A) params (if used):}
#'       \code{A_PARENT_INTENSITY} (parents per area; optional),
#'       \code{A_MEAN_OFFSPRING} (mean children per parent),
#'       \code{A_CLUSTER_SCALE} (Gaussian sd of offspring displacement; map units).
#'     \item \emph{Strauss/Geyer (others) params (if used):}
#'       \code{OTHERS_R} (interaction radius),
#'       \code{OTHERS_GAMMA} (interaction parameter; \code{< 1} inhibition, \code{> 1} attraction),
#'       \code{OTHERS_S} (Geyer saturation count; ignored by Strauss),
#'       \code{OTHERS_BETA} (baseline intensity/multiplier; optional).
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

  # --- small helpers ---------------------------------------------------
  .in_domain <- function(xy, poly) {
    if (NROW(xy) == 0) {
      return(logical(0))
    }
    pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]), coords = c("x", "y"), crs = sf::st_crs(poly))
    as.logical(sf::st_within(pts, poly, sparse = FALSE)[, 1])
  }
  .sample_uniform_in_domain <- function(n, poly, max_tries = n * 50) {
    bb <- sf::st_bbox(poly)
    res <- matrix(NA_real_, 0, 2)
    tries <- 0L
    while (nrow(res) < n && tries < max_tries) {
      need <- n - nrow(res)
      cand <- cbind(
        stats::runif(need, bb["xmin"], bb["xmax"]),
        stats::runif(need, bb["ymin"], bb["ymax"])
      )
      keep <- .in_domain(cand, poly)
      if (any(keep)) res <- rbind(res, cand[keep, , drop = FALSE])
      tries <- tries + 1L
    }
    if (nrow(res) > n) res <- res[seq_len(n), , drop = FALSE]
    res
  }
  .simulate_thomas_points <- function(n, poly, mu = 10, sigma = 1, kappa = NA_real_) {
    if (n <= 0) {
      return(matrix(numeric(0), ncol = 2))
    }
    area <- as.numeric(sf::st_area(sf::st_as_sf(poly)[1, ]))
    if (!is.finite(kappa) || kappa <= 0) kappa <- max(1e-6, n / (mu * area))
    # parents
    n_par <- stats::rpois(1L, lambda = kappa * area)
    par_xy <- .sample_uniform_in_domain(max(1L, n_par), poly)
    # offspring
    out <- matrix(numeric(0), ncol = 2)
    if (nrow(par_xy) > 0) {
      # Poisson(mu) per parent
      k_vec <- stats::rpois(nrow(par_xy), mu)
      if (sum(k_vec) > 0) {
        ox <- rep(par_xy[, 1], times = k_vec) + stats::rnorm(sum(k_vec), 0, sigma)
        oy <- rep(par_xy[, 2], times = k_vec) + stats::rnorm(sum(k_vec), 0, sigma)
        off <- cbind(ox, oy)
        keep <- .in_domain(off, poly)
        if (any(keep)) out <- off[keep, , drop = FALSE]
      }
    }
    if (nrow(out) == 0) out <- .sample_uniform_in_domain(n, poly)
    if (nrow(out) > n) out <- out[sample.int(nrow(out), n), , drop = FALSE]
    out
  }
  .simulate_strauss_points <- function(n, poly, r, s = 0.7) {
    # sequential inhibition with hard-core radius r_eff ~ r * (1 - 0.5*(1 - s))
    if (n <= 0) {
      return(matrix(numeric(0), ncol = 2))
    }
    bb <- sf::st_bbox(poly)
    r_eff <- max(1e-6, as.numeric(r) * (0.5 + 0.5 * s)) # s in (0,1]; smaller s => larger r_eff
    out <- matrix(NA_real_, 0, 2)
    tries <- 0L
    max_tries <- 2000L + 50L * n
    while (nrow(out) < n && tries < max_tries) {
      cand <- cbind(
        stats::runif(1, bb["xmin"], bb["xmax"]),
        stats::runif(1, bb["ymin"], bb["ymax"])
      )
      if (.in_domain(cand, poly)) {
        ok <- TRUE
        if (nrow(out) > 0) {
          d2 <- colSums((t(out) - as.numeric(cand))^2)
          ok <- min(d2) >= r_eff^2
        }
        if (ok) out <- rbind(out, cand)
      }
      tries <- tries + 1L
    }
    if (nrow(out) < n) {
      # top up with uniform if infeasible
      need <- n - nrow(out)
      out <- rbind(out, .sample_uniform_in_domain(need, poly))
    }
    out
  }
  .simulate_geyer_points <- function(n, poly, r, gamma = 1.5, s = 2) {
    if (n <= 0) {
      return(matrix(numeric(0), ncol = 2))
    }
    # start from Poisson candidates; accept with prob ~ gamma^{min(m, s)}
    cand <- .sample_uniform_in_domain(3L * n, poly)
    if (nrow(cand) == 0) {
      return(cand)
    }
    out <- matrix(NA_real_, 0, 2)
    for (i in seq_len(nrow(cand))) {
      xy <- cand[i, , drop = FALSE]
      m <- 0L
      if (nrow(out) > 0) {
        d2 <- colSums((t(out) - as.numeric(xy))^2)
        m <- sum(d2 <= (as.numeric(r)^2))
      }
      acc_prob <- (as.numeric(gamma))^(min(m, as.integer(s)))
      acc_prob <- max(0, min(1, acc_prob / (1 + acc_prob))) # squash to (0,1)
      if (stats::runif(1) < acc_prob) out <- rbind(out, xy)
      if (nrow(out) >= n) break
    }
    if (nrow(out) < n) {
      top <- .sample_uniform_in_domain(n - nrow(out), poly)
      out <- rbind(out, top)
    }
    out
  }
  .as_sfc_points <- function(xy, crs) {
    if (NROW(xy) == 0) {
      return(sf::st_sfc(crs = crs))
    }
    sf::st_sfc(lapply(seq_len(NROW(xy)), function(i) sf::st_point(as.numeric(xy[i, 1:2]))), crs = crs)
  }

  # --- 1) Target abundances --------------------------------------------
  abund <- generate_fisher_log_series(
    P$N_SPECIES, P$N_INDIVIDUALS,
    P$DOMINANT_FRACTION, P$FISHER_ALPHA, P$FISHER_X
  )
  n_A <- unname(abund["A"] %||% 0)
  other_species <- rep(names(abund)[names(abund) != "A"], times = abund[names(abund) != "A"])

  # --- 2) Environmental grid (for later joins/weights) -----------------
  env_grid <- create_environmental_gradients(domain, P$SAMPLING_RESOLUTION, P$ENVIRONMENTAL_NOISE)
  env_sf <- sf::st_as_sf(env_grid, coords = c("x", "y"), crs = sf::st_crs(domain))

  # --- 3) Simulate A / others -----------------------------------------
  crs_dom <- sf::st_crs(domain)

  # Dominant A
  proc_A <- tolower(P$SPATIAL_PROCESS_A %||% "poisson")
  xy_A <- switch(proc_A,
    "poisson" = .sample_uniform_in_domain(n_A, domain),
    "thomas" = .simulate_thomas_points(n_A, domain,
      mu    = as.numeric(P$A_MEAN_OFFSPRING %||% 10),
      sigma = as.numeric(P$A_CLUSTER_SCALE %||% 1),
      kappa = as.numeric(P$A_PARENT_INTENSITY %||% NA_real_)
    ),
    # fallback
    .sample_uniform_in_domain(n_A, domain)
  )

  # Others as one pool of points; species assigned later by counts
  n_all_others <- length(other_species)
  proc_O <- tolower(P$SPATIAL_PROCESS_OTHERS %||% "poisson")
  xy_O <- switch(proc_O,
    "poisson" = .sample_uniform_in_domain(n_all_others, domain),
    "thomas" = .simulate_thomas_points(n_all_others, domain,
      mu    = as.numeric(P$OTHERS_MU %||% 10),
      sigma = as.numeric(P$OTHERS_SIGMA %||% 1),
      kappa = as.numeric(P$OTHERS_BETA %||% NA_real_)
    ), # treat as kappa if given
    "strauss" = .simulate_strauss_points(n_all_others, domain,
      r = as.numeric(P$OTHERS_R %||% 1),
      s = as.numeric(P$OTHERS_S %||% 0.7)
    ),
    "geyer" = .simulate_geyer_points(n_all_others, domain,
      r = as.numeric(P$OTHERS_R %||% 1),
      gamma = as.numeric(P$OTHERS_GAMMA %||% 1.5),
      s = as.numeric(P$OTHERS_S %||% 2)
    ),
    .sample_uniform_in_domain(n_all_others, domain)
  )

  # --- 4) Build sf with consistent schema ------------------------------
  sfc_A <- .as_sfc_points(xy_A, crs_dom)
  sfc_O <- .as_sfc_points(xy_O, crs_dom)

  pts_A <- if (length(sfc_A) > 0) {
    sf::st_sf(species = factor(rep("A", length(sfc_A))), geometry = sfc_A)
  } else {
    sf::st_sf(species = factor(character(0)), geometry = sf::st_sfc(crs = crs_dom))
  }
  pts_O <- if (length(sfc_O) > 0) {
    sf::st_sf(species = factor(other_species, levels = LETTERS[1:P$N_SPECIES]), geometry = sfc_O)
  } else {
    sf::st_sf(
      species = factor(character(0), levels = LETTERS[1:P$N_SPECIES]),
      geometry = sf::st_sfc(crs = crs_dom)
    )
  }

  pts <- rbind(pts_A, pts_O)
  if (nrow(pts) == 0) {
    return(sf::st_sf(
      species = factor(character(0), levels = LETTERS[1:P$N_SPECIES]),
      geometry = sf::st_sfc(crs = crs_dom)
    ))
  }

  # --- 5) Importance reweighting by environment (keeps counts) ---------
  if (!is.null(P$GRADIENT) && NROW(P$GRADIENT) > 0) {
    tmp <- sf::st_join(pts, env_sf, join = sf::st_nearest_feature)
    resel_idx <- integer(0)
    sp_levels <- unique(as.character(pts$species))
    for (sp in sp_levels) {
      ids <- which(tmp$species == sp)
      q <- length(ids)
      if (q == 0) next

      probs <- rep(1, q)
      if (sp %in% P$GRADIENT$species) {
        row <- P$GRADIENT[P$GRADIENT$species == sp, , drop = FALSE][1, ]
        env_col <- switch(row$gradient,
          temperature = "temperature",
          elevation   = "elevation",
          rainfall    = "rainfall"
        )
        vals <- tmp[[env_col]][ids]
        probs <- exp(-((vals - row$optimum)^2) / (2 * row$tol^2))
        if (!all(is.finite(probs)) || all(probs <= 0)) probs <- rep(1, q)
      }
      resel_idx <- c(resel_idx, sample(ids, size = q, prob = probs, replace = FALSE))
    }
    pts <- pts[sort(resel_idx), , drop = FALSE]
  }

  # --- 6) Final env join and return -----------------------------------
  pts <- sf::st_join(pts, env_sf, join = sf::st_nearest_feature)
  pts$species <- as.character(pts$species)
  pts
}
