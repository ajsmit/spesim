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
  n_dominant <- round(n_individuals * dominant_fraction)
  n_remaining <- n_individuals - n_dominant

  if (n_species == 1L) {
    # With a single species, all individuals are assigned to A.
    out <- c(as.integer(n_individuals))
    names(out) <- "A"
    return(out)
  }

  ranks <- 2:n_species
  rel <- alpha * (x^ranks) / ranks
  abund <- if (sum(rel) > 0) {
    round(rel / sum(rel) * n_remaining)
  } else {
    rep(0L, length(rel))
  }

  all_abund <- c(n_dominant, abund)
  adj <- n_individuals - sum(all_abund)
  if (length(all_abund) >= 2L) {
    all_abund[2L] <- all_abund[2L] + adj
  } else {
    all_abund[1L] <- all_abund[1L] + adj
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
#'   \item{\code{"poisson"}}{Homogeneous Poisson inside \code{domain}.}
#'   \item{\code{"thomas"}}{Thomas (Neyman-Scott) cluster process; uses parent intensity
#'   (derived if missing), mean offspring, and Gaussian cluster scale. Implemented internally.}
#'   \item{\code{"strauss"}}{Inhibition via a sequential-inhibition surrogate using
#'   interaction radius \code{OTHERS_R} and inhibition strength \code{OTHERS_S} in \eqn{(0,1]}.
#'   Lower \code{OTHERS_S} increases inhibition (a larger effective hard-core distance).}
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

  # --- small helpers ---------------------------------------------------
  .in_domain <- function(xy, poly) {
    if (NROW(xy) == 0) {
      return(logical(0))
    }
    pts <- sf::st_as_sf(
      data.frame(x = xy[, 1], y = xy[, 2]),
      coords = c("x", "y"),
      crs = sf::st_crs(poly)
    )
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
      if (any(keep)) {
        res <- rbind(res, cand[keep, , drop = FALSE])
      }
      tries <- tries + 1L
    }
    if (nrow(res) > n) {
      res <- res[seq_len(n), , drop = FALSE]
    }
    res
  }
  .simulate_thomas_points <- function(
    n,
    poly,
    mu = 10,
    sigma = 1,
    kappa = NA_real_
  ) {
    if (n <= 0) {
      return(matrix(numeric(0), ncol = 2))
    }
    area <- as.numeric(sf::st_area(sf::st_as_sf(poly)[1, ]))
    if (!is.finite(kappa) || kappa <= 0) {
      kappa <- max(1e-6, n / (mu * area))
    }
    # parents
    n_par <- stats::rpois(1L, lambda = kappa * area)
    par_xy <- .sample_uniform_in_domain(max(1L, n_par), poly)
    # offspring
    out <- matrix(numeric(0), ncol = 2)
    if (nrow(par_xy) > 0) {
      # Poisson(mu) per parent
      k_vec <- stats::rpois(nrow(par_xy), mu)
      if (sum(k_vec) > 0) {
        ox <- rep(par_xy[, 1], times = k_vec) +
          stats::rnorm(sum(k_vec), 0, sigma)
        oy <- rep(par_xy[, 2], times = k_vec) +
          stats::rnorm(sum(k_vec), 0, sigma)
        off <- cbind(ox, oy)
        keep <- .in_domain(off, poly)
        if (any(keep)) out <- off[keep, , drop = FALSE]
      }
    }
    if (nrow(out) == 0) {
      out <- .sample_uniform_in_domain(n, poly)
    }
    if (nrow(out) > n) {
      out <- out[sample.int(nrow(out), n), , drop = FALSE]
    }
    out
  }
  .simulate_strauss_points <- function(n, poly, r, s = 0.7) {
    # Sequential inhibition surrogate with an effective hard-core radius r_eff.
    # We treat s in (0,1] as an *inhibition strength* parameter: smaller s => stronger inhibition.
    if (n <= 0) {
      return(matrix(numeric(0), ncol = 2))
    }
    bb <- sf::st_bbox(poly)
    s <- as.numeric(s)
    if (!is.finite(s)) s <- 0.7
    s <- max(1e-6, min(1, s))
    # Map s -> r_eff such that s=1 gives r_eff=r, and s->0 increases r_eff up to ~1.5*r.
    r_eff <- max(1e-6, as.numeric(r) * (1.5 - 0.5 * s))
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
      if (stats::runif(1) < acc_prob) {
        out <- rbind(out, xy)
      }
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
    sf::st_sfc(
      lapply(seq_len(NROW(xy)), function(i) {
        sf::st_point(as.numeric(xy[i, 1:2]))
      }),
      crs = crs
    )
  }
  .resolve_env_col <- function(gname, cols) {
    g <- as.character(gname %||% "")
    if (!nzchar(g)) return(NA_character_)
    if (g %in% cols) return(g)
    legacy <- switch(g,
      temperature = "temperature_C",
      elevation = "elevation_m",
      rainfall = "rainfall_mm",
      NA_character_
    )
    if (!is.na(legacy) && legacy %in% cols) return(legacy)
    cand <- grep(paste0("^", g, "(_|$)"), cols, value = TRUE)
    if (length(cand)) return(cand[1])
    NA_character_
  }
  .linear_domain_state <- function(poly, P) {
    d_type <- tolower(as.character(P$DOMAIN_TYPE %||% "polygon"))
    d_type <- gsub("-", "_", d_type)
    enabled <- d_type %in% c("network", "coastline")
    bb <- sf::st_bbox(poly)
    axis <- tolower(as.character(P$LINEAR_AXIS %||% "x"))
    if (!axis %in% c("x", "y")) axis <- "x"
    wrap <- isTRUE(P$LINEAR_WRAP)
    jitter_sd <- as.numeric(P$LINEAR_JITTER_SD %||% 0)
    if (!is.finite(jitter_sd) || jitter_sd < 0) jitter_sd <- 0
    list(
      enabled = enabled,
      type = d_type,
      axis = axis,
      wrap = wrap,
      jitter_sd = jitter_sd,
      xmin = as.numeric(bb["xmin"]),
      xmax = as.numeric(bb["xmax"]),
      ymin = as.numeric(bb["ymin"]),
      ymax = as.numeric(bb["ymax"])
    )
  }
  .xy_to_linear_pos <- function(xy, st) {
    if (NROW(xy) == 0) return(numeric(0))
    if (!isTRUE(st$enabled)) return(rep(NA_real_, NROW(xy)))
    if (identical(st$axis, "x")) {
      den <- st$xmax - st$xmin
      if (!is.finite(den) || den == 0) den <- 1
      p <- (xy[, 1] - st$xmin) / den
    } else {
      den <- st$ymax - st$ymin
      if (!is.finite(den) || den == 0) den <- 1
      p <- (xy[, 2] - st$ymin) / den
    }
    if (isTRUE(st$wrap)) {
      p %% 1
    } else {
      pmax(0, pmin(1, p))
    }
  }
  .linear_pos_to_xy <- function(pos, st, poly) {
    if (length(pos) == 0L) return(matrix(numeric(0), ncol = 2))
    p <- as.numeric(pos)
    p[!is.finite(p)] <- 0.5
    if (isTRUE(st$wrap)) {
      p <- p %% 1
    } else {
      p <- pmax(0, pmin(1, p))
    }
    n <- length(p)
    out <- matrix(NA_real_, nrow = n, ncol = 2)
    ymid <- (st$ymin + st$ymax) / 2
    xmid <- (st$xmin + st$xmax) / 2
    if (identical(st$axis, "x")) {
      out[, 1] <- st$xmin + p * (st$xmax - st$xmin)
      out[, 2] <- ymid + stats::rnorm(n, 0, st$jitter_sd)
    } else {
      out[, 2] <- st$ymin + p * (st$ymax - st$ymin)
      out[, 1] <- xmid + stats::rnorm(n, 0, st$jitter_sd)
    }
    # ensure points remain inside polygonal domain
    bad <- which(!.in_domain(out, poly))
    if (length(bad) > 0) {
      out[bad, ] <- .sample_uniform_in_domain(length(bad), poly)
    }
    out
  }
  .constrain_xy_to_linear <- function(xy, st, poly) {
    if (!isTRUE(st$enabled) || NROW(xy) == 0) return(xy)
    p <- .xy_to_linear_pos(xy, st)
    .linear_pos_to_xy(p, st, poly)
  }
  .simulate_neutral_recruitment <- function(poly, P) {
    J <- as.integer(P$N_INDIVIDUALS %||% 0L)
    S <- as.integer(P$N_SPECIES %||% 1L)
    spp <- LETTERS[seq_len(max(1L, S))]
    crs_dom <- sf::st_crs(poly)

    if (!is.finite(J) || J < 0L) J <- 0L
    if (J == 0L) {
      return(sf::st_sf(species = character(0), geometry = sf::st_sfc(crs = crs_dom)))
    }

    family <- tolower(as.character(P$MODEL_FAMILY %||% "manual"))
    family <- gsub("-", "_", family)
    hybrid_mode <- identical(family, "hybrid")
    linear_st <- .linear_domain_state(poly, P)

    m <- as.numeric(P$NEUTRAL_M %||% 0.1)
    nu <- as.numeric(P$NEUTRAL_NU %||% 0.0)
    m <- min(1, max(0, ifelse(is.finite(m), m, 0.1)))
    nu <- min(1, max(0, ifelse(is.finite(nu), nu, 0.0)))

    kernel <- tolower(as.character(P$DISPERSAL_KERNEL %||% "gaussian"))
    kernel <- gsub("-", "_", kernel)
    if (!kernel %in% c("gaussian", "exponential", "power_law")) kernel <- "gaussian"
    scale <- as.numeric(P$DISPERSAL_SCALE %||% 0.5)
    if (!is.finite(scale) || scale <= 0) scale <- 0.5
    alpha <- as.numeric(P$DISPERSAL_ALPHA %||% 2)
    if (!is.finite(alpha) || alpha <= 0) alpha <- 2
    dir_bias <- as.numeric(P$DISPERSAL_DIRECTION_BIAS %||% 0)
    if (!is.finite(dir_bias)) dir_bias <- 0
    dir_bias <- max(-1, min(1, dir_bias))
    lambda_env <- as.numeric(P$HYBRID_ENV_WEIGHT %||% 1)
    if (!is.finite(lambda_env) || lambda_env < 0) lambda_env <- 1

    meta_model <- as.character(P$NEUTRAL_META_MODEL %||% P$SAD_MODEL %||% "zsm")
    meta_counts <- generate_sad(
      n_species = length(spp),
      n_individuals = max(J, length(spp)),
      model = meta_model,
      dominant_fraction = P$DOMINANT_FRACTION,
      alpha = P$FISHER_ALPHA,
      x = P$FISHER_X,
      k = P$GEOMETRIC_K,
      exponent = P$ZIPF_EXPONENT,
      q = P$ZIPF_Q,
      meanlog = P$LOGNORMAL_MEANLOG,
      sdlog = P$LOGNORMAL_SDLOG,
      shape = P$POIGAMMA_SHAPE,
      rate = P$POIGAMMA_RATE,
      theta = P$ZSM_THETA,
      m = P$ZSM_M,
      sad = P$SAD_VECTOR
    )
    meta_prob <- as.numeric(meta_counts)
    if (!all(is.finite(meta_prob)) || sum(meta_prob) <= 0) {
      meta_prob <- rep(1 / length(spp), length(spp))
    } else {
      meta_prob <- meta_prob / sum(meta_prob)
    }

    drv <- unique(c(as.character(P$ENV_DRIVERS %||% character(0)),
                    as.character(P$GRADIENT$gradient %||% character(0))))
    env_grid <- create_environmental_gradients(
      poly,
      P$SAMPLING_RESOLUTION,
      P$ENVIRONMENTAL_NOISE,
      covariates = P$ENV_COVARIATES %||% NULL,
      drivers = drv
    )
    env_sf <- sf::st_as_sf(env_grid, coords = c("x", "y"), crs = crs_dom)
    gx <- sort(unique(env_grid$x))
    gy <- sort(unique(env_grid$y))
    ix <- match(env_grid$x, gx)
    iy <- match(env_grid$y, gy)
    .to_grid_index <- function(v, grid_vals) {
      if (length(grid_vals) <= 1L) return(1L)
      mids <- (grid_vals[-1L] + grid_vals[-length(grid_vals)]) / 2
      findInterval(v, mids) + 1L
    }
    .build_grid <- function(colname) {
      mat <- matrix(NA_real_, nrow = length(gy), ncol = length(gx))
      mat[cbind(iy, ix)] <- as.numeric(env_grid[[colname]])
      mat
    }
    temp_mat <- .build_grid("temperature")
    elev_mat <- .build_grid("elevation")
    rain_mat <- .build_grid("rainfall")

    grad_lookup <- list()
    if (hybrid_mode && !is.null(P$GRADIENT) && NROW(P$GRADIENT) > 0) {
      for (i in seq_len(NROW(P$GRADIENT))) {
        row <- P$GRADIENT[i, , drop = FALSE]
        grad_lookup[[as.character(row$species)]] <- list(
          gradient = as.character(row$gradient),
          optimum = as.numeric(row$optimum),
          tol = as.numeric(row$tol)
        )
      }
    }

    .sample_immigrant_species <- function() {
      sample(spp, size = 1L, prob = meta_prob)
    }
    .step_magnitude <- function() {
      if (kernel == "gaussian") {
        abs(stats::rnorm(1L, 0, scale))
      } else if (kernel == "exponential") {
        stats::rexp(1L, rate = 1 / scale)
      } else {
        u <- stats::runif(1L)
        scale * ((1 - u)^(-1 / alpha) - 1)
      }
    }
    .propose_displacement <- function(parent_xy) {
      ang <- stats::runif(1L, 0, 2 * pi)
      r <- .step_magnitude()
      dx <- r * cos(ang)
      dy <- r * sin(ang)
      c(parent_xy[1] + dx, parent_xy[2] + dy)
    }
    .propose_linear_pos <- function(parent_pos) {
      mag <- .step_magnitude()
      axis_range <- if (identical(linear_st$axis, "x")) linear_st$xmax - linear_st$xmin else linear_st$ymax - linear_st$ymin
      if (!is.finite(axis_range) || axis_range <= 0) axis_range <- 1
      step <- mag / axis_range
      p_plus <- (1 + dir_bias) / 2
      sgn <- if (stats::runif(1) < p_plus) 1 else -1
      out <- parent_pos + sgn * step
      if (isTRUE(linear_st$wrap)) {
        out %% 1
      } else {
        max(0, min(1, out))
      }
    }
    .is_in_domain <- function(xy) {
      .in_domain(matrix(as.numeric(xy), ncol = 2), poly)[1]
    }
    .env_accept_prob <- function(sp, xy) {
      if (!hybrid_mode || lambda_env == 0 || is.null(grad_lookup[[sp]])) return(1)
      g <- grad_lookup[[sp]]
      env_col <- .resolve_env_col(g$gradient, names(env_grid))
      if (is.na(env_col)) return(1)
      ixx <- .to_grid_index(as.numeric(xy[1]), gx)
      iyy <- .to_grid_index(as.numeric(xy[2]), gy)
      val <- if (env_col == "temperature") {
        temp_mat[iyy, ixx]
      } else if (env_col == "elevation") {
        elev_mat[iyy, ixx]
      } else if (env_col == "rainfall") {
        rain_mat[iyy, ixx]
      } else {
        pt <- sf::st_as_sf(
          data.frame(x = as.numeric(xy[1]), y = as.numeric(xy[2])),
          coords = c("x", "y"),
          crs = crs_dom
        )
        j <- sf::st_nearest_feature(pt, env_sf)
        as.numeric(env_sf[[env_col]][j])
      }
      if (!is.finite(val) || !is.finite(g$tol) || g$tol <= 0) return(1)
      w <- exp(-((val - g$optimum)^2) / (2 * g$tol^2))
      w <- max(0, min(1, w^lambda_env))
      w
    }

    xy <- matrix(NA_real_, nrow = J, ncol = 2)
    linear_pos <- rep(NA_real_, J)
    spp_out <- character(J)

    for (i in seq_len(J)) {
      is_immigrant <- (i == 1L) || (stats::runif(1L) < m)
      parent_idx <- NA_integer_
      sp_i <- NA_character_
      if (is_immigrant) {
        sp_i <- .sample_immigrant_species()
      } else {
        parent_idx <- sample.int(i - 1L, 1L)
        sp_i <- spp_out[parent_idx]
        if (stats::runif(1L) < nu) {
          sp_i <- .sample_immigrant_species()
        }
      }

      accepted <- FALSE
      for (att in seq_len(40L)) {
        cand_xy <- if (is_immigrant || !is.finite(parent_idx)) {
          if (isTRUE(linear_st$enabled)) {
            p_i <- stats::runif(1)
            linear_pos[i] <- p_i
            .linear_pos_to_xy(p_i, linear_st, poly)
          } else {
            .sample_uniform_in_domain(1L, poly)
          }
        } else {
          if (isTRUE(linear_st$enabled)) {
            p_i <- .propose_linear_pos(linear_pos[parent_idx])
            linear_pos[i] <- p_i
            .linear_pos_to_xy(p_i, linear_st, poly)
          } else {
            cxy <- .propose_displacement(xy[parent_idx, ])
            if (.is_in_domain(cxy)) matrix(cxy, ncol = 2) else matrix(numeric(0), ncol = 2)
          }
        }
        if (NROW(cand_xy) == 0L) next
        p_acc <- .env_accept_prob(sp_i, cand_xy[1, ])
        if (stats::runif(1L) <= p_acc) {
          xy[i, ] <- cand_xy[1, ]
          spp_out[i] <- sp_i
          accepted <- TRUE
          break
        }
      }

      if (!accepted) {
        fallback <- if (isTRUE(linear_st$enabled)) {
          p_i <- stats::runif(1)
          linear_pos[i] <- p_i
          .linear_pos_to_xy(p_i, linear_st, poly)
        } else {
          .sample_uniform_in_domain(1L, poly)
        }
        xy[i, ] <- fallback[1, ]
        spp_out[i] <- sp_i
      }
    }

    out <- sf::st_sf(
      species = as.character(spp_out),
      geometry = .as_sfc_points(xy, crs_dom)
    )
    if (isTRUE(linear_st$enabled)) {
      out$linear_pos <- linear_pos
    }
    sf::st_join(out, env_sf, join = sf::st_nearest_feature)
  }

  family <- tolower(as.character(P$MODEL_FAMILY %||% "manual"))
  family <- gsub("-", "_", family)
  if (family %in% c("neutral_hubbell_like", "hybrid")) {
    return(.simulate_neutral_recruitment(domain, P))
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
  xy_A <- switch(
    proc_A,
    "poisson" = .sample_uniform_in_domain(n_A, domain),
    "thomas" = .simulate_thomas_points(
      n_A,
      domain,
      mu = as.numeric(P$A_MEAN_OFFSPRING %||% 10),
      sigma = as.numeric(P$A_CLUSTER_SCALE %||% 1),
      kappa = as.numeric(P$A_PARENT_INTENSITY %||% NA_real_)
    ),
    # fallback
    .sample_uniform_in_domain(n_A, domain)
  )

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
      xy_cand <- switch(proc_A,
        "poisson" = .sample_uniform_in_domain(n_cand, domain),
        "thomas" = .simulate_thomas_points(n_cand, domain,
          mu    = as.numeric(P$A_MEAN_OFFSPRING %||% 10),
          sigma = as.numeric(P$A_CLUSTER_SCALE %||% 1),
          kappa = as.numeric(P$A_PARENT_INTENSITY %||% NA_real_)
        ),
        .sample_uniform_in_domain(n_cand, domain)
      )
      if (NROW(xy_cand) < n_cand) {
        xy_cand <- rbind(xy_cand, .sample_uniform_in_domain(n_cand - NROW(xy_cand), domain))
      }

      cand_sf <- sf::st_sf(species = "A", geometry = .as_sfc_points(xy_cand, crs_dom))
      cand_sf <- sf::st_join(cand_sf, env_sf, join = sf::st_nearest_feature)
      vals <- cand_sf[[env_col]]

      w <- exp(-((vals - rowA$optimum)^2) / (2 * rowA$tol^2))
      if (!all(is.finite(w)) || all(w <= 0)) w <- rep(1, length(w))

      sel <- sample.int(length(w), size = n_A, replace = FALSE, prob = w)
      xy_A <- sf::st_coordinates(cand_sf)[sel, , drop = FALSE]
    }
  }

  # Others as one pool of points; species assigned later by counts
  n_all_others <- length(other_species)
  proc_O <- tolower(P$SPATIAL_PROCESS_OTHERS %||% "poisson")
  xy_O <- switch(
    proc_O,
    "poisson" = .sample_uniform_in_domain(n_all_others, domain),
    "thomas" = .simulate_thomas_points(
      n_all_others,
      domain,
      mu = as.numeric(P$OTHERS_MU %||% 10),
      sigma = as.numeric(P$OTHERS_SIGMA %||% 1),
      kappa = as.numeric(P$OTHERS_BETA %||% NA_real_)
    ), # treat as kappa if given
    "strauss" = .simulate_strauss_points(
      n_all_others,
      domain,
      r = as.numeric(P$OTHERS_R %||% 1),
      s = as.numeric(P$OTHERS_S %||% 0.7)
    ),
    "geyer" = .simulate_geyer_points(
      n_all_others,
      domain,
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
