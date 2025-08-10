#' Check for (optional) spatstat availability
#'
#' These simulations no longer require spatstat. This helper remains as a
#' no-op stub to avoid breaking any internal calls that previously relied
#' on it. If you later re-enable spatstat-powered diagnostics, you can
#' repurpose this function to warn or fail accordingly.
#'
#' @keywords internal
.require_spatstat <- function() {
  invisible(TRUE)
}

# ---------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------

#' @keywords internal
.nz <- function(x, y) if (!is.null(x)) x else y

#' Build an sf POINT layer from a numeric matrix of XY
#' @keywords internal
.as_points_sf <- function(xy, crs) {
  if (length(xy) == 0 || NROW(xy) == 0) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = crs)))
  }
  sfc <- sf::st_sfc(lapply(seq_len(NROW(xy)), function(i) sf::st_point(xy[i, 1:2])), crs = crs)
  sf::st_sf(geometry = sfc)
}

#' Keep only points that fall inside the domain polygon
#' @keywords internal
.clip_to_domain <- function(pts_sf, domain) {
  if (nrow(pts_sf) == 0) {
    return(pts_sf)
  }
  inside <- as.logical(sf::st_within(pts_sf, domain, sparse = FALSE)[, 1])
  pts_sf[inside, , drop = FALSE]
}

# ---------------------------------------------------------------------
# Simulators (spatstat-free)
# ---------------------------------------------------------------------

#' Simulate point locations via a Thomas (Gaussian Neyman–Scott) process
#'
#' Parent–offspring clustering without external dependencies:
#' parents are sampled uniformly inside \code{domain}; each parent
#' produces \eqn{\text{Pois}(\mu)} offspring displaced by isotropic
#' Gaussian noise with sd \code{sigma}. Offspring are clipped to the
#' domain. The result is thinned/replicated to match \code{n_target}.
#'
#' @param domain sf polygon/multipolygon.
#' @param n_target Integer target number of points.
#' @param kappa Optional numeric parent intensity (parents per unit area).
#'   If \code{NULL}, a value is derived to hit \code{n_target} given
#'   \code{mu}.
#' @param mu Mean offspring per parent.
#' @param sigma Cluster scale (Gaussian sd, in map units).
#'
#' @return sf POINT layer (no attributes).
#' @keywords internal
simulate_points_thomas <- function(domain, n_target,
                                   kappa = NULL, mu = 10, sigma = 1) {
  crs_dom <- sf::st_crs(domain)
  if (n_target <= 0) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = crs_dom)))
  }

  areaW <- as.numeric(sf::st_area(sf::st_union(domain)))
  if (!is.finite(areaW) || areaW <= 0) stop("simulate_points_thomas: invalid domain area.")

  if (is.null(kappa)) {
    # Roughly n ≈ (kappa * area) * mu  =>  kappa ≈ n / (area * mu)
    kappa <- max(1e-8, n_target / (areaW * max(mu, 1e-8)))
  }

  # Number of parents ~ Poisson(kappa * area)
  n_par <- stats::rpois(1L, lambda = kappa * areaW)
  n_par <- max(1L, n_par)

  parents <- sf::st_sample(domain, size = n_par, type = "random")
  par_xy <- sf::st_coordinates(parents)

  # Offspring per parent ~ Poisson(mu)
  k_off <- stats::rpois(n_par, mu)
  tot_off <- sum(k_off)

  if (tot_off == 0) {
    # fallback: uniform Poisson
    pts <- sf::st_sample(domain, size = n_target, type = "random")
    return(sf::st_sf(geometry = pts))
  }

  # Build offspring coordinates
  # Repeat parents by k_off and add Gaussian noise
  rep_idx <- rep(seq_len(n_par), times = k_off)
  base_xy <- par_xy[rep_idx, , drop = FALSE]
  jitter <- cbind(stats::rnorm(tot_off, sd = sigma), stats::rnorm(tot_off, sd = sigma))
  off_xy <- base_xy + jitter

  pts_sf <- .as_points_sf(off_xy, crs_dom)
  pts_sf <- .clip_to_domain(pts_sf, domain)

  # Harmonize to target size
  n_now <- nrow(pts_sf)
  if (n_now == 0) {
    pts <- sf::st_sample(domain, size = n_target, type = "random")
    return(sf::st_sf(geometry = pts))
  }
  if (n_now > n_target) {
    pts_sf <- pts_sf[sample.int(n_now, n_target), , drop = FALSE]
  } else if (n_now < n_target) {
    add_idx <- sample.int(n_now, n_target - n_now, replace = TRUE)
    pts_sf <- rbind(pts_sf, pts_sf[add_idx, , drop = FALSE])
  }
  pts_sf
}

#' Simulate locations via a Strauss-like (soft inhibition) process
#'
#' Implements a simple sequential acceptance sampler: propose candidate
#' points uniformly inside \code{domain}. Let \code{k} be the number of
#' already-accepted points within radius \code{r}. Accept with probability
#' \eqn{\gamma^{k}} (with \eqn{0 < \gamma \le 1}); otherwise reject.
#' Repeat until \code{n_target} points or a max iteration budget.
#'
#' @param domain sf polygon/multipolygon.
#' @param n_target Integer target number of points.
#' @param beta Unused placeholder (kept for API compatibility).
#' @param gamma Inhibition strength (0–1). Smaller values increase inhibition.
#' @param r Interaction radius (map units).
#'
#' @return sf POINT layer (no attributes).
#' @keywords internal
simulate_points_strauss <- function(domain, n_target,
                                    beta = NULL, gamma = 0.2, r = 1) {
  crs_dom <- sf::st_crs(domain)
  if (n_target <= 0) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = crs_dom)))
  }

  acc <- matrix(NA_real_, 0, 2)
  max_iter <- max(1e4, 200L * n_target)
  iter <- 0L

  while (NROW(acc) < n_target && iter < max_iter) {
    iter <- iter + 1L
    cand <- sf::st_sample(domain, size = 1, type = "random")
    xy <- as.numeric(sf::st_coordinates(cand))[1:2]

    if (NROW(acc) == 0) {
      accept <- TRUE
    } else {
      dx <- acc[, 1] - xy[1]
      dy <- acc[, 2] - xy[2]
      k <- sum(dx * dx + dy * dy <= r * r)
      p_acc <- gamma^k
      accept <- (stats::runif(1) < p_acc)
    }

    if (accept) {
      acc <- rbind(acc, xy)
    }
  }

  .as_points_sf(acc, crs_dom)
}

#' Simulate locations via a Geyer-like saturation process (spatstat-free)
#'
#' Implements a simple sequential acceptance sampler: propose candidate
#' points uniformly inside \code{domain}. Let \code{k} be the number of
#' already-accepted points within radius \code{r}. Accept with probability
#' \eqn{\gamma^{\min(k, s)}}, where \eqn{\gamma} is the interaction parameter
#' and \eqn{s} is the saturation count. For \eqn{\gamma > 1}, the rule induces
#' clustering up to \eqn{s} neighbours (values \eqn{< 1} yield inhibition).
#'
#' @param domain sf polygon/multipolygon.
#' @param n_target Integer target number of points.
#' @param beta Unused placeholder (kept for API compatibility).
#' @param gamma Interaction parameter (numeric). \eqn{\gamma > 1} encourages
#'   clustering; \eqn{\gamma < 1} encourages inhibition.
#' @param r Interaction radius (map units).
#' @param sat Saturation count (non-negative integer).
#'
#' @return sf POINT layer (no attributes).
#' @keywords internal
simulate_points_geyer <- function(domain, n_target,
                                  beta = NULL, gamma = 1.5, r = 1, sat = 2) {
  crs_dom <- sf::st_crs(domain)
  if (n_target <= 0) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = crs_dom)))
  }

  acc <- matrix(NA_real_, 0, 2)
  max_iter <- max(1e4, 200L * n_target)
  iter <- 0L
  sat <- max(0L, floor(sat))

  while (NROW(acc) < n_target && iter < max_iter) {
    iter <- iter + 1L
    cand <- sf::st_sample(domain, size = 1, type = "random")
    xy <- as.numeric(sf::st_coordinates(cand))[1:2]

    if (NROW(acc) == 0) {
      accept <- TRUE
    } else {
      dx <- acc[, 1] - xy[1]
      dy <- acc[, 2] - xy[2]
      k <- sum(dx * dx + dy * dy <= r * r)
      p_acc <- gamma^(min(k, sat))
      # Bound acceptance probability to [0, 1] for numerical safety
      p_acc <- max(0, min(1, p_acc))
      accept <- (stats::runif(1) < p_acc)
    }

    if (accept) {
      acc <- rbind(acc, xy)
    }
  }

  if (NROW(acc) == 0) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = crs_dom)))
  }

  # If we overshot or undershot, thin/replicate to match n_target
  if (NROW(acc) > n_target) {
    acc <- acc[sample.int(NROW(acc), n_target), , drop = FALSE]
  } else if (NROW(acc) < n_target) {
    acc <- rbind(acc, acc[sample.int(NROW(acc), n_target - NROW(acc), replace = TRUE), , drop = FALSE])
  }

  sfc <- sf::st_sfc(lapply(seq_len(NROW(acc)), function(i) sf::st_point(acc[i, 1:2])), crs = crs_dom)
  sf::st_sf(geometry = sfc)
}

# ---------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------

#' Dispatch simulation by process kind
#'
#' A small internal router that calls one of the simulators:
#' \itemize{
#'   \item \code{"poisson"} — homogeneous Poisson (uniform) via \code{sf::st_sample()}
#'   \item \code{"thomas"}  — parent–offspring clustering
#'   \item \code{"strauss"} — soft inhibition (sequential acceptance)
#'   \item \code{"geyer"}   — saturation model (sequential acceptance)
#' }
#'
#' Arguments for a given process can be supplied via \code{args} using either
#' the generic names (\code{beta}, \code{mu}, \code{sigma}, \code{gamma}, \code{r}, \code{sat})
#' or the package’s config-style names (\code{A_PARENT_INTENSITY}, \code{A_MEAN_OFFSPRING},
#' \code{A_CLUSTER_SCALE}, \code{OTHERS_BETA}, \code{OTHERS_GAMMA}, \code{OTHERS_R}, \code{OTHERS_S}).
#'
#' @param kind Character scalar: one of \code{"poisson"}, \code{"thomas"},
#'   \code{"strauss"}, \code{"geyer"} (case-insensitive).
#' @param domain sf polygon/multipolygon sampling window.
#' @param n_target Integer target number of points.
#' @param args Named list of extra parameters, as described above.
#'
#' @return sf POINT layer (no attributes).
#' @keywords internal
simulate_points_dispatch <- function(kind, domain, n_target, args = list()) {
  kind <- tolower(.nz(kind, "poisson"))
  crs_dom <- sf::st_crs(domain)

  use_fast_thomas <- .has_cpp_thomas()
  use_fast_strauss <- .has_cpp_strauss()

  if (n_target <= 0L) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = crs_dom)))
  }

  # Fast path: Poisson via sf::st_sample
  if (kind == "poisson") {
    pts <- sf::st_sample(domain, size = n_target, type = "random")
    return(sf::st_sf(geometry = pts))
  }

  # Helper to build sfc POINT from a 2-col matrix/data.frame
  .xy_to_sfc <- function(xy, crs) {
    if (is.null(xy) || length(xy) == 0L) {
      return(sf::st_sfc(crs = crs))
    }
    if (is.data.frame(xy)) xy <- as.matrix(xy)
    stopifnot(is.matrix(xy), ncol(xy) >= 2)
    sf::st_sfc(lapply(seq_len(nrow(xy)), function(i) sf::st_point(xy[i, 1:2])), crs = crs)
  }

  switch(kind,
    thomas = {
      # Prefer fast Rcpp sampler if available; fall back to spatstat version
      if (isTRUE(exists("rthomas_bbox_cpp", mode = "function"))) {
        bb <- sf::st_bbox(domain)
        mu <- .nz(args$A_MEAN_OFFSPRING, .nz(args$mu, 10))
        sigma <- .nz(args$A_CLUSTER_SCALE, .nz(args$sigma, 1))
        kappa <- .nz(args$A_PARENT_INTENSITY, args$beta) # may be NULL -> pass NA

        # ---- IMPORTANT: positional call (no names) to match Rcpp signature ----
        xy <- rthomas_bbox_cpp(
          as.integer(n_target),
          unname(bb["xmin"]), unname(bb["xmax"]),
          unname(bb["ymin"]), unname(bb["ymax"]),
          as.numeric(mu),
          as.numeric(sigma),
          if (is.null(kappa)) NA_real_ else as.numeric(kappa)
        )
        sfc <- .xy_to_sfc(xy, crs_dom)
        return(sf::st_sf(geometry = sfc))
      } else {
        simulate_points_thomas(
          domain, n_target,
          kappa = .nz(args$A_PARENT_INTENSITY, args$beta),
          mu = .nz(args$A_MEAN_OFFSPRING, .nz(args$mu, 10)),
          sigma = .nz(args$A_CLUSTER_SCALE, .nz(args$sigma, 1))
        )
      }
    },
    strauss = {
      if (.has_cpp_strauss()) {
        simulate_points_strauss_fast(
          domain, n_target,
          r = .nz(args$OTHERS_R, .nz(args$r, 1)),
          gamma = .nz(args$OTHERS_GAMMA, .nz(args$gamma, 0.2))
        )
      } else {
        simulate_points_strauss( # your existing slow fallback (spatstat or R)
          domain, n_target,
          beta = .nz(args$OTHERS_BETA, args$beta),
          gamma = .nz(args$OTHERS_GAMMA, .nz(args$gamma, 0.2)),
          r = .nz(args$OTHERS_R, .nz(args$r, 1))
        )
      }
    },
    geyer = {
      if (kind == "geyer" && .has_cpp_geyer()) {
        message("[spesim] Geyer: using fast Rcpp engine")
        return(simulate_points_geyer_fast(
          domain, n_target,
          r = .nz(args$OTHERS_R, .nz(args$r, 1)),
          gamma = .nz(args$OTHERS_GAMMA, .nz(args$gamma, 1.5)),
          sat = as.integer(.nz(args$OTHERS_S, .nz(args$sat, 2))),
          sweeps = .nz(args$sweeps, 2000),
          burnin = .nz(args$burnin, 200),
          thin = .nz(args$thin, 1)
        ))
      } else {
        # keep your existing spatstat-based or R fallback
        simulate_points_geyer(
          domain, n_target,
          beta = .nz(args$OTHERS_BETA, args$beta),
          gamma = .nz(args$OTHERS_GAMMA, .nz(args$gamma, 1.5)),
          r = .nz(args$OTHERS_R, .nz(args$r, 1)),
          sat = .nz(args$OTHERS_S, .nz(args$sat, 2))
        )
      }
    },
    stop("simulate_points_dispatch: unknown spatial process: ", kind)
  )
}
