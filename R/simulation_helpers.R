# Internal simulation helper functions, extracted from
# generate_heterogeneous_distribution() for clarity and maintenance.

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Build a fast point-in-polygon tester by pre-extracting ring coordinates from
# the sf polygon once.  Returns a closure that accepts an nx2 matrix of (x, y)
# points and returns a logical vector.  Use this in tight loops where `poly`
# is fixed to avoid repeated sf/GEOS overhead.
.build_pip_tester <- function(poly) {
  coords <- sf::st_coordinates(poly)
  rx <- coords[, "X"]
  ry <- coords[, "Y"]
  function(xy) {
    if (NROW(xy) == 0L) return(logical(0L))
    pip_cpp(as.numeric(xy[, 1L]), as.numeric(xy[, 2L]), rx, ry)
  }
}

# Replaces the old sf::st_as_sf + sf::st_within implementation.  Extracts ring
# coordinates from `poly` on each call and delegates to the C++ ray-caster.
# For repeated calls with the same polygon, prefer .build_pip_tester() to
# amortise the coordinate extraction cost.
.in_domain <- function(xy, poly) {
  if (NROW(xy) == 0L) return(logical(0L))
  coords <- sf::st_coordinates(poly)
  pip_cpp(as.numeric(xy[, 1L]), as.numeric(xy[, 2L]),
          coords[, "X"], coords[, "Y"])
}

# Accept an optional pre-built pip_fn to avoid re-extracting ring coords on
# every call when sampling many batches from the same polygon.
.sample_uniform_in_domain <- function(n, poly, max_tries = n * 50,
                                      pip_fn = NULL) {
  if (is.null(pip_fn)) pip_fn <- .build_pip_tester(poly)
  bb <- sf::st_bbox(poly)
  res <- matrix(NA_real_, 0, 2)
  tries <- 0L
  while (nrow(res) < n && tries < max_tries) {
    need <- n - nrow(res)
    cand <- cbind(
      stats::runif(need, bb["xmin"], bb["xmax"]),
      stats::runif(need, bb["ymin"], bb["ymax"])
    )
    keep <- pip_fn(cand)
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
  # Pre-build pip tester once to avoid per-call sf/CRS overhead in the loop.
  pip_fn <- .build_pip_tester(poly)
  out <- matrix(NA_real_, 0, 2)
  tries <- 0L
  max_tries <- 2000L + 50L * n
  while (nrow(out) < n && tries < max_tries) {
    cand <- cbind(
      stats::runif(1, bb["xmin"], bb["xmax"]),
      stats::runif(1, bb["ymin"], bb["ymax"])
    )
    if (pip_fn(cand)) {
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
    out <- rbind(out, .sample_uniform_in_domain(need, poly, pip_fn = pip_fn))
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
  if (NROW(xy) == 0L) return(sf::st_sfc(crs = crs))
  # Vectorised path: much faster than lapply(st_point(...)) for large xy.
  sf::st_as_sf(
    data.frame(x = xy[, 1L], y = xy[, 2L]),
    coords = c("x", "y"), crs = crs
  )$geometry
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

.linear_pos_to_xy <- function(pos, st, poly, pip_fn = NULL) {
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
  if (is.null(pip_fn)) pip_fn <- .build_pip_tester(poly)
  bad <- which(!pip_fn(out))
  if (length(bad) > 0) {
    out[bad, ] <- .sample_uniform_in_domain(length(bad), poly, pip_fn = pip_fn)
  }
  out
}

.constrain_xy_to_linear <- function(xy, st, poly) {
  if (!isTRUE(st$enabled) || NROW(xy) == 0) return(xy)
  p <- .xy_to_linear_pos(xy, st)
  .linear_pos_to_xy(p, st, poly)
}
