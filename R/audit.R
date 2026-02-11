#' Audit whether a spesim run matches the intended qualitative regime
#'
#' @description
#' spesim is designed for teaching and method testing: you often want to know
#' whether the *realised* simulation actually exhibits the qualitative regime you
#' think you requested (e.g., clustering vs inhibition; strong vs weak
#' environmental filtering; heavy boundary exclusion in a sampling scheme).
#'
#' `spesim_audit()` computes lightweight, dependency-minimal diagnostics that can
#' be used as a conceptual "tutor": it reports what the simulator instantiated,
#' and quantifies key patterns.
#'
#' - **Spatial structure:** nearest-neighbour (NN) diagnostics per species,
#'   comparing the observed mean NN distance to the CSR expectation.
#' - **Environmental filtering:** checks whether quadrat abundances for
#'   gradient-responsive species peak near their configured optima.
#' - **Sampling design:** reports boundary exclusion / rejection rates for the
#'   chosen quadrat placement scheme (when available).
#'
#' For more formal point-process diagnostics (Ripley's K, pair correlation
#' functions, edge corrections in irregular domains), use `spatstat` tooling on
#' the exported point patterns.
#'
#' @param res A `spesim_result` list returned by [spesim_run()] or
#'   [run_spatial_simulation()].
#' @param species Character vector of species labels to audit for spatial
#'   structure. Default `"all"` audits all species present.
#' @param nn_k Integer; the neighbour order used for NN distances. Default 1.
#' @param diagnostics Character vector controlling which diagnostics to compute.
#'   Default `c("nn", "spatstat")` will try to compute both nearest-neighbour
#'   and spatstat-based diagnostics; if spatstat is not installed, those
#'   diagnostics are skipped with a message. Use `"nn"` only for a lightweight
#'   audit.
#'
#' @return A list of class `spesim_audit` with components:
#'
#' - `spatial`: data.frame of NN diagnostics per species
#' - `spatial_spatstat`: data.frame of edge-corrected diagnostics (if requested)
#' - `filtering`: data.frame of environment-abundance checks (may be empty)
#' - `sampling`: a list describing quadrat placement rejection/exclusion metrics
#' - `regime`: green/amber/red classification (see [spesim_regime()])
#'
#' @export
spesim_audit <- function(res, species = "all", nn_k = 1, diagnostics = c("nn", "spatstat")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (is.null(res$species_dist) || is.null(res$P) || is.null(res$domain)) {
    stop("`res` must contain at least: P, domain, species_dist")
  }

  diagnostics <- unique(as.character(diagnostics))

  spatial <- if ("nn" %in% diagnostics) {
    audit_spatial_structure(res$species_dist, res$domain, species = species, nn_k = nn_k)
  } else {
    data.frame()
  }

  spatial_spatstat <- data.frame()
  if ("spatstat" %in% diagnostics) {
    if (requireNamespace("spatstat.geom", quietly = TRUE) && requireNamespace("spatstat.explore", quietly = TRUE)) {
      spatial_spatstat <- spesim_spatstat_diagnostics(res$species_dist, res$domain, species = species)
    } else {
      message("spesim_audit: spatstat diagnostics requested but spatstat is not installed; skipping.")
    }
  }

  filtering <- audit_environmental_filtering(res)
  sampling <- audit_sampling_scheme(res$quadrats %||% NULL)

  out <- list(
    spatial = spatial,
    spatial_spatstat = spatial_spatstat,
    filtering = filtering,
    sampling = sampling
  )
  out$regime <- spesim_regime(res, audit = out)

  class(out) <- "spesim_audit"
  out
}

#' @export
print.spesim_audit <- function(x, ...) {
  cat("<spesim_audit>\n")

  if (!is.null(x$regime) && is.list(x$regime)) {
    cat("\nRegime classification (teaching-friendly):\n")
    # print in a compact way
    for (nm in names(x$regime)) {
      if (nm %in% c("spatial", "filtering", "sampling")) {
        r <- x$regime[[nm]]
        if (is.list(r) && !is.null(r$grade)) {
          cat(sprintf("  %s: %s", nm, r$grade))
          if (!is.null(r$reason)) cat(sprintf(" - %s", r$reason))
          cat("\n")
        }
      }
    }
  }

  cat("\nSpatial structure (nearest-neighbour diagnostics):\n")
  if (nrow(x$spatial) == 0) {
    cat("  (no NN spatial diagnostics available)\n")
  } else {
    keep <- c("species", "n", "nn_mean", "nn_expected_csr", "nn_ratio", "qualitative")
    keep <- keep[keep %in% names(x$spatial)]
    print(x$spatial[, keep, drop = FALSE], row.names = FALSE)
  }

  cat("\nSpatial structure (spatstat; edge-corrected):\n")
  if (is.null(x$spatial_spatstat) || nrow(x$spatial_spatstat) == 0) {
    cat("  (no spatstat diagnostics available)\n")
  } else {
    keep <- c("species", "n", "Lmr_mean", "g_max", "g_min", "qualitative")
    keep <- keep[keep %in% names(x$spatial_spatstat)]
    print(x$spatial_spatstat[, keep, drop = FALSE], row.names = FALSE)
  }

  cat("\nEnvironmental filtering checks:\n")
  if (nrow(x$filtering) == 0) {
    cat("  (no gradient-responsive species configured)\n")
  } else {
    keep <- c("species", "gradient", "n_sites", "occupancy_sites", "rho", "slope", "qualitative")
    keep <- keep[keep %in% names(x$filtering)]
    print(x$filtering[, keep, drop = FALSE], row.names = FALSE)
  }

  cat("\nSampling scheme diagnostics:\n")
  if (length(x$sampling) == 0) {
    cat("  (no sampling audit attached to quadrats)\n")
  } else {
    cat(sprintf("  scheme: %s\n", x$sampling$scheme %||% "?"))
    if (!is.null(x$sampling$summary_lines)) {
      for (ln in x$sampling$summary_lines) cat("  ", ln, "\n", sep = "")
    }
  }

  invisible(x)
}

# ---- internals (exported for reuse, but still lightweight) -------------------

#' Audit spatial structure via nearest-neighbour distances
#'
#' @param species_dist `sf` points with a `species` column.
#' @param domain `sf` polygon used only for area estimation.
#' @param species Species labels to include; `"all"` uses all present.
#' @param nn_k neighbour order (1 = nearest neighbour).
#'
#' @return data.frame with NN summary per species.
#' @export
audit_spatial_structure <- function(species_dist, domain, species = "all", nn_k = 1) {
  if (!inherits(species_dist, "sf")) stop("`species_dist` must be an sf object")
  if (!"species" %in% names(species_dist)) stop("`species_dist` must have a `species` column")

  spp_all <- sort(unique(as.character(species_dist$species)))
  spp <- if (identical(species, "all")) spp_all else intersect(as.character(species), spp_all)

  if (length(spp) == 0) {
    return(data.frame())
  }

  # Area estimate (irregular polygon); used for CSR expectation.
  area <- as.numeric(sf::st_area(sf::st_union(domain)))
  if (!is.finite(area) || area <= 0) area <- NA_real_

  out <- lapply(spp, function(s) {
    pts <- species_dist[species_dist$species == s, ]
    n <- nrow(pts)
    if (n < (nn_k + 1)) {
      return(data.frame(
        species = s, n = n,
        nn_mean = NA_real_, nn_median = NA_real_, nn_cv = NA_real_,
        nn_expected_csr = NA_real_, nn_ratio = NA_real_,
        qualitative = "too_few_points"
      ))
    }

    xy <- sf::st_coordinates(pts)

    # Use RANN (already an import) for kNN distances.
    # nn2 returns self as first neighbour (distance 0), so we take nn_k+1.
    knn <- RANN::nn2(data = xy, query = xy, k = nn_k + 1)
    d <- knn$nn.dists[, nn_k + 1]

    nn_mean <- mean(d)
    nn_median <- stats::median(d)
    nn_cv <- stats::sd(d) / nn_mean

    # CSR expectation for mean nearest-neighbour distance in 2D Poisson:
    # E[r] = 1 / (2 * sqrt(lambda)), with lambda = n / area.
    if (is.finite(area)) {
      lambda <- n / area
      nn_expected <- 1 / (2 * sqrt(lambda))
      ratio <- nn_mean / nn_expected
    } else {
      nn_expected <- NA_real_
      ratio <- NA_real_
    }

    qualitative <- NA_character_
    if (is.finite(ratio)) {
      qualitative <- if (ratio < 0.9) "more_clustered_than_CSR" else if (ratio > 1.1) "more_inhibited_than_CSR" else "CSR_like"
    }

    data.frame(
      species = s,
      n = n,
      nn_mean = nn_mean,
      nn_median = nn_median,
      nn_cv = nn_cv,
      nn_expected_csr = nn_expected,
      nn_ratio = ratio,
      qualitative = qualitative
    )
  })

  do.call(rbind, out)
}

#' Audit environmental filtering (optima/tolerances vs realised abundances)
#'
#' @param res Simulation results list.
#'
#' @return A data.frame; one row per gradient-responsive species, with columns:
#' \describe{
#'   \item{species}{Species label.}
#'   \item{gradient}{Gradient name (e.g. temperature/elevation/rainfall).}
#'   \item{n_sites}{Number of sampled quadrats with finite environment + abundance.}
#'   \item{occupancy_sites}{Number of quadrats with abundance > 0 (a simple sparsity check).}
#'   \item{rho}{Spearman correlation between abundance and closeness-to-optimum (higher is more consistent with filtering).}
#'   \item{slope}{Slope from a simple linear model of abundance vs distance-from-optimum (sign is diagnostic only).}
#'   \item{qualitative}{Teaching-friendly label (e.g. abundance_peaks_near_optimum, weak_or_no_signal, too_sparse_to_assess).}
#' }
#' @export
audit_environmental_filtering <- function(res) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (is.null(res$P$GRADIENT) || nrow(res$P$GRADIENT) == 0) {
    return(data.frame())
  }
  if (is.null(res$abund_matrix) || is.null(res$quadrats) || is.null(res$env_gradients)) {
    return(data.frame())
  }

  # compute per-quadrat environment (natural units)
  E <- calculate_quadrat_environment(res$env_gradients, res$quadrats, sf::st_crs(res$domain))

  # abundance table: ensure site keyed
  A <- res$abund_matrix

  # helper: convert optimum/tol (0..1) to natural units to match E
  .opt_to_units <- function(opt, gname) {
    switch(gname,
      "temperature" = opt * 30 - 2,
      "elevation"   = opt * 2000,
      "rainfall"    = opt * 700 + 200,
      NA_real_
    )
  }
  .tol_to_units <- function(tol, gname) {
    switch(gname,
      "temperature" = tol * 30,
      "elevation"   = tol * 2000,
      "rainfall"    = tol * 700,
      NA_real_
    )
  }
  .env_col <- function(gname) {
    switch(gname,
      "temperature" = "temperature_C",
      "elevation"   = "elevation_m",
      "rainfall"    = "rainfall_mm",
      NA_character_
    )
  }

  G <- res$P$GRADIENT
  rows <- lapply(seq_len(nrow(G)), function(i) {
    s <- as.character(G$species[i])
    g <- as.character(G$gradient[i])
    env_col <- .env_col(g)
    if (is.na(env_col) || !(s %in% colnames(A))) {
      return(NULL)
    }

    opt_u <- .opt_to_units(as.numeric(G$optimum[i]), g)
    tol_u <- .tol_to_units(as.numeric(G$tol[i]), g)
    if (!is.finite(tol_u) || tol_u <= 0) tol_u <- NA_real_

    df <- dplyr::left_join(E, A, by = c("site" = "site"))
    env <- df[[env_col]]
    abund <- df[[s]]

    ok <- is.finite(env) & is.finite(abund)
    env <- env[ok]; abund <- abund[ok]

    occupancy_sites <- sum(abund > 0, na.rm = TRUE)

    if (length(env) < 5 || occupancy_sites < 2) {
      return(data.frame(
        species = s,
        gradient = g,
        n_sites = length(env),
        occupancy_sites = occupancy_sites,
        rho = NA_real_,
        slope = NA_real_,
        qualitative = if (length(env) < 5) "too_few_sites" else "too_sparse_to_assess"
      ))
    }

    # Standardise distance from optimum by tolerance where possible.
    if (is.finite(opt_u) && is.finite(tol_u)) {
      z <- abs(env - opt_u) / tol_u
    } else if (is.finite(opt_u)) {
      z <- abs(env - opt_u)
    } else {
      # if no optimum in units, we can only report correlation with env
      z <- env
    }

    # Expectation under filtering: abundance decreases as distance-from-optimum increases.
    # We report Spearman rho (robust) and a simple linear slope.
    rho <- suppressWarnings(stats::cor(-z, abund, method = "spearman"))
    fit <- try(stats::lm(abund ~ z), silent = TRUE)
    slope <- if (inherits(fit, "try-error")) NA_real_ else as.numeric(stats::coef(fit)["z"])

    qualitative <- NA_character_
    if (is.finite(rho)) {
      qualitative <- if (occupancy_sites < 5) {
        "too_sparse_to_assess"
      } else if (rho > 0.2) {
        "abundance_peaks_near_optimum"
      } else if (rho < -0.2) {
        "opposite_of_expected"
      } else {
        "weak_or_no_signal"
      }
    }

    data.frame(
      species = s,
      gradient = g,
      n_sites = length(env),
      occupancy_sites = occupancy_sites,
      rho = rho,
      slope = slope,
      qualitative = qualitative
    )
  })

  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

#' Audit quadrat placement: boundary exclusion and rejection metrics
#'
#' @param quadrats An `sf` object returned by a `place_quadrats_*()` function.
#'
#' @return A list; empty if no audit info is attached.
#' @export
audit_sampling_scheme <- function(quadrats) {
  if (is.null(quadrats)) return(list())
  aud <- .validate_placement_audit(quadrats, stop_on_error = FALSE)
  if (length(aud) == 0) return(list())

  # Provide a couple of compact summary lines for printing.
  summary_lines <- character(0)
  if (!is.null(aud$n_requested) && !is.null(aud$n_returned)) {
    summary_lines <- c(summary_lines, sprintf("requested: %s | returned: %s", aud$n_requested, aud$n_returned))
  }
  if (!is.null(aud$attempts)) {
    summary_lines <- c(summary_lines, sprintf("attempts: %s", aud$attempts))
  }
  if (!is.null(aud$reject_boundary)) {
    summary_lines <- c(summary_lines, sprintf("rejected (boundary): %s", aud$reject_boundary))
  }
  if (!is.null(aud$reject_overlap)) {
    summary_lines <- c(summary_lines, sprintf("rejected (overlap): %s", aud$reject_overlap))
  }
  if (!is.null(aud$candidates_total) && !is.null(aud$candidates_valid)) {
    rate <- 1 - (aud$candidates_valid / aud$candidates_total)
    summary_lines <- c(summary_lines, sprintf("boundary exclusion: %.1f%% (%s/%s valid)", 100 * rate, aud$candidates_valid, aud$candidates_total))
  }

  if (!is.null(aud$safe_area_fraction) && is.finite(aud$safe_area_fraction)) {
    summary_lines <- c(summary_lines, sprintf("effective sampling frame (area fraction): %.1f%%", 100 * aud$safe_area_fraction))
  }

  aud$summary_lines <- summary_lines
  aud
}
