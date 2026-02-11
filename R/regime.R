#' Classify whether a spesim run matches the intended qualitative regime
#'
#' @description
#' This helper turns the numeric outputs of [spesim_audit()] into a compact,
#' teaching-friendly judgement (green/amber/red) about whether the realised
#' simulation matches the qualitative regime the user believes they are
#' simulating.
#'
#' The defaults are deliberately conservative: when sample sizes are too small
#' or signals are weak, the result tends toward **amber** rather than over-
#' confident green/red.
#'
#' @param res A `spesim_result` list (as returned by [spesim_run()]).
#' @param audit Optional list returned by [spesim_audit()]. If `NULL`, the
#'   function will compute an audit internally.
#' @param thresholds Optional named list overriding default thresholds. See
#'   Details.
#'
#' @details
#' Default thresholds (you can override via `thresholds=`):
#'
#' - Spatial NN ratio vs CSR (per species; uses dominant A if present):
#'   - `clustered`: nn_ratio < 0.85
#'   - `CSR_like`: 0.85 ≤ nn_ratio ≤ 1.15
#'   - `inhibited`: nn_ratio > 1.15
#'   - minimum points for classification: `min_n = 50`
#'
#' - Filtering (per gradient-responsive species):
#'   - minimum occupied sites: `min_occ = 5`
#'   - `good`: Spearman rho(-|z|, abundance) > 0.2 AND slope(abund~z) < 0
#'   - `opposite`: rho < -0.2 OR slope > 0
#'
#' - Sampling constraint:
#'   - if boundary exclusion estimate > 50% → amber/red (scheme-dependent)
#'   - if random placement attempts are close to the cap (≥90% of max attempts)
#'     → amber
#'
#' @return A named list with sub-lists `spatial`, `filtering`, `sampling`. Each
#' contains:
#' 
#' - `grade`: one of `"green"`, `"amber"`, `"red"`, `"unknown"`
#' - `reason`: short explanation
#' - `metrics`: compact numeric context
#'
#' @export
spesim_regime <- function(res, audit = NULL, thresholds = list()) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (is.null(audit)) {
    audit <- spesim_audit(res)
  }

  thr <- list(
    spatial = list(min_n = 50, clustered = 0.85, inhibited = 1.15),
    filtering = list(min_occ = 5, rho_good = 0.2, rho_bad = -0.2),
    sampling = list(boundary_exclusion_amber = 0.5, attempt_frac_amber = 0.9)
  )
  # shallow override
  for (nm in names(thresholds)) {
    thr[[nm]] <- utils::modifyList(thr[[nm]], thresholds[[nm]])
  }

  # ---- spatial ----------------------------------------------------------------
  sp <- audit$spatial
  spatial_out <- list(grade = "unknown", reason = "no spatial audit", metrics = list())
  if (is.data.frame(sp) && nrow(sp) > 0) {
    # Prefer dominant A if present; else most abundant species.
    sp2 <- sp
    sp2 <- sp2[is.finite(sp2$nn_ratio) & sp2$n >= thr$spatial$min_n, , drop = FALSE]
    if (nrow(sp2) > 0) {
      if ("A" %in% sp2$species) {
        row <- sp2[sp2$species == "A", ][1, , drop = FALSE]
      } else {
        row <- sp2[order(-sp2$n), ][1, , drop = FALSE]
      }
      r <- as.numeric(row$nn_ratio)
      # interpret relative to CSR
      qual <- if (r < thr$spatial$clustered) "clustered" else if (r > thr$spatial$inhibited) "inhibited" else "CSR_like"
      spatial_out$metrics <- list(species = as.character(row$species), n = as.integer(row$n), nn_ratio = r)
      if (qual == "CSR_like") {
        spatial_out$grade <- "green"
        spatial_out$reason <- sprintf("NN ratio suggests CSR-like spacing for %s (%.2f)", row$species, r)
      } else {
        # Not necessarily bad; it depends on what user intended.
        spatial_out$grade <- "amber"
        spatial_out$reason <- sprintf("NN ratio suggests %s for %s (%.2f)", qual, row$species, r)
      }
    } else {
      spatial_out$grade <- "unknown"
      spatial_out$reason <- sprintf("too few points per species (need n >= %d)", thr$spatial$min_n)
    }
  }

  # ---- filtering --------------------------------------------------------------
  f <- audit$filtering
  filtering_out <- list(grade = "unknown", reason = "no filtering audit", metrics = list())
  if (is.data.frame(f) && nrow(f) > 0) {
    # Add occupancy if present, else fall back on n_sites.
    occ <- if ("occupancy_sites" %in% names(f)) f$occupancy_sites else f$n_sites
    ok <- is.finite(f$rho) & is.finite(f$slope) & is.finite(occ) & occ >= thr$filtering$min_occ

    if (!any(ok)) {
      filtering_out$grade <- "unknown"
      filtering_out$reason <- sprintf("too sparse to assess (need occupancy >= %d)", thr$filtering$min_occ)
    } else {
      f2 <- f[ok, , drop = FALSE]
      # strongest signal by |rho|
      f2 <- f2[order(-abs(f2$rho)), , drop = FALSE]
      row <- f2[1, , drop = FALSE]
      rho <- as.numeric(row$rho)
      slope <- as.numeric(row$slope)
      filtering_out$metrics <- list(species = as.character(row$species), gradient = as.character(row$gradient), rho = rho, slope = slope)

      good <- (rho > thr$filtering$rho_good) && (slope < 0)
      opposite <- (rho < thr$filtering$rho_bad) || (slope > 0)

      if (good) {
        filtering_out$grade <- "green"
        filtering_out$reason <- sprintf("abundance peaks near optimum for %s (%s)", row$species, row$gradient)
      } else if (opposite) {
        filtering_out$grade <- "red"
        filtering_out$reason <- sprintf("filtering signal opposite of expected for %s (%s)", row$species, row$gradient)
      } else {
        filtering_out$grade <- "amber"
        filtering_out$reason <- sprintf("weak/unclear filtering signal for %s (%s)", row$species, row$gradient)
      }
    }
  }

  # ---- sampling ---------------------------------------------------------------
  s <- audit$sampling
  sampling_out <- list(grade = "unknown", reason = "no sampling audit", metrics = list())
  if (is.list(s) && length(s) > 0) {
    sampling_out$metrics$scheme <- s$scheme %||% NA_character_

    # Boundary exclusion estimate from candidates (tiled/systematic) OR
    # from explicit rejection counts (random).
    boundary_excl <- NA_real_

    if (!is.null(s$candidates_total) && !is.null(s$candidates_valid) && s$candidates_total > 0) {
      boundary_excl <- 1 - (s$candidates_valid / s$candidates_total)
    } else if (!is.null(s$attempts) && !is.null(s$reject_boundary) && s$attempts > 0) {
      boundary_excl <- s$reject_boundary / s$attempts
    }

    sampling_out$metrics$boundary_exclusion <- boundary_excl

    # Attempts near cap indicates difficulty.
    attempt_frac <- NA_real_
    if (!is.null(s$attempts) && !is.null(s$max_attempts) && s$max_attempts > 0) {
      attempt_frac <- s$attempts / s$max_attempts
    }
    sampling_out$metrics$attempt_frac <- attempt_frac

    if (is.finite(boundary_excl) && boundary_excl > thr$sampling$boundary_exclusion_amber) {
      sampling_out$grade <- "amber"
      sampling_out$reason <- sprintf("high boundary exclusion (~%.0f%%)", 100 * boundary_excl)
    } else if (is.finite(attempt_frac) && attempt_frac >= thr$sampling$attempt_frac_amber) {
      sampling_out$grade <- "amber"
      sampling_out$reason <- sprintf("placement close to attempt cap (%.0f%% of max)", 100 * attempt_frac)
    } else {
      sampling_out$grade <- "green"
      sampling_out$reason <- "sampling placement constraints appear modest"
    }
  }

  list(
    spatial = spatial_out,
    filtering = filtering_out,
    sampling = sampling_out
  )
}
