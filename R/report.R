#' Build a human-readable multi-section report from simulation outputs
#'
#' @description
#' Compiles high-level textual summaries of the simulated system and sampling
#' results, including: environmental gradient ranges and responsive species;
#' pairwise correlations among gradients; a species-abundance distribution
#' (with dominant/gradient-responsive tagging); per-quadrat alpha diversity
#' snapshots; diversity partitioning (alpha, beta, gamma; Shannon and Simpson);
#' simple spatial autocorrelation of richness; and a goodness-of-fit check for
#' Fisher's log-series. The report is returned as a single character string with
#' section headers and formatted values suitable for appending to a log file.
#'
#' @param res A named list produced by the simulator with at least the elements:
#' \describe{
#'   \item{\code{P}}{Parameter list returned by \code{\link{load_config}}(). Must
#'     include keys such as \code{N_INDIVIDUALS}, \code{N_SPECIES},
#'     \code{DOMINANT_FRACTION}, \code{FISHER_ALPHA}, \code{FISHER_X}, and, if used,
#'     a tibble \code{GRADIENT} with columns \code{species}, \code{gradient}
#'     (one of "temperature","elevation","rainfall"), \code{optimum}, and \code{tol}.}
#'   \item{\code{env_gradients}}{A data frame with columns \code{temperature_C},
#'     \code{elevation_m}, and \code{rainfall_mm} containing the gridded
#'     environmental fields used for context.}
#'   \item{\code{species_dist}}{An \code{sf} POINT layer of individuals with a
#'     \code{species} column. Used for tallies and alpha snapshots.}
#'   \item{\code{quadrats}}{An \code{sf} POLYGON layer with \code{quadrat_id};
#'     used to compute per-site alpha summaries.}
#'   \item{\code{abund_matrix}}{A site \eqn{\times} species abundance table (first column
#'     \code{site}; remaining columns are species counts), typically returned by
#'     \code{\link{create_abundance_matrix}}().}
#'   \item{\code{site_coords}}{A data frame with columns \code{site}, \code{x},
#'     \code{y} giving quadrat centroids in the same CRS used for analysis.}
#' }
#'
#' @details
#' Sections produced:
#' \itemize{
#'   \item \strong{Environmental Gradients:} min/max/range for temperature (deg C),
#'         elevation (m), rainfall (mm), with a short pattern note and a list of
#'         gradient-responsive species including their optima/tolerances rendered
#'         in natural units.
#'   \item \strong{Gradient Correlations:} pairwise Pearson correlations among
#'         the three gradients and a brief interpretation (orthogonal vs. correlated).
#'   \item \strong{Species Abundance Distribution:} ranked counts with percentage
#'         and labels for dominant and gradient-responsive species.
#'   \item \strong{Spatial Alpha Diversity:} per-quadrat species richness and a
#'         compact species list (with counts) for each quadrat.
#'   \item \strong{Diversity Partitioning:} mean alpha richness (+/-SE), Shannon's
#'         H', Simpson's 1-D, gamma richness, Whittaker and additive beta, mean
#'         pairwise Sorensen dissimilarity (computed as binary Bray-Curtis on
#'         presence-absence), and simple abundance dispersion metrics.
#'   \item \strong{Spatial Autocorrelation:} Pearson correlation between
#'         inter-quadrat distances and differences in richness (a Mantel-style
#'         proxy) with p-value and interpretation.
#'   \item \strong{Fisher Log-series Validation:} RMSE, R^2, max residual between
#'         observed ranks and theoretical log-series abundances; compares configured
#'         \code{FISHER_ALPHA} with \code{\link[vegan]{fisher.alpha}} estimated from
#'         the data.
#'   \item \strong{Computation Notes:} which baseline point-process models were
#'         requested for the dominant and other species (\code{SPATIAL_PROCESS_A},
#'         \code{SPATIAL_PROCESS_OTHERS}), and whether the fast Rcpp engines were
#'         used when applicable:
#'         \itemize{
#'           \item Thomas: \code{rthomas_bbox_cpp} (fast) or spatstat-based fallback
#'           \item Strauss: \code{rstrauss_bbox_cpp} (fast) or spatstat-based fallback
#'         }
#' }
#'
#' Internally, this routine relies on base summaries, \pkg{sf} for spatial
#' intersections, and \pkg{vegan} for diversity indices. It avoids tidy-evaluation
#' in favor of explicit column access to keep dependencies minimal within a
#' non-interactive reporting context.
#'
#' @return A single character scalar containing the full report text. No files
#' are written by this function; callers typically append the string to a log
#' or include it in the simulation report sink.
#'
#' @section Requirements:
#' The \code{res} object must be consistent (all components correspond to the
#' same simulation and coordinate reference system). Missing or empty inputs
#' will lead to sections being populated with \code{NA} summaries or informative
#' defaults where possible; irrecoverable inconsistencies raise errors.
#'
#' @seealso
#' \code{\link{run_spatial_simulation}}, \code{\link{create_abundance_matrix}},
#' \code{\link{calculate_quadrat_environment}}, \code{\link[vegan]{diversity}},
#' \code{\link[vegan]{vegdist}}, \code{\link[vegan]{fisher.alpha}}
#'
#' @examples
#' \dontrun{
#' res <- list(
#'   P = P, env_gradients = env, species_dist = spp,
#'   quadrats = quads, abund_matrix = A, site_coords = C
#' )
#' cat(generate_full_report(res))
#' }
generate_full_report <- function(res) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

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
  .fmt_grad_line <- function(rows, gname, units_label) {
    if (nrow(rows) == 0) {
      return("None")
    }
    paste(sprintf(
      "%s (opt=%.2f [%.1f %s], tol=%.2f [%.1f %s])",
      rows$species, rows$optimum, .opt_to_units(rows$optimum, gname), units_label,
      rows$tol, .tol_to_units(rows$tol, gname), units_label
    ), collapse = ", ")
  }

  env_report <- c("\nEnvironmental Gradients:")
  temp_range <- range(res$env_gradients$temperature_C, na.rm = TRUE)
  elev_range <- range(res$env_gradients$elevation_m, na.rm = TRUE)
  rain_range <- range(res$env_gradients$rainfall_mm, na.rm = TRUE)

  if (!is.null(res$P$GRADIENT)) {
    G <- res$P$GRADIENT
    temp_species_str <- .fmt_grad_line(G[G$gradient == "temperature", , drop = FALSE], "temperature", "deg C")
    elev_species_str <- .fmt_grad_line(G[G$gradient == "elevation", , drop = FALSE], "elevation", "m")
    rain_species_str <- .fmt_grad_line(G[G$gradient == "rainfall", , drop = FALSE], "rainfall", "mm")
  } else {
    temp_species_str <- "None"
    elev_species_str <- "None"
    rain_species_str <- "None"
  }

  env_report <- c(
    env_report,
    sprintf("  Temperature: %.1f-%.1f deg C (range: %.1f deg C)", temp_range[1], temp_range[2], diff(temp_range)),
    "    Pattern: Diagonal (NW cool -> SE warm)",
    paste0("    Responsive species: ", temp_species_str),
    sprintf("  Elevation: %.0f-%.0f m (range: %.0f m)", elev_range[1], elev_range[2], diff(elev_range)),
    "    Pattern: Central peak (mountain-like topology)",
    paste0("    Responsive species: ", elev_species_str),
    sprintf("  Rainfall: %.0f-%.0f mm (range: %.0f mm)", rain_range[1], rain_range[2], diff(rain_range)),
    "    Pattern: Perpendicular (NE dry -> SW wet)",
    paste0("    Responsive species: ", rain_species_str)
  )

  cor_mat <- cor(res$env_gradients[, c("temperature_C", "elevation_m", "rainfall_mm")], use = "complete.obs")
  interp <- if (max(abs(cor_mat[upper.tri(cor_mat)])) < 0.3) "Gradients are approximately orthogonal (low correlation)" else "Some gradient correlation detected"
  corr_report <- c(
    "\nGradient Correlations:",
    sprintf("  Temperature-Elevation: r=%.3f", cor_mat[1, 2]),
    sprintf("  Temperature-Rainfall:  r=%.3f", cor_mat[1, 3]),
    sprintf("  Elevation-Rainfall:    r=%.3f", cor_mat[2, 3]),
    paste0("  Interpretation: ", interp)
  )

  sad <- as.data.frame(table(res$species_dist$species))
  colnames(sad) <- c("Species", "Count")
  sad <- dplyr::arrange(sad, dplyr::desc(Count))
  sad <- dplyr::mutate(
    sad,
    Percent = 100 * Count / sum(Count),
    Role = dplyr::case_when(
      !is.null(res$P$GRADIENT) & Species %in% res$P$GRADIENT$species ~ {
        g <- res$P$GRADIENT$gradient[match(Species, res$P$GRADIENT$species)]
        sprintf("[%s-RESPONSIVE]", toupper(g))
      },
      Species == "A" ~ "[DOMINANT - clustered]",
      TRUE ~ "[SUBORDINATE]"
    )
  )
  sad_report <- c(
    "\nSpecies Abundance Distribution:",
    sprintf("  Species %s: %3d individuals (%5.1f%%) %s", sad$Species, sad$Count, sad$Percent, sad$Role)
  )

  alpha_report <- c("\nSpatial Distribution of Alpha Diversity:")
  for (i in 1:nrow(res$quadrats)) {
    spp_in_q <- suppressWarnings(sf::st_intersection(res$species_dist, res$quadrats[i, ]))
    if (nrow(spp_in_q) > 0) {
      abunds <- as.data.frame(table(spp_in_q$species))
      abunds_str <- paste(sprintf("%s(%s)", abunds$Var1, abunds$Freq), collapse = ", ")
      alpha_report <- c(
        alpha_report,
        sprintf("  Quadrat %2d: alpha = %2d species | N = %3d individuals", i, dplyr::n_distinct(spp_in_q$species), nrow(spp_in_q)),
        sprintf("    Species: %s", abunds_str)
      )
    } else {
      alpha_report <- c(
        alpha_report,
        sprintf("  Quadrat %2d: alpha = %2d species | N = %3d individuals", i, 0, 0),
        "    Species: None"
      )
    }
  }

  abund_data <- res$abund_matrix[, setdiff(colnames(res$abund_matrix), "site"), drop = FALSE]
  richness_data <- data.frame(richness = rowSums(abund_data > 0), n_ind = rowSums(abund_data))
  shannon_H <- vegan::diversity(abund_data, index = "shannon")
  simpson_D <- vegan::diversity(abund_data, index = "simpson")
  mean_alpha <- mean(richness_data$richness)
  se_alpha <- stats::sd(richness_data$richness) / sqrt(nrow(richness_data))
  gamma_div <- dplyr::n_distinct(res$species_dist$species)
  beta_whittaker <- gamma_div / mean_alpha
  beta_additive <- gamma_div - mean_alpha
  pa_matrix <- (as.matrix(abund_data) > 0) * 1
  mean_sorensen <- mean(vegan::vegdist(pa_matrix, method = "bray", binary = TRUE))

  div_report <- c(
    "\nDiversity Partitioning:",
    sprintf("Alpha (mean local richness): %.2f +/- %.2f SE", mean_alpha, se_alpha),
    sprintf("Shannon's H' (mean): %.3f +/- %.3f SE", mean(shannon_H), stats::sd(shannon_H) / sqrt(length(shannon_H))),
    sprintf("Simpson's (1-D, mean): %.3f +/- %.3f SE", mean(simpson_D), stats::sd(simpson_D) / sqrt(length(simpson_D))),
    sprintf("Gamma (regional species pool): %d species", gamma_div),
    sprintf("Beta (Whittaker): %.2f", beta_whittaker),
    sprintf("Beta (additive): %.2f", beta_additive),
    sprintf("Mean pairwise beta (Sorensen): %.3f", mean_sorensen),
    sprintf("Mean quadrat abundance: %.1f +/- %.1f", mean(richness_data$n_ind), stats::sd(richness_data$n_ind)),
    sprintf("Abundance variation (CV): %.3f", stats::sd(richness_data$n_ind) / mean(richness_data$n_ind))
  )

  ## --- Robust Mantel-style summary ------------------------------------------
  space_dist <- stats::dist(res$site_coords[, c("x", "y")])
  richness_dist <- stats::dist(richness_data$richness)

  mantel_interp <- "  No significant spatial autocorrelation."
  mantel_r <- NA_real_
  mantel_p <- NA_real_

  if (nrow(res$site_coords) >= 4 && var(richness_data$richness) > 0) {
    mt <- suppressWarnings(stats::cor.test(space_dist, richness_dist, method = "pearson"))
    mantel_r <- as.numeric(mt$estimate)
    mantel_p <- as.numeric(mt$p.value)
    if (is.finite(mantel_p) && mantel_p < 0.05) {
      mantel_interp <- if (is.finite(mantel_r) && mantel_r > 0) {
        "  Significant positive autocorrelation."
      } else {
        "  Significant negative autocorrelation."
      }
    }
  }

  mantel_r_str <- if (is.finite(mantel_r)) {
    sprintf("Spatial autocorrelation in richness (Mantel r): %.3f", mantel_r)
  } else {
    "Spatial autocorrelation in richness (Mantel r): NA"
  }

  mantel_p_str <- if (is.finite(mantel_p)) {
    sprintf("Statistical significance: p = %.3f", mantel_p)
  } else {
    "Statistical significance: p = NA"
  }

  spat_report <- c(
    "\nSpatial Autocorrelation Analysis:",
    mantel_r_str,
    mantel_p_str,
    mantel_interp,
    sprintf("  Mean inter-quadrat distance: %.2f", mean(space_dist)),
    sprintf("  Richness variance: %.2f", var(richness_data$richness))
  )
  ## --------------------------------------------------------------------------

  obs_abund <- sort(table(res$species_dist$species), decreasing = TRUE)
  n_rem <- res$P$N_INDIVIDUALS * (1 - res$P$DOMINANT_FRACTION)
  ranks <- 2:res$P$N_SPECIES
  rel_abund_theory <- res$P$FISHER_ALPHA * (res$P$FISHER_X^ranks) / ranks
  theory_abund <- c(res$P$N_INDIVIDUALS * res$P$DOMINANT_FRACTION, rel_abund_theory / sum(rel_abund_theory) * n_rem)
  len_obs <- length(obs_abund)
  len_theory <- length(theory_abund)
  max_len <- max(len_obs, len_theory)
  length(obs_abund) <- max_len
  length(theory_abund) <- max_len
  obs_abund[is.na(obs_abund)] <- 0
  residuals <- obs_abund - theory_abund
  rmse <- sqrt(mean(residuals^2, na.rm = TRUE))
  r_squared <- if (var(obs_abund, na.rm = TRUE) > 0) stats::cor(obs_abund, theory_abund, use = "complete.obs")^2 else 1
  eff_alpha <- vegan::fisher.alpha(table(res$species_dist$species))

  fisher_report <- c(
    "\nFisher's Log Series Model Validation:",
    sprintf("  RMSE: %.2f", rmse),
    sprintf("  R-squared: %.3f", r_squared),
    sprintf("  Max residual: %.1f", max(abs(residuals), na.rm = TRUE)),
    sprintf("  Specified alpha: %.2f", res$P$FISHER_ALPHA),
    sprintf("  Effective alpha from data: %.2f", eff_alpha)
  )

  ## --- Computation Notes -----------------------------------------------------
  proc_A <- tolower(as.character(res$P$SPATIAL_PROCESS_A %||% "poisson"))
  proc_O <- tolower(as.character(res$P$SPATIAL_PROCESS_OTHERS %||% "poisson"))
  fast_thomas_available <- .has_cpp_thomas()
  fast_strauss_available <- .has_cpp_strauss()
  fast_geyer_available <- .has_cpp_geyer()

  note_A <- sprintf(
    "  Dominant (A) process: %s%s",
    proc_A,
    if (identical(proc_A, "thomas")) {
      if (fast_thomas_available) " -- fast Rcpp engine detected" else " -- spatstat-based fallback"
    } else if (identical(proc_A, "strauss")) {
      if (fast_strauss_available) " -- fast Rcpp engine detected" else " -- spatstat-based fallback"
    } else {
      ""
    }
  )

  note_O <- sprintf(
    "  Others process: %s%s",
    proc_O,
    if (identical(proc_O, "thomas")) {
      if (fast_thomas_available) " -- fast Rcpp engine detected" else " -- spatstat-based fallback"
    } else if (identical(proc_O, "strauss")) {
      if (fast_strauss_available) " -- fast Rcpp engine detected" else " -- spatstat-based fallback"
    } else if (identical(proc_O, "geyer")) {
      if (fast_geyer_available) " -- fast Rcpp engine detected" else " -- spatstat-based fallback"
    } else {
      ""
    }
  )

  comp_notes <- c("\nComputation Notes:", note_A, note_O)
  ## --------------------------------------------------------------------------

  full_report <- c(
    "========== ANALYSIS REPORT ==========",
    env_report,
    corr_report,
    sad_report,
    alpha_report,
    div_report,
    spat_report,
    fisher_report,
    comp_notes,
    "\nSIMULATION COMPLETED SUCCESSFULLY."
  )
  paste(full_report, collapse = "\n")
}
