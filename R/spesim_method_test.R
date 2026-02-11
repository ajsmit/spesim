#' Run a method-testing bundle (simulate + audit + standard outputs)
#'
#' @description
#' Convenience helper for teaching and **method testing**.
#'
#' `spesim_method_test()` runs a simulation, computes a lightweight audit of the
#' realised regime (via [spesim_audit()]), and optionally returns a small set of
#' standard analysis tables and plots.
#'
#' Compared to [spesim_demo()], this helper is oriented toward *auditing and
#' comparing settings* (e.g., sampling schemes) rather than a curated "basic →
#' advanced" teaching progression.
#'
#' @param init_file Optional path to an init file. If `NULL`, uses the packaged
#'   example `examples/spesim_init_basic.txt`.
#' @param P Optional parameter list as returned by [load_config()]. If supplied,
#'   it takes precedence over `init_file`.
#' @param seed Optional integer seed. If provided, overrides any seed in `P`.
#' @param write_outputs Logical. Passed to [run_spatial_simulation()]. Default
#'   `FALSE` (fast; no files written).
#' @param output_prefix Optional output prefix passed to
#'   [run_spatial_simulation()]. Only used when `write_outputs = TRUE`.
#' @param diagnostics Diagnostics to compute in [spesim_audit()]. Default `"nn"`
#'   for deterministic, dependency-minimal behaviour.
#' @param make_tables Logical; if `TRUE`, include standard analysis tables.
#' @param make_plots Logical; if `TRUE`, include standard plots.
#' @param ... Additional arguments passed to [run_spatial_simulation()].
#'
#' @return A list with elements:
#' * `res`: results list returned by [run_spatial_simulation()].
#' * `audit`: audit list returned by [spesim_audit()].
#' * `tables` (optional): named list of analysis tables (when `make_tables = TRUE`).
#' * `plots` (optional): named list of plots (when `make_plots = TRUE`).
#'
#' @examples
#' \dontrun{
#' x <- spesim_method_test(make_plots = TRUE)
#' x$audit
#' x$plots$spatial
#' }
#'
#' @export
spesim_method_test <- function(
  init_file = NULL,
  P = NULL,
  seed = NULL,
  write_outputs = FALSE,
  output_prefix = NULL,
  diagnostics = "nn",
  make_tables = TRUE,
  make_plots = TRUE,
  ...
) {
  if (is.null(P)) {
    if (is.null(init_file)) {
      init_file <- system.file("examples/spesim_init_basic.txt", package = "spesim")
    }
    P <- load_config(init_file)
  }

  if (!is.null(seed)) {
    P$SEED <- as.integer(seed)
    set.seed(P$SEED)
  }

  res <- run_spatial_simulation(
    P = P,
    write_outputs = write_outputs,
    output_prefix = output_prefix,
    ...
  )

  out <- list(
    res = res,
    audit = spesim_audit(res, diagnostics = diagnostics)
  )

  if (isTRUE(make_tables)) {
    out$tables <- list(
      rank_abundance = calculate_rank_abundance(res$species_dist, res$P),
      occupancy_abundance = calculate_occupancy_abundance(res$abund_matrix),
      species_area = calculate_species_area(res$abund_matrix),
      distance_decay = calculate_distance_decay(res$abund_matrix, res$site_coords),
      rarefaction = calculate_rarefaction(res$abund_matrix)
    )
  }

  if (isTRUE(make_plots)) {
    plots <- list(
      quadrats = plot_quadrats(res$domain, res$quadrats, title = "Quadrats"),
      spatial = plot_spatial_sampling(res$domain, res$species_dist, res$quadrats, res$P)
    )

    # Only add the classic plots if we computed their tables
    if (!is.null(out$tables)) {
      plots$rank_abundance <- plot_rank_abundance(out$tables$rank_abundance)
      plots$occupancy_abundance <- plot_occupancy_abundance(out$tables$occupancy_abundance)
      plots$species_area <- plot_species_area(out$tables$species_area)
      plots$distance_decay <- plot_distance_decay(out$tables$distance_decay)
      plots$rarefaction <- plot_rarefaction(out$tables$rarefaction)
    }

    # If gradient species exist, add the filtering sanity plot too
    if (!is.null(res$P$GRADIENT) && nrow(res$P$GRADIENT) > 0) {
      plots$filtering_response <- plot_filtering_response(res)
    }

    out$plots <- plots
  }

  out
}
