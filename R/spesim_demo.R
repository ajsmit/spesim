#' Run a compact spesim demonstration workflow
#'
#' @description
#' Convenience helper for teaching and quick exploration.
#' It runs a small simulation (in-memory by default) and optionally returns a
#' set of ready-made plots.
#'
#' The goal is to provide a single, discoverable entry point for the common
#' “basic → intermediate → advanced” workflow described in the workflow vignette.
#'
#' @param level Character scalar. One of `"basic"`, `"intermediate"`, or
#'   `"advanced"`.
#' @param init_file Optional path to an init file. If `NULL`, uses the packaged
#'   example `examples/spesim_init_basic.txt`.
#' @param P Optional parameter list as returned by [load_config()]. If supplied,
#'   it takes precedence over `init_file`.
#' @param seed Optional integer seed. If provided, overrides any seed in `P`.
#'   If `NULL`, the demo uses `P$SEED` when available (so results are reproducible).
#' @param write_outputs Logical. Passed to [spesim_run()]. Default `FALSE` (fast;
#'   no files written).
#' @param output_prefix Optional output prefix passed to [spesim_run()]. Only
#'   used when `write_outputs = TRUE`.
#' @param make_plots Logical. If `TRUE`, returns a list of ggplot/patchwork
#'   objects in `out$plots`.
#' @param ... Additional arguments passed to [spesim_run()].
#'
#' @return A list with elements:
#' 
#' * `res`: the `spesim_result` returned by [spesim_run()].
#' * `plots` (optional): named list of plots (when `make_plots = TRUE`).
#' * `tables` (optional): named list of intermediate analysis tables (when
#'   `level != "basic"`).
#'
#' @examples
#' \dontrun{
#' # Basic (fast): quadrats + individuals
#' x <- spesim_demo("basic", make_plots = TRUE)
#' x$plots$quadrats
#'
#' # Intermediate: adds classic teaching plots
#' y <- spesim_demo("intermediate", make_plots = TRUE)
#' y$plots$rank_abundance
#'
#' # Advanced: adds the full multi-panel diagnostic
#' z <- spesim_demo("advanced", make_plots = TRUE)
#' z$plots$advanced_panel
#' }
#'
#' @export
spesim_demo <- function(
  level = c("basic", "intermediate", "advanced"),
  init_file = NULL,
  P = NULL,
  seed = NULL,
  write_outputs = FALSE,
  output_prefix = NULL,
  make_plots = FALSE,
  ...
) {
  level <- match.arg(level)

  if (is.null(P)) {
    if (is.null(init_file)) {
      init_file <- system.file("examples/spesim_init_basic.txt", package = "spesim")
    }
    P <- load_config(init_file)
  }

  seed_eff <- if (!is.null(seed)) {
    as.integer(seed)
  } else if (!is.null(P$SEED) && is.finite(P$SEED)) {
    as.integer(P$SEED)
  } else {
    NULL
  }

  # Keep demo fast unless explicitly advanced
  if (identical(level, "basic")) {
    P$ADVANCED_ANALYSIS <- FALSE
  }

  res <- spesim_run(
    config = P,
    seed = seed_eff,
    write_outputs = write_outputs,
    output_prefix = output_prefix,
    ...
  )

  out <- list(res = res)

  if (!identical(level, "basic")) {
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

    if (!identical(level, "basic")) {
      plots$rank_abundance <- plot_rank_abundance(out$tables$rank_abundance)
      plots$occupancy_abundance <- plot_occupancy_abundance(out$tables$occupancy_abundance)
      plots$species_area <- plot_species_area(out$tables$species_area)
      plots$distance_decay <- plot_distance_decay(out$tables$distance_decay)
      plots$rarefaction <- plot_rarefaction(out$tables$rarefaction)
    }

    if (identical(level, "advanced")) {
      plots$advanced_panel <- generate_advanced_panel(res)
    }

    out$plots <- plots
  }

  out
}
