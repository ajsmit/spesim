#' Write a standard output bundle for a spesim result
#'
#' @description
#' Internal helper used by [spesim_run()] when `write_outputs = TRUE`.
#' Writes CSV tables, a 2×2 spatial panel, an optional advanced-panel figure,
#' and a human-readable report.
#'
#' This function is separated from the simulation engine so that the engine can
#' be tested without touching the filesystem.
#'
#' @param res A `spesim_result`.
#' @param output_prefix Base path prefix for outputs. A timestamp suffix is
#'   appended. If `NULL`, defaults to `res$P$OUTPUT_PREFIX` or
#'   `"simulation_output"`.
#' @param timestamp POSIXct; used to generate the timestamp suffix.
#' @param include_audit Logical; passed to [generate_full_report()].
#' @param audit_top_n Integer; passed to [generate_full_report()].
#' @param quiet Logical; suppress progress messages.
#'
#' @return A named list of file paths invisibly.
#'
#' @keywords internal
#' @noRd
spesim_write_outputs <- function(res,
                                output_prefix = NULL,
                                timestamp = Sys.time(),
                                include_audit = TRUE,
                                audit_top_n = 6,
                                quiet = FALSE) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (is.null(res) || !is.list(res)) stop("`res` must be a list.")
  if (is.null(res$P)) stop("`res$P` is required.")

  output_prefix <- output_prefix %||% res$P$OUTPUT_PREFIX %||% "simulation_output"
  stamp <- format(timestamp, "%Y%m%d_%H%M%S")
  stamped_prefix <- paste0(output_prefix, "_", stamp)
  output_dir <- dirname(stamped_prefix)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Paths
  f_abund  <- paste0(stamped_prefix, "_abundances.csv")
  f_env    <- paste0(stamped_prefix, "_environments.csv")
  f_coord  <- paste0(stamped_prefix, "_quadrat_centroids.csv")
  f_panel  <- paste0(stamped_prefix, "_fig_panel.png")
  f_adv    <- paste0(stamped_prefix, "_fig_advanced_panel.png")
  f_report <- paste0(stamped_prefix, "_report.txt")

  if (!quiet) message("spesim: writing outputs to ", normalizePath(output_dir))

  # Tables
  utils::write.csv(res$abund_matrix, f_abund, row.names = FALSE)
  utils::write.csv(res$site_env, f_env, row.names = FALSE)
  utils::write.csv(res$site_coords, f_coord, row.names = FALSE)

  # 2x2 spatial panel
  p1 <- plot_spatial_sampling(res$domain, res$species_dist, res$quadrats, res$P)
  env_cols <- setdiff(names(res$env_gradients), c("x", "y"))
  env_cols <- env_cols[vapply(res$env_gradients[env_cols], is.numeric, logical(1))]
  preferred <- c("temperature_C", "elevation_m", "rainfall_mm")
  picks <- unique(c(preferred[preferred %in% env_cols], env_cols))
  picks <- utils::head(picks, 3)
  if (!length(picks)) picks <- "temperature_C"
  while (length(picks) < 3) picks <- c(picks, picks[length(picks)])
  p2 <- plot_spatial_sampling(res$domain, res$species_dist, res$quadrats, res$P, TRUE, res$env_gradients, picks[1])
  p3 <- plot_spatial_sampling(res$domain, res$species_dist, res$quadrats, res$P, TRUE, res$env_gradients, picks[2])
  p4 <- plot_spatial_sampling(res$domain, res$species_dist, res$quadrats, res$P, TRUE, res$env_gradients, picks[3])

  panel_plot <- patchwork::wrap_plots(list(
    p1 + ggplot2::theme(legend.position = "right"),
    p2 + ggplot2::theme(legend.position = "right"),
    p3 + ggplot2::theme(legend.position = "right"),
    p4 + ggplot2::theme(legend.position = "right")
  ), ncol = 2) + patchwork::plot_layout(guides = "collect")

  ggplot2::ggsave(f_panel, panel_plot, width = 14, height = 10, dpi = 300, bg = "white")

  # Advanced panel (optional)
  wrote_adv <- FALSE
  if (isTRUE(res$P$ADVANCED_ANALYSIS)) {
    ap <- generate_advanced_panel(res)
    ggplot2::ggsave(f_adv, ap, width = 12, height = 14, dpi = 300, bg = "white")
    wrote_adv <- TRUE
  }

  # Report text
  report_lines <- c(
    "========== INITIALISING SPATIAL SAMPLING SIMULATION ==========",
    "========== SIMULATION PARAMETERS ==========",
    paste(capture.output(dput(res$P[!names(res$P) %in% "QUADRAT_SIZE"])), collapse = "\n"),
    "",
    generate_full_report(res, include_audit = include_audit, audit_top_n = audit_top_n),
    "",
    paste0("Outputs saved to: ", normalizePath(output_dir))
  )
  writeLines(report_lines, f_report)

  paths <- list(
    abundances_csv = f_abund,
    environments_csv = f_env,
    quadrat_centroids_csv = f_coord,
    panel_png = f_panel,
    advanced_panel_png = if (wrote_adv) f_adv else NULL,
    report_txt = f_report,
    output_dir = output_dir,
    prefix = stamped_prefix
  )

  invisible(paths)
}
