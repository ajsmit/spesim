#' Run the full spatial sampling simulation and write outputs
#'
#' @description
#' End‑to‑end runner that:
#' \enumerate{
#'   \item reads a text configuration via \code{\link{load_config}()},
#'   \item optionally reads a separate interactions config via \code{\link{load_interactions}()},
#'   \item simulates a heterogeneous community (environmental filtering, dominant clustering, local interactions),
#'   \item places sampling quadrats using the requested scheme,
#'   \item compiles site × species tables and per‑site environment,
#'   \item saves data tables and a 2×2 figure panel to disk,
#'   \item (optionally) computes an advanced analysis panel, and
#'   \item appends a human‑readable text report with key results.
#' }
#'
#' @details
#' \strong{Configuration file} — See \code{\link{load_config}()} for supported keys. At minimum, the file
#' should specify community size (\code{N_INDIVIDUALS}), number of species (\code{N_SPECIES}), gradient
#' assignments for any responsive species (\code{GRADIENT_SPECIES}, \code{GRADIENT_ASSIGNMENTS}), and
#' sampling options (\code{SAMPLING_SCHEME}, quadrat settings).
#'
#' \strong{Interactions file} — If provided, \code{\link{load_interactions}()} parses a neighborhood
#' radius and either a full matrix or edgelist of interaction coefficients. If omitted or unreadable,
#' interactions default to neutral (radius = 0, matrix of 1s).
#'
#' \strong{Sampling schemes} — The function supports:
#' \itemize{
#'   \item \code{“random”} (non‑overlapping rectangles placed uniformly),
#'   \item \code{“tiled”} (systematic grid cells fully contained in the domain),
#'   \item \code{“systematic”} (regular point lattice with rectangles centered on grid points),
#'   \item \code{“transect”} (parallel lines with evenly spaced quadrats),
#'   \item \code{“voronoi”} (cells with largest inscribed circles that fit the quadrat).
#' }
#'
#' \strong{Side effects and outputs} — Files are written adjacent to \code{output_prefix} with an
#' appended timestamp. Filenames include:
#' \itemize{
#'   \item \code{_abundances.csv}: site × species matrix (rows = quadrats),
#'   \item \code{_environments.csv}: per‑quadrat environment summaries,
#'   \item \code{_quadrat_centroids.csv}: centroids of sampled quadrats,
#'   \item \code{_fig_panel.png}: 2×2 panel (species map + 3 environment overlays),
#'   \item \code{_fig_advanced_panel.png}: optional multi‑plot analysis panel,
#'   \item \code{_report.txt}: appended textual report of settings and results.
#' }
#'
#' A sink is opened to append both standard output and messages to the report; it is
#' safely closed on exit even if errors occur. The function returns \code{NULL} (invisibly).
#'
#' @param init_file Character scalar. Path to a simulation config text file (e.g. \code{“simul_init.txt”}).
#'   Lines are parsed as \code{KEY = value} with \code{#} comments ignored. Supports numerics, logicals,
#'   simple vectors (comma‑separated or \code{c(…)}) and named pairs (\code{name:value}). See
#'   \code{\link{load_config}()} for full specification. Default: \code{“simul_init.txt”}.
#'
#' @param interactions_file \emph{Optional} character scalar. Path to a separate interactions config.
#'   Recognized keys include \code{INTERACTION_RADIUS}, \code{MATRIX_CSV}, \code{EDGELIST_CSV}, and \code{AUTO}.
#'   If \code{NULL} (default) or not found, interactions are disabled (radius = 0, all coefficients = 1).
#'
#' @param output_prefix Character scalar used as the base for all output filenames. A timestamp
#'   (\code{YYYYMMDD_HHMMSS}) is appended automatically to avoid overwriting previous runs.
#'   May include a directory component; it will be created if missing. Default: \code{“simulation_output”}.
#'
#' @return \code{NULL}, invisibly. All primary results are written to disk; key objects are also
#'   constructed internally to build plots and the textual report.
#'
#' @section Reproducibility:
#' The random seed is taken from \code{SEED} in the config (default 77) and set within
#' \code{\link{load_config}()}. For full reproducibility across systems, ensure identical versions
#' of \pkg{sf}, GEOS/PROJ, and random number generator defaults.
#'
#' @section Errors and warnings:
#' Most errors are caught and appended to the report under “ERROR in core simulation”.
#' Geometry warnings like “attribute variables are assumed to be spatially constant”
#' originate from \pkg{sf} and are benign for plotting/aggregation here.
#'
#' @seealso \code{\link{load_config}()}, \code{\link{load_interactions}()},
#'   \code{\link{plot_spatial_sampling}()}, \code{\link{generate_full_report}()},
#'   \code{\link{place_quadrats}}, \code{\link{place_quadrats_tiled}},
#'   \code{\link{place_quadrats_systematic}}, \code{\link{place_quadrats_transect}},
#'   \code{\link{place_quadrats_voronoi}}
#'
#' @examples
#' \dontrun{
#' # Minimal run with defaults (reads “simul_init.txt” in working directory)
#' run_spatial_simulation()
#'
#' # Custom configs and output prefix
#' run_spatial_simulation(
#'   init_file = “config/my_sim.txt”,
#'   interactions_file = “config/interactions.txt”,
#'   output_prefix = “out/my_run”
#' )
#' }
#' @export
run_spatial_simulation <- function(init_file = "simul_init.txt",
                                   interactions_file = NULL,
                                   output_prefix = "simulation_output") {
  # --- 1) Setup (paths, report sink) -----------------------------------------
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_prefix <- paste0(output_prefix, "_", timestamp)
  output_dir <- dirname(output_prefix)
  report_path <- paste0(output_prefix, "_report.txt")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  con <- NULL
  on.exit(
    {
      num_sinks <- sink.number()
      if (num_sinks > 0) for (i in seq_len(num_sinks)) sink()
      try(close(con), silent = TRUE)
    },
    add = TRUE
  )

  con <- file(report_path, open = "wt")
  sink(con, type = "output")
  sink(con, type = "message")

  # --- 2) Core simulation block ----------------------------------------------
  results_list <- NULL
  tryCatch(
    {
      P <- load_config(init_file)

      # Load interactions from separate file
      I <- load_interactions(interactions_file, P$N_SPECIES)
      P$INTERACTION_RADIUS <- I$radius
      P$INTERACTION_MATRIX <- I$matrix

      cat("========== SIMULATION PARAMETERS ==========\n")
      cat(paste(capture.output(dput(P[!names(P) %in% "QUADRAT_SIZE"])), collapse = "\n"), "\n\n")

      cat("========== RUNNING SIMULATION ==========\n")
      domain <- create_sampling_domain()
      species_dist <- generate_heterogeneous_distribution(domain, P)

      quadrats <- switch(P$SAMPLING_SCHEME,
        "random"     = place_quadrats(domain, P$N_QUADRATS, P$QUADRAT_SIZE),
        "tiled"      = place_quadrats_tiled(domain, P$N_QUADRATS, P$QUADRAT_SIZE),
        "systematic" = place_quadrats_systematic(domain, P$N_QUADRATS, P$QUADRAT_SIZE),
        "transect"   = place_quadrats_transect(domain, P$N_TRANSECTS, P$N_QUADRATS_PER_TRANSECT, P$QUADRAT_SIZE, P$TRANSECT_ANGLE),
        "voronoi"    = place_quadrats_voronoi(domain, P$N_QUADRATS, P$QUADRAT_SIZE, P$VORONOI_SEED_FACTOR),
        stop("Invalid 'SAMPLING_SCHEME'.")
      )
      if (nrow(quadrats) == 0) stop("Failed to place any quadrats.")

      all_species <- LETTERS[1:P$N_SPECIES]
      env_gradients <- create_environmental_gradients(domain, P$SAMPLING_RESOLUTION, P$ENVIRONMENTAL_NOISE)
      abund_matrix <- create_abundance_matrix(species_dist, quadrats, all_species)
      site_env <- calculate_quadrat_environment(env_gradients, quadrats, sf::st_crs(domain))

      site_coords <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(quadrats))) |>
        as.data.frame() |>
        dplyr::mutate(site = quadrats$quadrat_id) |>
        dplyr::select(site, x = X, y = Y)

      # --- 3) Save outputs ------------------------------------------------------
      cat("\n========== SAVING OUTPUT FILES ==========\n")
      f_abund <- paste0(output_prefix, "_abundances.csv")
      utils::write.csv(abund_matrix, f_abund, row.names = FALSE)

      f_env <- paste0(output_prefix, "_environments.csv")
      utils::write.csv(site_env, f_env, row.names = FALSE)

      f_coord <- paste0(output_prefix, "_quadrat_centroids.csv")
      utils::write.csv(site_coords, f_coord, row.names = FALSE)

      f_png <- paste0(output_prefix, "_fig_panel.png")

      # Build panel plot without patchwork "&" operator
      p1 <- plot_spatial_sampling(domain, species_dist, quadrats, P)
      p2 <- plot_spatial_sampling(domain, species_dist, quadrats, P, TRUE, env_gradients, "temperature_C")
      p3 <- plot_spatial_sampling(domain, species_dist, quadrats, P, TRUE, env_gradients, "elevation_m")
      p4 <- plot_spatial_sampling(domain, species_dist, quadrats, P, TRUE, env_gradients, "rainfall_mm")

      # Apply a consistent legend position (add theme, then tweak legend)
      plots <- list(p1, p2, p3, p4)
      plots <- lapply(plots, function(p) p + ggplot2::theme_get() + ggplot2::theme(legend.position = "right"))

      panel_plot <- patchwork::wrap_plots(plots, ncol = 2) +
        patchwork::plot_layout(guides = "collect")

      ggplot2::ggsave(f_png, panel_plot, width = 14, height = 10, dpi = 300, bg = "white")
      cat(sprintf("Main panel figure saved as %s\n", f_png))
      cat("Data tables saved.\n")

      results_list <- list(
        P = P,
        domain = domain,
        species_dist = species_dist,
        quadrats = quadrats,
        env_gradients = env_gradients,
        abund_matrix = abund_matrix,
        site_coords = site_coords
      )
    },
    error = function(e) {
      cat("\nERROR in core simulation:\n", e$message, "\n", sep = "")
    }
  )

  # --- 4) Advanced analysis (optional) ---------------------------------------
  if (!is.null(results_list) && isTRUE(results_list$P$ADVANCED_ANALYSIS)) {
    tryCatch(
      {
        cat("\n========== ADVANCED ANALYSIS ==========\n")
        f_adv_panel <- paste0(output_prefix, "_fig_advanced_panel.png")
        advanced_panel_plot <- generate_advanced_panel(results_list)
        ggplot2::ggsave(f_adv_panel, advanced_panel_plot, width = 12, height = 14, dpi = 300, bg = "white")
        cat(sprintf("Advanced analysis panel saved as %s\n", f_adv_panel))
      },
      error = function(e) {
        cat("\nERROR during advanced analysis:\n", e$message, "\n", sep = "")
      }
    )
  }

  # --- 5) Final reporting -----------------------------------------------------
  if (!is.null(results_list)) {
    cat("\n", generate_full_report(results_list), "\n", sep = "")
    cat("Outputs saved to: ", normalizePath(output_dir), "\n", sep = "")
  }

  invisible(NULL)
}
