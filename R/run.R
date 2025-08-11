#' Run the full spatial sampling simulation
#'
#' @description
#' End-to-end runner that can be called in two ways:
#'
#' **A. File-driven (legacy / CLI-friendly)**
#' Provide `init_file` (and optionally `interactions_file`). The function reads
#' configuration via [load_config()], simulates the community, places quadrats,
#' builds tables and figures, writes them to disk under `output_prefix`
#' (timestamped), and returns a results list.
#' **Interactions** can be supplied (i) via `interactions_file`, (ii) by putting
#' `INTERACTIONS_FILE = path/to/file` in the init file, or (iii) by embedding
#' inline rules with `INTERACTIONS_EDGELIST = c("A,B-D,0.8", "C,A,1.2", "E,*,0.95")`
#' directly in the init file.
#'
#' **B. Programmatic (as used in the README)**
#' Provide an in-memory parameter list `P` (typically from [load_config()])
#' and an `sf` polygon `domain` (or let the function create one). Interactions
#' are resolved from arguments or fields inside `P`: you may pass an external
#' `interactions_file`, set `P$INTERACTIONS_FILE` (pointer), or set
#' `P$INTERACTIONS_EDGELIST` (inline rules). By default `write_outputs = TRUE`,
#' so the same files are written; set `write_outputs = FALSE` to skip writing and
#' just get the results list back.
#'
#' Internally it:
#' \enumerate{
#'   \item resolves parameters (`P`) and optional interspecific interactions,
#'   \item simulates a heterogeneous community (point process + environment),
#'   \item places sampling quadrats using the requested scheme,
#'   \item compiles site \eqn{\times} species tables and per-site environment,
#'   \item (optionally) writes a 2\eqn{\times}2 map panel and an advanced multi-plot panel,
#'   \item (optionally) appends a human-readable text report, and
#'   \item returns a named results list.
#' }
#'
#' @param init_file Character path to a text configuration (KEY = value). If
#'   supplied, this takes precedence over `P`/`domain` inputs (file-driven mode).
#'   See [load_config()] for supported keys.
#' @param interactions_file Optional character path for a separate interactions
#'   config. If `NULL`, the function will look inside the init file (or `P`) for
#'   either `INTERACTIONS_FILE` (pointer) or `INTERACTIONS_EDGELIST` (inline rules).
#'   Falls back to neutral interactions (all 1s; radius 0) if nothing is provided.
#' @param output_prefix Base path for outputs (timestamp appended). If `NULL`,
#'   derives from `P$OUTPUT_PREFIX` or defaults to `"simulation_output"`.
#'   Ignored when `write_outputs = FALSE`.
#' @param domain Optional `sf` polygon (study area). If `NULL`, a default domain
#'   is created by [create_sampling_domain()]. Ignored in file-driven mode
#'   when `init_file` is provided (a domain is still created unless you supply
#'   your own).
#' @param P Optional fully materialized parameter list (often from
#'   [load_config()]). If `init_file` is **not** given, `P` is used as-is.
#' @param write_outputs Logical; write CSVs/figures/report to disk? Default `TRUE`.
#' @param interactions_validate Logical; validate the resolved interaction matrix
#'   with [validate_interactions()] (default `TRUE`).
#' @param interactions_strict Logical; if validation finds problems, stop with
#'   an error (`TRUE`) or just warn (`FALSE`). Default `TRUE`.
#' @param interactions_print Logical; pretty-print a compact summary of the
#'   interaction matrix with [print_interactions()] (default `TRUE`).
#' @param interactions_top_n Integer; cap the number of non-1.0 entries shown by
#'   the pretty-printer. Default `20`.
#'
#' @return
#' A named list with components:
#' \itemize{
#'   \item `P` -- resolved parameters
#'   \item `domain` -- `sf` polygon(s) of the study area
#'   \item `species_dist` -- `sf` points with `species` and attached env columns
#'   \item `quadrats` -- `sf` polygons of sampled quadrats
#'   \item `env_gradients` -- data.frame grid of environmental fields
#'   \item `abund_matrix` -- site \eqn{\times} species table (first column `site`)
#'   \item `site_coords` -- data.frame with `site, x, y` (quadrat centroids)
#' }
#' When `write_outputs = TRUE`, files are written to disk under `output_prefix`.
#'
#' @section Files written (when `write_outputs = TRUE`):
#' \itemize{
#'   \item `*_abundances.csv` -- site \eqn{\times} species matrix
#'   \item `*_environments.csv` -- mean environment per quadrat
#'   \item `*_quadrat_centroids.csv` -- centroids table
#'   \item `*_fig_panel.png` -- 2\eqn{\times}2 spatial panel
#'   \item `*_fig_advanced_panel.png` -- advanced diagnostics (if enabled)
#'   \item `*_report.txt` -- textual analysis report
#' }
#'
#' @seealso
#' [load_config()], [load_interactions()], [load_interactions_inline()],
#' [validate_interactions()], [print_interactions()],
#' [generate_heterogeneous_distribution()], [create_abundance_matrix()],
#' [calculate_quadrat_environment()], [plot_spatial_sampling()],
#' [generate_advanced_panel()], [generate_full_report()]
#'
#' @examples
#' \dontrun{
#' ## A) File-driven
#' run_spatial_simulation(
#'   init_file = "inst/examples/spesim_init_complete.txt",
#'   interactions_file = NULL,
#'   output_prefix = "out/run"
#' )
#'
#' ## B) Programmatic (as in README)
#' dom <- create_sampling_domain()
#' init <- system.file("examples", "spesim_init_complete.txt", package = "spesim")
#' P <- load_config(init)
#' res <- run_spatial_simulation(
#'   domain = dom, P = P,
#'   interactions_file = NULL,
#'   write_outputs = FALSE
#' )
#' str(res$abund_matrix)
#' }
#' @export
run_spatial_simulation <- function(init_file = NULL,
                                   interactions_file = NULL,
                                   output_prefix = NULL,
                                   domain = NULL,
                                   P = NULL,
                                   write_outputs = TRUE,
                                   interactions_validate = TRUE,
                                   interactions_strict = TRUE,
                                   interactions_print = TRUE,
                                   interactions_top_n = 20) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # ---- tiny helper: resolve a relative path against an anchor file ----------
  .is_abs_path <- function(path) {
    grepl("^(/|[A-Za-z]:[\\/])", path)
  }
  .resolve_rel <- function(path, anchor_file) {
    if (is.null(path) || !nzchar(path)) {
      return(NULL)
    }
    if (.is_abs_path(path)) {
      return(path)
    }
    file.path(dirname(anchor_file), path)
  }

  # --- 0) Resolve configuration source ---------------------------------------
  file_mode <- !is.null(init_file)
  if (file_mode) {
    P <- load_config(init_file)
  } else {
    if (is.null(P)) stop("Provide either `init_file` or an in-memory parameter list `P`.")
  }

  # Output prefix (only matters if writing)
  output_prefix <- output_prefix %||% P$OUTPUT_PREFIX %||% "simulation_output"
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  stamped_prefix <- paste0(output_prefix, "_", timestamp)
  output_dir <- dirname(stamped_prefix)

  # Domain
  domain <- domain %||% create_sampling_domain()

  # --- Interactions: precedence & new inline/pointer support -----------------
  # 1) explicit function argument wins
  src_file <- interactions_file

  # 2) if not given, check init/P for a pointer to a file
  if (is.null(src_file) && !is.null(P$INTERACTIONS_FILE)) {
    if (file_mode) {
      src_file <- .resolve_rel(as.character(P$INTERACTIONS_FILE), init_file)
    } else {
      src_file <- as.character(P$INTERACTIONS_FILE)
    }
  }

  # 3) load from one of: file | inline edgelist | fallback neutral
  if (!is.null(src_file)) {
    I <- load_interactions(src_file, P$N_SPECIES)
  } else if (!is.null(P$INTERACTIONS_EDGELIST)) {
    I <- load_interactions_inline(
      rules = as.character(P$INTERACTIONS_EDGELIST),
      n_species = P$N_SPECIES,
      radius = as.numeric(P$INTERACTION_RADIUS %||% 0)
    )
  } else {
    # neutral defaults
    I <- list(
      radius = as.numeric(P$INTERACTION_RADIUS %||% 0),
      matrix = matrix(1, P$N_SPECIES, P$N_SPECIES,
        dimnames = list(LETTERS[1:P$N_SPECIES], LETTERS[1:P$N_SPECIES])
      )
    )
  }

  # optional: validate/print
  if (isTRUE(interactions_validate)) {
    validate_interactions(
      I,
      spp_names = LETTERS[1:P$N_SPECIES],
      stop_on_error = isTRUE(interactions_strict)
    )
  }
  if (isTRUE(interactions_print)) {
    print_interactions(I, digits = 3, top_n = interactions_top_n)
  }

  # assign onto P after successful resolution
  P$INTERACTION_RADIUS <- I$radius
  P$INTERACTION_MATRIX <- I$matrix

  # --- 1) Optional report sink ------------------------------------------------
  report_path <- paste0(stamped_prefix, "_report.txt")
  con <- NULL
  if (isTRUE(write_outputs)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    con <- file(report_path, open = "wt")
    sink(con, type = "output")
    sink(con, type = "message")
    on.exit(
      {
        num_sinks <- sink.number()
        if (num_sinks > 0) for (i in seq_len(num_sinks)) sink()
        try(close(con), silent = TRUE)
      },
      add = TRUE
    )
  }

  # --- 2) Run core simulation -------------------------------------------------
  results_list <- NULL
  tryCatch(
    {
      cat("========== INITIALISING SPATIAL SAMPLING SIMULATION ==========\n")
      cat("========== SIMULATION PARAMETERS ==========\n")
      cat(paste(capture.output(dput(P[!names(P) %in% "QUADRAT_SIZE"])), collapse = "\n"), "\n\n")

      cat("========== RUNNING SIMULATION ==========\n")
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

      # --- 3) Save outputs (optional) ------------------------------------------
      if (isTRUE(write_outputs)) {
        cat("\n========== SAVING OUTPUT FILES ==========\n")
        f_abund <- paste0(stamped_prefix, "_abundances.csv")
        f_env <- paste0(stamped_prefix, "_environments.csv")
        f_coord <- paste0(stamped_prefix, "_quadrat_centroids.csv")
        utils::write.csv(abund_matrix, f_abund, row.names = FALSE)
        utils::write.csv(site_env, f_env, row.names = FALSE)
        utils::write.csv(site_coords, f_coord, row.names = FALSE)

        # 2x2 spatial panel
        p1 <- plot_spatial_sampling(domain, species_dist, quadrats, P)
        p2 <- plot_spatial_sampling(domain, species_dist, quadrats, P, TRUE, env_gradients, "temperature_C")
        p3 <- plot_spatial_sampling(domain, species_dist, quadrats, P, TRUE, env_gradients, "elevation_m")
        p4 <- plot_spatial_sampling(domain, species_dist, quadrats, P, TRUE, env_gradients, "rainfall_mm")
        panel_plot <- patchwork::wrap_plots(list(
          p1 + ggplot2::theme(legend.position = "right"),
          p2 + ggplot2::theme(legend.position = "right"),
          p3 + ggplot2::theme(legend.position = "right"),
          p4 + ggplot2::theme(legend.position = "right")
        ), ncol = 2) + patchwork::plot_layout(guides = "collect")
        f_panel <- paste0(stamped_prefix, "_fig_panel.png")
        ggplot2::ggsave(f_panel, panel_plot, width = 14, height = 10, dpi = 300, bg = "white")
        cat(sprintf("Main panel figure saved as %s\n", f_panel))
        cat("Data tables saved.\n")
      }

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
  if (!is.null(results_list) && isTRUE(results_list$P$ADVANCED_ANALYSIS) && isTRUE(write_outputs)) {
    tryCatch(
      {
        cat("\n========== ADVANCED ANALYSIS ==========\n")
        f_adv <- paste0(stamped_prefix, "_fig_advanced_panel.png")
        ap <- generate_advanced_panel(results_list)
        ggplot2::ggsave(f_adv, ap, width = 12, height = 14, dpi = 300, bg = "white")
        cat(sprintf("Advanced analysis panel saved as %s\n", f_adv))
      },
      error = function(e) {
        cat("\nERROR during advanced analysis:\n", e$message, "\n", sep = "")
      }
    )
  }

  # --- 5) Text report (optional) ---------------------------------------------
  if (!is.null(results_list) && isTRUE(write_outputs)) {
    cat("\n", generate_full_report(results_list), "\n", sep = "")
    cat("Outputs saved to: ", normalizePath(dirname(stamped_prefix)), "\n", sep = "")
  }

  # Return results
  results_list
}
