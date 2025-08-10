#' Setup a fresh spesim project: main init, interactions init, and SxS matrix CSV
#'
#' This creates:
#'   1) {out_dir}/{spesim_init}           -- main config (NO interaction keys)
#'   2) {out_dir}/{interactions_init}     -- interaction config pointing to...
#'   3) {out_dir}/{interactions_matrix}   -- SxS matrix of 1s with A.. species headers
#'
#' You’ll pass N_SPECIES here once; all files will match it.
#'
#' @param n_species Integer number of species (>=1)
#' @param out_dir Output directory (created if missing)
#' @param spesim_init File name for main init
#' @param interactions_init File name for interactions init
#' @param interactions_matrix File name for the SxS matrix CSV
#' @param interaction_radius Default radius; 0 disables interactions
#' @param overwrite Logical; if FALSE, won't clobber existing files
#' @param seed RNG seed written into main init for reproducibility
#' @return Invisibly, a list with file paths
setup_spesim_project <- function(
    n_species,
    out_dir = "spesim_project",
    spesim_init = "spesim_init.txt",
    interactions_init = "interactions_init.txt",
    interactions_matrix = "interactions_matrix.csv",
    interaction_radius = 2.0,
    overwrite = FALSE,
    seed = 77) {
  # ---- validation ----
  stopifnot(is.numeric(n_species), length(n_species) == 1, is.finite(n_species), n_species >= 1)
  n_species <- as.integer(n_species)

  # Make dir
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Resolve paths
  p_spec <- file.path(out_dir, spesim_init)
  p_inter <- file.path(out_dir, interactions_init)
  p_matrix <- file.path(out_dir, interactions_matrix)

  # Overwrite guard
  check_write <- function(path) {
    if (file.exists(path) && !overwrite) {
      stop("File exists and overwrite=FALSE: ", path)
    }
  }
  check_write(p_spec)
  check_write(p_inter)
  check_write(p_matrix)

  # ---- build the SxS matrix of 1s with A.. species labels ----
  spp <- utils::head(LETTERS, n_species)
  M <- matrix(1.0,
    nrow = n_species, ncol = n_species,
    dimnames = list(spp, spp)
  )
  utils::write.csv(M, p_matrix, row.names = TRUE, quote = FALSE)

  # ---- write main init (NO interaction keys here) ----
  # Keep these sane/editable. You can trim or expand as you like.
  main_init_lines <- c(
    sprintf("SEED = %d", seed),
    sprintf("N_SPECIES = %d", n_species),
    "N_INDIVIDUALS = 2000",
    "DOMINANT_FRACTION = 0.30",
    "FISHER_ALPHA = 3.0",
    "FISHER_X = 0.95",
    "",
    "# Gradient assignments (optional demo; can be edited later)",
    sprintf("GRADIENT_SPECIES = c(%s)", paste(sprintf('"%s"', spp), collapse = ", ")),
    sprintf(
      "GRADIENT_ASSIGNMENTS = c(%s)",
      paste(rep('"temperature"', n_species), collapse = ", ")
    ),
    "GRADIENT_OPTIMA = 0.5",
    "GRADIENT_TOLERANCE = 0.12",
    "SAMPLING_RESOLUTION = 50",
    "ENVIRONMENTAL_NOISE = 0.05",
    "",
    "# Quadrat sampling",
    'SAMPLING_SCHEME = "random"', # random|systematic|transect|tiled|voronoi
    "N_QUADRATS = 20",
    'QUADRAT_SIZE_OPTION = "medium"', # small|medium|large
    "N_TRANSECTS = 1",
    "N_QUADRATS_PER_TRANSECT = 8",
    "TRANSECT_ANGLE = 90",
    "VORONOI_SEED_FACTOR = 10",
    "",
    "# Plot styling",
    "POINT_SIZE = 0.2",
    "POINT_ALPHA = 1.0",
    "QUADRAT_ALPHA = 0.05",
    'BACKGROUND_COLOUR = "white"',
    'FOREGROUND_COLOUR = "#22223b"',
    'QUADRAT_COLOUR = "black"',
    "",
    "ADVANCED_ANALYSIS = FALSE"
    # NOTE: no INTERACTION_* keys here by design
  )
  writeLines(main_init_lines, con = p_spec, useBytes = TRUE)

  # ---- write interactions init that points at the matrix we just made ----
  inter_lines <- c(
    sprintf("INTERACTION_RADIUS = %g", interaction_radius),
    sprintf("MATRIX_CSV = %s", basename(p_matrix)),
    "AUTO = FALSE"
    # (You can add EDGELIST_CSV later if you want that mode.)
  )
  writeLines(inter_lines, con = p_inter, useBytes = TRUE)

  message("Wrote:\n  - ", p_spec, "\n  - ", p_inter, "\n  - ", p_matrix)

  invisible(list(
    spesim_init = normalizePath(p_spec),
    interactions_init = normalizePath(p_inter),
    interactions_matrix = normalizePath(p_matrix)
  ))
}
