#' Run a complete spesim simulation (recommended)
#'
#' @description
#' High-level, teaching-friendly wrapper around [run_spatial_simulation()].
#' It supports both file-driven and programmatic workflows, adds optional
#' progress messaging, and returns an S3 result object with stable components.
#'
#' Use this as the main entry point in workshops and teaching material.
#'
#' @param config Either a path to an init file (character scalar), or an in-memory
#'   parameter list `P` (typically from [load_config()]).
#' @param domain Optional `sf` polygon study area (see [create_sampling_domain()]).
#'   If `NULL`, a default domain is created.
#' @param interactions_file Optional path to an interactions config file.
#' @param output_prefix Base output prefix (timestamp is appended) when
#'   `write_outputs = TRUE`.
#' @param write_outputs Logical; write CSVs/figures/report to disk? Default `FALSE`
#'   for a smoother interactive experience.
#' @param seed Optional integer. If supplied, overrides `P$SEED` (and is applied
#'   via `set.seed()` before simulation).
#' @param quiet Logical; suppress most messages.
#' @param ... Passed through to [run_spatial_simulation()].
#'
#' @return An object of class `spesim_result` (a named list) with components:
#' 
#' - `P`: resolved parameters
#' - `domain`: study area (`sf`)
#' - `species_dist`: simulated individuals (`sf` points)
#' - `quadrats`: sampling quadrats (`sf` polygons; class `spesim_quadrats`)
#' - `env_gradients`: environmental grid
#' - `abund_matrix`: site × species abundances
#' - `site_coords`: quadrat centroids
#'
#' @seealso [spesim_demo()], [load_config()], [run_spatial_simulation()]
#' @export
spesim_run <- function(config,
                       domain = NULL,
                       interactions_file = NULL,
                       output_prefix = NULL,
                       write_outputs = FALSE,
                       seed = NULL,
                       quiet = FALSE,
                       ...) {
  if (is.character(config) && length(config) == 1L) {
    init_file <- config
    P <- NULL
  } else if (is.list(config)) {
    init_file <- NULL
    P <- config
  } else {
    stop("`config` must be either a path to an init file, or an in-memory parameter list P.")
  }

  if (!quiet) {
    message("spesim: running simulation")
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("`seed` must be a single finite integer-like value.")
    }
    seed <- as.integer(seed)
    # If P exists, override; otherwise load_config will set from file.
    if (!is.null(P)) P$SEED <- seed
    set.seed(seed)
  }

  res <- run_spatial_simulation(
    init_file = init_file,
    interactions_file = interactions_file,
    output_prefix = output_prefix,
    domain = domain,
    P = P,
    write_outputs = write_outputs,
    ...
  )

  if (is.null(res)) {
    stop("Simulation failed; see messages above.")
  }

  # tag quadrats for downstream methods
  if (!is.null(res$quadrats)) {
    class(res$quadrats) <- unique(c("spesim_quadrats", class(res$quadrats)))
  }

  class(res) <- unique(c("spesim_result", class(res)))
  res
}


#' @export
print.spesim_result <- function(x, ...) {
  cat("<spesim_result>\n")
  if (!is.null(x$P)) {
    P <- x$P
    cat(sprintf("  species: %s\n", P$N_SPECIES %||% "?"))
    cat(sprintf("  individuals: %s\n", P$N_INDIVIDUALS %||% "?"))
    cat(sprintf("  sampling scheme: %s\n", P$SAMPLING_SCHEME %||% "?"))
    cat(sprintf("  n_quadrats: %s\n", P$N_QUADRATS %||% "?"))
    if (!is.null(P$SEED)) cat(sprintf("  seed: %s\n", P$SEED))
  }
  if (!is.null(x$species_dist)) cat(sprintf("  points: %s\n", nrow(x$species_dist)))
  if (!is.null(x$quadrats)) cat(sprintf("  quadrats: %s\n", nrow(x$quadrats)))
  invisible(x)
}


#' @export
summary.spesim_result <- function(object, ...) {
  x <- object
  out <- list(
    n_species = x$P$N_SPECIES %||% NA_integer_,
    n_individuals = x$P$N_INDIVIDUALS %||% NA_integer_,
    sampling_scheme = x$P$SAMPLING_SCHEME %||% NA_character_,
    n_quadrats = nrow(x$quadrats %||% data.frame()),
    seed = x$P$SEED %||% NA_integer_
  )
  class(out) <- "summary_spesim_result"
  out
}

#' @export
print.summary_spesim_result <- function(x, ...) {
  cat("spesim result summary\n")
  cat(sprintf("  species         : %s\n", x$n_species))
  cat(sprintf("  individuals     : %s\n", x$n_individuals))
  cat(sprintf("  sampling scheme : %s\n", x$sampling_scheme))
  cat(sprintf("  quadrats        : %s\n", x$n_quadrats))
  cat(sprintf("  seed            : %s\n", x$seed))
  invisible(x)
}
