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
#'   If `NULL`, a default domain is created. By default (no `seed` supplied),
#'   the default domain is *randomly generated per run*; if `seed` is supplied,
#'   the default domain is reproducible.
#' @param interactions_file Optional path to an interactions config file.
#' @param output_prefix Base output prefix (timestamp is appended) when
#'   `write_outputs = TRUE`.
#' @param write_outputs Logical; write CSVs/figures/report to disk? Default `FALSE`
#'   for a smoother interactive experience.
#' @param seed Optional integer. If supplied, it (i) overrides `P$SEED` for the
#'   simulation, and (ii) is also used to generate a reproducible default
#'   domain when `domain = NULL`.
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
.auto_domain_seed <- local({
  nonce <- 0L
  function() {
    nonce <<- nonce + 1L
    # Use time + a monotone counter so that consecutive calls within the same
    # second still generate different domains.
    t <- as.numeric(Sys.time())
    base <- as.integer((t * 1e6) %% .Machine$integer.max)
    as.integer((base + nonce) %% .Machine$integer.max)
  }
})

spesim_run <- function(config,
                       domain = NULL,
                       interactions_file = NULL,
                       output_prefix = NULL,
                       write_outputs = FALSE,
                       seed = NULL,
                       quiet = FALSE,
                       ...) {
  # Resolve seed (spesim_run seed controls BOTH the default domain and the
  # simulation RNG for reproducibility).
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("`seed` must be a single finite integer-like value.")
    }
    seed <- as.integer(seed)
  }

  # If no domain is supplied, create one.
  #
  # Behaviour:
  # - seed provided  -> domain is reproducible across runs
  # - seed NULL      -> domain varies across calls (even though the simulator
  #                     itself may reseed internally via P$SEED)
  if (is.null(domain)) {
    if (!is.null(seed)) {
      set.seed(seed)
    } else {
      set.seed(.auto_domain_seed())
    }
    domain <- create_sampling_domain()
  }

  # Resolve config source. If a file path is given, we intentionally load it
  # here (programmatic mode) so we can override P$SEED when seed is provided,
  # and so we can always pass an explicit domain into run_spatial_simulation().
  init_file <- NULL
  if (is.character(config) && length(config) == 1L) {
    init_file <- config
    P <- load_config(init_file)

    # If init refers to a relative interactions file, resolve it relative to the
    # init file (mimics run_spatial_simulation(file_mode) behaviour).
    if (!is.null(P$INTERACTIONS_FILE) && nzchar(P$INTERACTIONS_FILE)) {
      is_abs <- grepl("^(/|[A-Za-z]:[\\/])", P$INTERACTIONS_FILE)
      if (!is_abs) {
        P$INTERACTIONS_FILE <- file.path(dirname(init_file), P$INTERACTIONS_FILE)
      }
    }
  } else if (is.list(config)) {
    P <- config
  } else {
    stop("`config` must be either a path to an init file, or an in-memory parameter list P.")
  }

  if (!quiet) {
    message("spesim: running simulation")
  }

  # Override simulation seed if requested.
  if (!is.null(seed)) {
    P$SEED <- seed
  }

  res <- run_spatial_simulation(
    init_file = NULL,
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
