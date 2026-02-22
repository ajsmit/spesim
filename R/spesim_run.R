# Internal helper: produce a different seed on each call.
# Used to generate a varying default domain for spesim_run() when no seed is
# supplied.
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

# Internal: load and validate config from a file path or an in-memory list.
# Keeps the type-dispatch and error message in one place so spesim_run() body
# does not have to branch on the config type.
.resolve_run_config <- function(config) {
  if (is.character(config) && length(config) == 1L) {
    load_config(config)   # handles parsing, defaults, validation, path resolution
  } else if (is.list(config)) {
    config
  } else {
    stop("`config` must be either a path to an init file, or an in-memory parameter list P.")
  }
}

#' Run a complete spesim simulation (recommended)
#'
#' @description
#' `spesim_run()` is the **main public entry point** for the package.
#'
#' It supports both:
#' - **File-driven workflows**: pass a path to an init file (KEY = value)
#' - **Programmatic workflows**: pass an in-memory parameter list `P`
#'
#' Internally, `spesim_run()` orchestrates a full run by:
#' 1) resolving configuration (and optional seed override),
#' 2) creating (or accepting) a sampling domain,
#' 3) resolving interspecific interactions (file pointer / inline rules),
#' 4) calling the pure simulation engine, and
#' 5) optionally writing an output bundle (CSVs, figures, report).
#'
#' Compared to [run_spatial_simulation()], this function returns a stable S3
#' object (`spesim_result`) and keeps filesystem I/O optional.
#'
#' @param config Either a path to an init file (character scalar), or an
#'   in-memory parameter list `P` (typically from [load_config()]).
#' @param domain Optional `sf` polygon study area (see [create_sampling_domain()]).
#'   If `NULL`, a default domain is created.
#' @param interactions_file Optional path to an interactions config file.
#'   If `NULL`, the function will look for `INTERACTIONS_FILE` or
#'   `INTERACTIONS_EDGELIST` in `config`/`P`.
#' @param output_prefix Base output prefix (timestamp is appended) when
#'   `write_outputs = TRUE`.
#' @param write_outputs Logical; write CSVs/figures/report to disk? Default
#'   `FALSE` for a smoother interactive experience.
#' @param seed Optional integer. If supplied, it (i) overrides `P$SEED` for the
#'   simulation, and (ii) is also used to generate a reproducible default domain
#'   when `domain = NULL`.
#' @param quiet Logical; suppress most messages.
#' @param interactions_validate Logical; validate resolved interactions.
#' @param interactions_strict Logical; if validation finds problems, stop (`TRUE`)
#'   or warn (`FALSE`).
#' @param interactions_print Logical; pretty-print a compact interactions summary.
#' @param interactions_top_n Integer; cap the number of non-1.0 entries shown.
#' @param report_include_audit Logical; include the conceptual audit section in
#'   the written report (when `write_outputs = TRUE`).
#' @param report_audit_top_n Integer; rows to show in the audit section.
#' @param ... Reserved for future extensions. Currently unused.
#'
#' @return An object of class `spesim_result` (a named list) with components:
#' 
#' - `P`: resolved parameters
#' - `domain`: study area (`sf`)
#' - `species_dist`: simulated individuals (`sf` points)
#' - `quadrats`: sampling quadrats (`sf` polygons; class `spesim_quadrats`)
#' - `env_gradients`: environmental grid
#' - `abund_matrix`: site × species abundances
#' - `site_env`: per-quadrat mean environment
#' - `site_coords`: quadrat centroids
#' - `files` (optional): written file paths (when `write_outputs = TRUE`)
#'
#' @seealso [load_config()], [create_sampling_domain()],
#'   [generate_full_report()], [read_latest_report()]
#' @export
spesim_run <- function(config,
                       domain = NULL,
                       interactions_file = NULL,
                       output_prefix = NULL,
                       write_outputs = FALSE,
                       seed = NULL,
                       quiet = FALSE,
                       interactions_validate = TRUE,
                       interactions_strict = TRUE,
                       interactions_print = TRUE,
                       interactions_top_n = 20,
                       report_include_audit = TRUE,
                       report_audit_top_n = 6,
                       ...) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  dots <- list(...)
  if (length(dots)) {
    stop("Unused argument(s): ", paste(names(dots), collapse = ", "))
  }

  # Validate seed
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("`seed` must be a single finite integer-like value.")
    }
    seed <- as.integer(seed)
  }

  # Resolve config: file path → load_config() (parsing, validation, path
  # resolution all handled there); list → used as-is.
  P <- .resolve_run_config(config)

  # Seed override: controls both simulation reproducibility and default domain.
  if (!is.null(seed)) {
    P$SEED <- seed
  }

  # Materialise derived fields for programmatic workflows (idempotent)
  P <- spesim_materialize_P(P)

  # Domain creation policy
  if (is.null(domain)) {
    if (!is.null(seed)) {
      set.seed(seed)
    } else {
      set.seed(.auto_domain_seed())
    }
    domain <- create_sampling_domain()
  }

  if (!quiet) message("spesim: running simulation")

  # Interactions
  ir <- spesim_resolve_interactions(
    P,
    interactions_file = interactions_file,
    interactions_validate = interactions_validate,
    interactions_strict = interactions_strict,
    interactions_print = interactions_print,
    interactions_top_n = interactions_top_n
  )
  P <- ir$P

  # Engine
  res <- spesim_engine(P, domain)
  res <- new_spesim_result(res)

  # Optional output bundle
  if (isTRUE(write_outputs)) {
    paths <- spesim_write_outputs(
      res,
      output_prefix = output_prefix,
      include_audit = report_include_audit,
      audit_top_n = report_audit_top_n,
      quiet = quiet
    )
    res$files <- paths
  }

  res
}
