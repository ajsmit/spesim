#' Run the full spatial sampling simulation (legacy wrapper)
#'
#' @description
#' `run_spatial_simulation()` is the historical end-to-end runner for **spesim**.
#'
#' As of **spesim 0.4.0**, the recommended entry point is [spesim_run()], which
#' consolidates file-driven and programmatic workflows and returns a stable S3
#' object (`spesim_result`).
#'
#' This function is retained for backwards compatibility and delegates to
#' [spesim_run()]. It will emit a deprecation warning.
#'
#' @note **Deprecated** since spesim 0.4.0. Use [spesim_run()] instead.
#'   Calling this function triggers `.Deprecated("spesim_run")`, which emits
#'   a warning and prints a message pointing to the replacement.
#'
#' @param init_file Character path to a text configuration (KEY = value). If
#'   supplied, the config is read via [load_config()].
#' @param interactions_file Optional character path for a separate interactions
#'   config. If `NULL`, interactions are resolved from `P`.
#' @param output_prefix Base path for outputs (timestamp appended). If `NULL`,
#'   derives from `P$OUTPUT_PREFIX`.
#' @param domain Optional `sf` polygon (study area). If `NULL`, a domain is
#'   created.
#' @param P Optional parameter list. Used when `init_file` is not supplied.
#' @param write_outputs Logical; write CSVs/figures/report to disk? Default `TRUE`.
#' @param interactions_validate,interactions_strict,interactions_print,interactions_top_n
#'   Passed to [spesim_run()].
#'
#' @return A `spesim_result`.
#'
#' @seealso [spesim_run()]
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
  .Deprecated("spesim_run", package = "spesim")

  if (!is.null(init_file)) {
    # File-driven mode (legacy): load once, then pass the list into spesim_run.
    P <- load_config(init_file)

    # Resolve relative interactions-file pointer like the old runner did.
    if (!is.null(P$INTERACTIONS_FILE) && nzchar(P$INTERACTIONS_FILE)) {
      is_abs <- grepl("^(/|[A-Za-z]:[\\/])", P$INTERACTIONS_FILE)
      if (!is_abs) {
        P$INTERACTIONS_FILE <- file.path(dirname(init_file), P$INTERACTIONS_FILE)
      }
    }
  } else {
    if (is.null(P)) stop("Provide either `init_file` or an in-memory parameter list `P`.")
  }

  # Preserve legacy behaviour: if domain is missing, seed domain generation from
  # P$SEED (so init-file workflows remain reproducible by default).
  if (is.null(domain)) {
    if (!is.null(P$SEED) && is.finite(P$SEED)) set.seed(as.integer(P$SEED))
    domain <- create_sampling_domain()
  }

  spesim_run(
    config = P,
    domain = domain,
    interactions_file = interactions_file,
    output_prefix = output_prefix,
    write_outputs = write_outputs,
    seed = NULL,
    quiet = FALSE,
    interactions_validate = interactions_validate,
    interactions_strict = interactions_strict,
    interactions_print = interactions_print,
    interactions_top_n = interactions_top_n
  )
}
