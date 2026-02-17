#' Construct a spesim_result
#'
#' @description
#' Internal constructor used by [spesim_run()] and legacy wrappers.
#'
#' @param x A named list with the standard components.
#'
#' @return The same list with class `spesim_result`.
#'
#' @keywords internal
#' @noRd
new_spesim_result <- function(x) {
  if (is.null(x) || !is.list(x)) stop("`x` must be a list")

  # Minimal required components
  req <- c("P", "domain", "species_dist", "quadrats", "env_gradients", "abund_matrix", "site_env", "site_coords")
  miss <- setdiff(req, names(x))
  if (length(miss)) {
    stop("spesim_result is missing required component(s): ", paste(miss, collapse = ", "))
  }

  # Tag quadrats for downstream methods
  if (!is.null(x$quadrats)) {
    class(x$quadrats) <- unique(c("spesim_quadrats", class(x$quadrats)))
  }

  class(x) <- unique(c("spesim_result", class(x)))
  x
}

#' @export
print.spesim_result <- function(x, ...) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  cat("<spesim_result>\n")
  if (!is.null(x$P)) {
    P <- x$P
    cat(sprintf("  species         : %s\n", P$N_SPECIES %||% "?"))
    cat(sprintf("  individuals     : %s\n", P$N_INDIVIDUALS %||% "?"))
    cat(sprintf("  sampling scheme : %s\n", P$SAMPLING_SCHEME %||% "?"))
    cat(sprintf("  n_quadrats      : %s\n", P$N_QUADRATS %||% "?"))
    if (!is.null(P$SEED)) cat(sprintf("  seed            : %s\n", P$SEED))
  }
  if (!is.null(x$species_dist)) cat(sprintf("  points          : %s\n", nrow(x$species_dist)))
  if (!is.null(x$quadrats)) cat(sprintf("  quadrats        : %s\n", nrow(x$quadrats)))
  invisible(x)
}

#' @export
summary.spesim_result <- function(object, ...) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

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
