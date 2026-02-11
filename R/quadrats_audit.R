# Quadrat placement audit helpers

#' Validate a placement_audit attribute
#'
#' @description
#' Quadrat placement functions in spesim attach a `placement_audit` attribute to
#' the returned `sf` object. This is used by [audit_sampling_scheme()] and the
#' report's conceptual audit section.
#'
#' This validator is intentionally permissive: it checks for a minimal common
#' schema and normalises field types where possible, but it does not error on
#' scheme-specific fields.
#'
#' @param x An object with (optional) `attr(x, "placement_audit")`.
#' @param stop_on_error Logical; stop on validation errors? Default `FALSE`.
#'
#' @return A normalised `placement_audit` list (possibly empty).
#' @keywords internal
.validate_placement_audit <- function(x, stop_on_error = FALSE) {
  aud <- attr(x, "placement_audit")
  if (is.null(aud)) return(list())
  if (!is.list(aud)) {
    msg <- "placement_audit must be a list"
    if (isTRUE(stop_on_error)) stop(msg) else warning(msg)
    return(list())
  }

  # minimal fields
  aud$scheme <- if (!is.null(aud$scheme)) as.character(aud$scheme) else NA_character_
  aud$n_requested <- if (!is.null(aud$n_requested)) as.integer(aud$n_requested) else NA_integer_
  aud$n_returned <- if (!is.null(aud$n_returned)) as.integer(aud$n_returned) else NA_integer_

  # optional standard fields
  .int_fields <- c(
    "attempts", "max_attempts",
    "reject_boundary", "reject_overlap",
    "candidates_total", "candidates_valid",
    "transects_total", "transects_with_segment",
    "seeds", "cells_total", "cells_suitable"
  )
  for (nm in .int_fields) {
    if (!is.null(aud[[nm]])) aud[[nm]] <- as.integer(aud[[nm]])
  }

  .dbl_fields <- c("safe_area_fraction", "boundary_exclusion")
  for (nm in .dbl_fields) {
    if (!is.null(aud[[nm]])) aud[[nm]] <- as.numeric(aud[[nm]])
  }

  if (!is.null(aud$safe_area_empty)) aud$safe_area_empty <- as.logical(aud$safe_area_empty)

  # sanity
  if (!is.na(aud$n_requested) && aud$n_requested < 0) {
    msg <- "placement_audit$n_requested must be >= 0"
    if (isTRUE(stop_on_error)) stop(msg) else warning(msg)
  }
  if (!is.na(aud$n_returned) && aud$n_returned < 0) {
    msg <- "placement_audit$n_returned must be >= 0"
    if (isTRUE(stop_on_error)) stop(msg) else warning(msg)
  }

  aud
}

#' Attach a placement_audit attribute to a quadrats sf object
#'
#' @param quadrats sf object.
#' @param audit list.
#'
#' @return quadrats with normalised `placement_audit` attribute.
#' @keywords internal
.attach_placement_audit <- function(quadrats, audit) {
  if (!is.list(audit)) stop("`audit` must be a list")
  attr(quadrats, "placement_audit") <- audit
  # normalise
  attr(quadrats, "placement_audit") <- .validate_placement_audit(quadrats, stop_on_error = FALSE)
  quadrats
}
