# Sampling frame / edge effects helpers

#' Estimate the effective sampling frame for quadrat centres
#'
#' @description
#' For an axis-aligned rectangle (quadrat) to lie fully within an irregular
#' polygon domain, the quadrat centre must lie inside a shrunk version of the
#' domain.
#'
#' Exactly, this shrink is an erosion in an L-infinity metric (half-width and
#' half-height constraints). For simplicity and robustness in teaching/method
#' testing, this helper uses a conservative Euclidean buffer based on the
#' half-diagonal of the quadrat.
#'
#' This provides an interpretable, single-number summary of boundary constraints:
#'
#' - `safe_area_fraction` = area(safe_domain) / area(domain)
#'
#' Values close to 1 mean boundary exclusion is mild; values near 0 mean that
#' only a small interior region can host fully-contained quadrats.
#'
#' @param domain An `sf` polygon/multipolygon.
#' @param quadrat_size Numeric length-2, c(width, height), in the same units as
#'   the domain CRS.
#'
#' @return A list with `safe_domain` (sf), `safe_area_fraction` (numeric), and
#'   `buffer_dist` (numeric half-diagonal).
#' @keywords internal
estimate_sampling_frame <- function(domain, quadrat_size) {
  stopifnot(length(quadrat_size) == 2)
  w <- as.numeric(quadrat_size[1]); h <- as.numeric(quadrat_size[2])
  if (!is.finite(w) || !is.finite(h) || w <= 0 || h <= 0) {
    stop("`quadrat_size` must be positive and finite")
  }

  buffer_dist <- sqrt(w^2 + h^2) / 2

  dom_u <- sf::st_union(domain)
  area_dom <- as.numeric(sf::st_area(dom_u))

  safe_domain <- suppressWarnings(sf::st_buffer(dom_u, -buffer_dist))

  if (inherits(safe_domain, "sfc") || inherits(safe_domain, "sfg")) {
    safe_sf <- sf::st_as_sf(safe_domain)
  } else {
    safe_sf <- sf::st_as_sf(sf::st_geometry(safe_domain))
  }

  # Handle empty/invalid safe domain
  safe_area <- 0
  if (!sf::st_is_empty(safe_domain) && length(safe_domain) > 0) {
    safe_area <- as.numeric(sf::st_area(sf::st_union(safe_domain)))
  }

  frac <- if (is.finite(area_dom) && area_dom > 0) safe_area / area_dom else NA_real_

  list(
    safe_domain = safe_sf,
    safe_area_fraction = frac,
    buffer_dist = buffer_dist
  )
}
