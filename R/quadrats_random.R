#' Random Non‑Overlapping Quadrat Placement
#'
#' @description
#' Places axis‑aligned rectangular quadrats at random locations inside
#' domain, rejecting candidates that either (i) are not fully contained or
#' (ii) overlap any previously accepted quadrat. A simple rejection‑sampling
#' loop is used with a finite attempt budget.
#'
#' @details
#' * Each candidate center is drawn uniformly from the domain’s bounding box,
#'   shrunk so the full rectangle of size quadrat_size = c(width, height)
#'   fits inside. The resulting polygon is accepted only if it is entirely
#'   within domain and non‑overlapping with all accepted quadrats.
#' * The loop stops when n_quadrats are placed or after
#'   n_quadrats * 100 attempts (a warning is issued if the target is not met).
#' * CRS: Use a projected CRS (e.g., meters). quadrat_size is
#'   interpreted in the CRS linear units.
#' * Reproducibility: Set a seed (e.g., set.seed(...)) before calling to
#'   obtain repeatable placements.
#'
#' @param domain An sf polygon/multipolygon defining the sampling region.
#'   Must carry a valid projected CRS suitable for linear units.
#' @param n_quadrats Integer (≥1). Target number of non‑overlapping quadrats to place.
#' @param quadrat_size Numeric vector of length 2, c(width, height), giving
#'   quadrat dimensions in the same units as the domain CRS. Both values
#'   must be positive and finite.
#'
#' @return An sf object (polygons) with:
#' \describe{
#'   \item{quadrat_id}{Sequential integer identifier for each placed quadrat.}
#'   \item{geometry}{Axis‑aligned rectangular polygon for the quadrat.}
#' }
#' If no valid placement is possible, returns an empty sf with the same CRS
#' and emits a warning.
#'
#' @seealso
#' \code{\link{place_quadrats_systematic}} (regular grid),
#' \code{\link{place_quadrats_tiled}} (systematic tiling by cell size),
#' \code{\link{place_quadrats_voronoi}} (Voronoi‑based centers),
#' \code{\link{place_quadrats_transect}} (parallel transects)
#'
#' @examples
#' \dontrun{
#' library(sf)
#' dom <- create_sampling_domain()
#' set.seed(42)
#' qs <- place_quadrats(dom, n_quadrats = 20, quadrat_size = c(1.5, 1.5))
#' plot(st_geometry(dom), col = “grey95”, border = “grey60”)
#' plot(st_geometry(qs), add = TRUE, border = “black”, lwd = 1)
#' }
place_quadrats <- function(domain, n_quadrats, quadrat_size) {
  bbox <- sf::st_bbox(domain)
  quadrats <- list()
  attempts <- 0
  max_attempts <- n_quadrats * 100

  while (length(quadrats) < n_quadrats && attempts < max_attempts) {
    attempts <- attempts + 1

    x_center <- runif(1, bbox["xmin"] + quadrat_size[1] / 2, bbox["xmax"] - quadrat_size[1] / 2)
    y_center <- runif(1, bbox["ymin"] + quadrat_size[2] / 2, bbox["ymax"] - quadrat_size[2] / 2)

    center_pt_sfc <- sf::st_sfc(sf::st_point(c(x_center, y_center)), crs = sf::st_crs(domain))
    new_quadrat_poly <- create_quadrat_from_center(center_pt_sfc, quadrat_size) |>
      sf::st_sfc(crs = sf::st_crs(domain))

    within_mat <- sf::st_within(new_quadrat_poly, domain, sparse = FALSE)
    is_within <- isTRUE(within_mat[, 1, drop = TRUE][1])

    if (is_within) {
      is_overlapping <- FALSE
      if (length(quadrats) > 0) {
        existing_quadrats_sfc <- do.call(c, quadrats)
        overlap_mat <- sf::st_intersects(new_quadrat_poly, existing_quadrats_sfc, sparse = FALSE)
        is_overlapping <- any(overlap_mat[1, , drop = TRUE])
      }
      if (!is_overlapping) {
        quadrats[[length(quadrats) + 1]] <- new_quadrat_poly
      }
    }
  }

  if (length(quadrats) == 0) {
    return(sf::st_sf(quadrat_id = integer(0), geometry = sf::st_sfc(crs = sf::st_crs(domain))))
  }

  sfc <- do.call(c, quadrats)
  sf::st_sf(quadrat_id = seq_along(sfc), geometry = sfc)
}


#' Create an Axis‑Aligned Quadrat Polygon from a Center
#'
#' @description
#' Internal helper that builds an axis‑aligned rectangular polygon centered
#' on center_point with dimensions size = c(width, height).
#'
#' @param center_point An sf point (single feature) providing the quadrat center.
#'   Must have the same CRS intended for the output polygon.
#' @param size Numeric vector of length 2, c(width, height), giving the
#'   rectangle’s dimensions in the CRS units. Both values must be positive and finite.
#'
#' @return An sf polygon geometry (sfg) representing the rectangle.
#'
#' @keywords internal
#' Create a Quadrat Polygon From Center
#' @param center_point An sf POINT (sfg or sfc of length 1).
#' @param size numeric length-2, c(width, height).
#' @return sfg POLYGON (single rectangle).
#' @keywords internal
create_quadrat_from_center <- function(center_point, size) {
  # Coerce to sfc length 1, preserve CRS if present
  crs0 <- try(sf::st_crs(center_point), silent = TRUE)
  cp_sfc <- suppressWarnings(sf::st_sfc(center_point, crs = if (inherits(crs0, "crs")) crs0 else NA))
  xy <- sf::st_coordinates(cp_sfc)

  # Ensure 1x2 matrix with X,Y
  if (is.null(dim(xy))) {
    xy <- matrix(as.numeric(xy), nrow = 1, ncol = 2,
                 dimnames = list(NULL, c("X", "Y")))
  } else {
    colnames(xy) <- c("X", "Y")
  }

  cx <- as.numeric(xy[1, "X"]); cy <- as.numeric(xy[1, "Y"])
  half_w <- as.numeric(size[1]) / 2
  half_h <- as.numeric(size[2]) / 2

  ring <- cbind(
    c(cx - half_w, cx + half_w, cx + half_w, cx - half_w, cx - half_w),
    c(cy - half_h, cy - half_h, cy + half_h, cy + half_h, cy - half_h)
  )
  sf::st_polygon(list(ring))
}
