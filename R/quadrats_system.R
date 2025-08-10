#' Systematic Grid Quadrat Placement
#'
#' @description
#' Places axis‑aligned rectangular quadrats at the centers of a regular grid
#' spanning the domain’s bounding box, and retains only those quadrats whose
#' full geometry lies inside the domain polygon. The grid density is chosen to
#' approximate n_quadrats total placements by respecting the domain’s
#' aspect ratio (so cells are roughly square in map space while matching the
#' requested count).
#'
#' @details
#' * The grid is defined over the domain’s bounding box using nx × ny grid
#'   centers, where nx ≈ sqrt(n_quadrats / aspect_ratio) and
#'   ny ≈ aspect_ratio × nx, with aspect_ratio = height / width of the
#'   bounding box. This produces roughly nx * ny ≈ n_quadrats centers.
#' * For each grid center, an axis‑aligned rectangle of size
#'   quadrat_size = c(width, height) is created. Quadrats are kept only if
#'   they are fully within the domain (using sf::st_within()).
#' * If the domain is very irregular or narrow, the number of retained quadrats
#'   can be less than n_quadrats. The function returns all valid quadrats
#'   and emits a warning if none fit.
#' * The domain should use a projected (planar) CRS so that distances and
#'   sizes are in linear units (e.g., meters). If using lon/lat, reproject first
#'   (e.g., to UTM).
#'
#' @param domain An sf polygon/multipolygon representing the sampling
#'   region. Must carry a valid projected CRS suitable for linear measurements.
#' @param n_quadrats Integer (≥1). Approximate number of quadrats to place.
#'   The actual number returned is the count of grid‑centered quadrats that
#'   fully fit inside domain.
#' @param quadrat_size Numeric vector of length 2, c(width, height), giving
#'   the quadrat dimensions in the same linear units as domain
#'   (e.g., meters). Both values must be positive and finite.
#'
#' @return An sf object (polygons) with:
#' \describe{
#'   \item{quadrat_id}{Sequential integer identifier for each placed quadrat.}
#'   \item{geometry}{Axis‑aligned rectangular polygon for the quadrat.}
#' }
#'
#' @seealso
#' \code{\link{place_quadrats}} (random, non‑overlapping),
#' \code{\link{place_quadrats_tiled}} (systematic tiling by cell size),
#' \code{\link{place_quadrats_transect}} (parallel transects),
#' \code{\link{place_quadrats_voronoi}} (Voronoi‑based centers)
#'
#' @examples
#' \dontrun{
#' library(sf)
#' dom <- create_sampling_domain()
#' qs  <- place_quadrats_systematic(
#'   domain = dom,
#'   n_quadrats = 24,
#'   quadrat_size = c(1.5, 1.5)
#' )
#' plot(st_geometry(dom), col = “grey95”, border = “grey60”)
#' plot(st_geometry(qs), add = TRUE, border = “black”, lwd = 1)
#' }
place_quadrats_systematic <- function(domain, n_quadrats, quadrat_size) {
  bbox <- sf::st_bbox(domain)
  aspect_ratio <- (bbox["ymax"] - bbox["ymin"]) / (bbox["xmax"] - bbox["xmin"])
  nx <- round(sqrt(n_quadrats / aspect_ratio))
  ny <- round(aspect_ratio * nx)
  candidate_centers <- sf::st_make_grid(domain, n = c(nx, ny), what = "centers")
  candidate_quadrats <- sf::st_sfc(lapply(candidate_centers, function(pt) create_quadrat_from_center(pt, quadrat_size)), crs = sf::st_crs(domain))
  within_mat <- sf::st_within(candidate_quadrats, domain, sparse = FALSE)
  is_within <- within_mat[, 1, drop = TRUE]
  valid_quadrats <- candidate_quadrats[is_within]
  if (length(valid_quadrats) == 0) {
    warning("Systematic sampling failed to place any quadrats.")
    return(sf::st_sf(quadrat_id = integer(0), geometry = sf::st_sfc(crs = sf::st_crs(domain))))
  }
  sf::st_sf(quadrat_id = seq_len(length(valid_quadrats)), geometry = valid_quadrats)
}
