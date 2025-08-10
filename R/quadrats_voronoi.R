#' Voronoi‑based Quadrat Placement
#'
#' @description
#' Places rectangular quadrats by first generating a dense set of random seed
#' points inside the domain, computing the Voronoi tessellation, and then
#' identifying those Voronoi cells whose largest inscribed circle is big
#' enough to contain the requested quadrat (by diagonal). For each suitable
#' cell, the quadrat is placed at the cell’s “pole of inaccessibility”
#' (center of the inscribed circle), which tends to yield well‑spaced sampling
#' locations even in irregular polygons.
#'
#' @details
#' * The function assumes a projected (planar) CRS for domain so that
#'   distances/areas are in linear units. If your data are in longitude/latitude,
#'   reproject (e.g., to UTM) before calling.
#' * The number of initial seeds is n_quadrats * voronoi_seed_factor. Larger
#'   values explore the domain more densely and typically find more suitable
#'   cells, at the cost of extra computation.
#' * Suitability is decided by comparing each cell’s inscribed‑circle radius
#'   to half the requested quadrat diagonal:
#'   radius >= sqrt(width^2 + height^2) / 2.
#' * If fewer than n_quadrats suitable cells exist, all suitable cells are
#'   used and a warning is issued.
#' * Quadrats are axis‑aligned rectangles centered at the selected cell
#'   centers. Because the cell centers can be near each other, axis‑aligned
#'   rectangles may still overlap slightly in tight spaces—even when circles
#'   do not. If strict non‑overlap is required, post‑filtering is recommended.
#' * Internally uses sf::st_voronoi() for tessellation and
#'   st_inscribed_circle() (provided by lwgeom) to compute maximal
#'   in‑circle centers/radii.
#'
#' @param domain An sf polygon/multipolygon defining the sampling
#'   region. Must carry a valid projected CRS suitable for linear measurements.
#' @param n_quadrats Integer (≥1). Target number of quadrats to place.
#'   The actual number returned may be smaller if the domain cannot accommodate
#'   enough suitable Voronoi cells.
#' @param quadrat_size Numeric vector of length 2, c(width, height), giving
#'   quadrat dimensions in the same linear units as domain (e.g., meters).
#'   Both values must be positive and finite.
#' @param voronoi_seed_factor Numeric (≥1 recommended). Multiplier controlling
#'   how many random seeds are used to build the Voronoi diagram:
#'   n_seeds = n_quadrats * voronoi_seed_factor. Use larger values (e.g., 8–20)
#'   in very irregular or narrow domains to improve coverage.
#'
#' @return An sf object (polygons) with:
#' \describe{
#'   \item{quadrat_id}{Sequential integer identifier of the placed quadrats.}
#'   \item{geometry}{Axis‑aligned rectangular polygon for each quadrat.}
#' }
#'
#' @seealso
#' \code{\link{place_quadrats}} (random, non‑overlapping),
#' \code{\link{place_quadrats_tiled}} (systematic tiling),
#' \code{\link{place_quadrats_systematic}} (grid centers within domain),
#' \code{\link{place_quadrats_transect}} (parallel transects)
#'
#' @examples
#' \dontrun{
#' library(sf)
#' set.seed(2)
#' dom <- create_sampling_domain()
#' qs  <- place_quadrats_voronoi(
#'   domain = dom,
#'   n_quadrats = 20,
#'   quadrat_size = c(1.5, 1.5),
#'   voronoi_seed_factor = 12
#' )
#' plot(st_geometry(dom), col = “grey95”, border = “grey60”)
#' plot(st_geometry(qs), add = TRUE, border = “black”, lwd = 1)
#' }
place_quadrats_voronoi <- function(domain, n_quadrats, quadrat_size, voronoi_seed_factor) {
  n_seeds <- n_quadrats * voronoi_seed_factor
  domain_union <- sf::st_union(domain)
  seed_points <- sf::st_sample(domain_union, size = n_seeds, type = "random")
  voronoi_polys <- sf::st_voronoi(sf::st_union(seed_points))
  voronoi_clipped <- sf::st_intersection(sf::st_cast(voronoi_polys), domain_union)
  inscribed_circles <- suppressWarnings(sf::st_inscribed_circle(voronoi_clipped))
  quadrat_half_diag <- sqrt(quadrat_size[1]^2 + quadrat_size[2]^2) / 2
  radii <- sqrt(sf::st_area(inscribed_circles) / pi)
  suitable_indices <- which(radii >= quadrat_half_diag)
  num_possible <- length(suitable_indices)
  if (num_possible == 0) {
    warning("Voronoi placement failed: no suitable cells.")
    return(sf::st_sf(quadrat_id = integer(0), geometry = sf::st_sfc(crs = sf::st_crs(domain))))
  }
  if (num_possible < n_quadrats) {
    warning(sprintf("Voronoi placement found %d suitable locations; using all.", num_possible))
    n_quadrats <- num_possible
  }
  sampled_indices <- sample(suitable_indices, size = n_quadrats)
  final_centers <- sf::st_centroid(inscribed_circles[sampled_indices])
  quadrat_list <- lapply(sf::st_geometry(final_centers), function(pt) create_quadrat_from_center(pt, quadrat_size))
  final_quadrats_sfc <- sf::st_sfc(quadrat_list, crs = sf::st_crs(domain))
  sf::st_sf(quadrat_id = 1:length(final_quadrats_sfc), geometry = final_quadrats_sfc)
}
