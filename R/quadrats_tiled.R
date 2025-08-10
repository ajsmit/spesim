#' Tiled (Systematic Cell) Quadrat Placement
#'
#' @description
#' Places non-overlapping rectangular quadrats using a regular grid that
#' tiles the extent of the input domain. Candidate cells are generated with
#' sf::st_make_grid() using the requested quadrat width and height. Cells that
#' fall entirely inside the domain (via sf::st_within(..., sparse = FALSE))
#' are kept as valid locations, from which up to n_quadrats are sampled at
#' random (set a seed for reproducibility). Because candidates come from a grid,
#' selected quadrats never overlap.
#'
#' @details
#' - The function assumes a planar (projected) CRS where distances and areas
#'   are in linear units; if your domain is in longitude/latitude, reproject it
#'   first (e.g., to UTM).
#' - quadrat_size is passed to sf::st_make_grid(cellsize = ...), so the two
#'   numbers are interpreted as width and height in the domain’s units.
#' - If fewer than n_quadrats valid cells fit entirely inside the domain,
#'   the function places as many as possible and emits a warning.
#' - Returned quadrats are labeled sequentially in the order they are sampled
#'   (quadrat_id = 1, 2, …, n_quadrats).
#'
#' @param domain An sf polygon or multipolygon object defining the
#'   sampling region. Must have a valid CRS suitable for linear measurements
#'   (projected, not geographic).
#' @param n_quadrats Integer (≥1). The target number of quadrats to place.
#'   If fewer valid grid cells exist inside domain, the function will place
#'   the maximum possible and adjust this value downward with a warning.
#' @param quadrat_size Numeric vector of length 2, c(width, height), giving
#'   the quadrat dimensions in the same units as domain (e.g., meters).
#'   Both values must be positive and finite.
#'
#' @return An sf object (polygons) with two columns:
#' \describe{
#'   \item{quadrat_id}{Sequential integer identifier of the placed quadrats.}
#'   \item{geometry}{Polygon geometry of each quadrat.}
#' }
#'
#' @seealso
#' \code{\link{place_quadrats}} (random non-overlapping),
#' \code{\link{place_quadrats_systematic}} (grid centers within domain),
#' \code{\link{place_quadrats_transect}} (along parallel transects),
#' \code{\link{place_quadrats_voronoi}} (via Voronoi cells)
#'
#' @examples
#' \dontrun{
#' library(sf)
#' set.seed(1)
#' dom <- create_sampling_domain()
#' qs  <- place_quadrats_tiled(dom, n_quadrats = 12, quadrat_size = c(1.5, 1.5))
#' plot(st_geometry(dom), col = “grey95”, border = “grey60”)
#' plot(st_geometry(qs), add = TRUE, border = “black”, lwd = 1)
#' }
#'
#' @export
place_quadrats_tiled <- function(domain, n_quadrats, quadrat_size) {
  candidate_grid <- sf::st_make_grid(domain, cellsize = quadrat_size, what = "polygons")
  within_mat <- sf::st_within(candidate_grid, domain, sparse = FALSE)
  inside <- if (is.matrix(within_mat)) drop(within_mat[, 1, drop = TRUE]) else as.logical(within_mat)
  valid_locations <- candidate_grid[inside]
  num_possible <- length(valid_locations)
  if (num_possible == 0) {
    warning("Systematic placement failed: no quadrats fit inside the domain.")
    return(sf::st_sf(quadrat_id = integer(0), geometry = sf::st_sfc(crs = sf::st_crs(domain))))
  }
  if (num_possible < n_quadrats) {
    warning(sprintf("Could only place %d of %d requested quadrats.", num_possible, n_quadrats))
    n_quadrats <- num_possible
  }
  sampled_indices <- sample.int(num_possible, size = n_quadrats)
  final_quadrats_sfc <- valid_locations[sampled_indices]
  sf::st_sf(quadrat_id = seq_len(n_quadrats), geometry = final_quadrats_sfc)
}
