#' Parallel‑Transect Quadrat Placement
#'
#' @description
#' Places axis‑aligned rectangular quadrats at evenly spaced positions along
#' a set of parallel transect lines crossing the domain at a specified
#' compass angle. The function first shrinks the domain by half the quadrat
#' diagonal to form a safe sampling area, ensuring that every returned
#' quadrat lies fully inside the original domain.
#'
#' @details
#' * The transect orientation is given by angle in compass degrees
#'   (0 = North, 90 = East, 180 = South, 270 = West).
#' * n_transects parallel lines are positioned across the domain’s extent,
#'   then each line is cropped to the safe area. Along each retained segment,
#'   n_quadrats_per_transect centers are placed by linear interpolation
#'   from segment start to end (a single quadrat uses the midpoint).
#' * Each center generates an axis‑aligned rectangular quadrat of size
#'   quadrat_size = c(width, height) in the domain’s units.
#' * Use a projected (planar) CRS for domain (e.g., meters) so that
#'   distances and sizes are meaningful. If the domain is too small relative to
#'   quadrat_size, the safe area may be empty and no quadrats are returned.
#'
#' @param domain An sf polygon/multipolygon defining the sampling region.
#'   Must carry a valid projected CRS suitable for linear units (e.g., meters).
#' @param n_transects Integer (≥1). Number of parallel transect lines to generate.
#' @param n_quadrats_per_transect Integer (≥1). Number of quadrats to place
#'   along each transect after clipping to the safe area.
#' @param quadrat_size Numeric vector of length 2, c(width, height), giving
#'   quadrat dimensions in the same units as the domain CRS. Both values
#'   must be positive and finite.
#' @param angle Numeric. Compass bearing in degrees for transect orientation
#'   (0 = North/up, 90 = East/right). Values outside [0, 360) are allowed and
#'   are interpreted modulo 360.
#'
#' @return An sf object (polygons) with:
#' \describe{
#'   \item{quadrat_id}{Sequential integer identifier for each placed quadrat.}
#'   \item{geometry}{Axis‑aligned rectangular polygon for the quadrat.}
#' }
#'
#' @seealso
#' \code{\link{place_quadrats}} (random, non‑overlapping),
#' \code{\link{place_quadrats_systematic}} (regular grid over bounding box),
#' \code{\link{place_quadrats_tiled}} (systematic tiling by cell size),
#' \code{\link{place_quadrats_voronoi}} (Voronoi‑based centers)
#'
#' @examples
#' \dontrun{
#' library(sf)
#' dom <- create_sampling_domain()
#' qs <- place_quadrats_transect(
#'   domain = dom,
#'   n_transects = 3,
#'   n_quadrats_per_transect = 6,
#'   quadrat_size = c(1.5, 1.5),
#'   angle = 45
#' )
#' plot(st_geometry(dom), col = “grey95”, border = “grey60”)
#' plot(st_geometry(qs), add = TRUE, border = “black”, lwd = 1)
#' }
place_quadrats_transect <- function(domain, n_transects, n_quadrats_per_transect, quadrat_size, angle) {
  buffer_dist <- sqrt(quadrat_size[1]^2 + quadrat_size[2]^2) / 2
  safe_domain <- sf::st_buffer(domain, -buffer_dist)
  if (sf::st_is_empty(safe_domain) || sf::st_area(safe_domain) == 0) {
    warning("Domain too small for given quadrat size.")
    return(sf::st_sf(quadrat_id = integer(0), geometry = sf::st_sfc(crs = sf::st_crs(domain))))
  }

  bbox <- sf::st_bbox(domain)
  center <- sf::st_centroid(sf::st_as_sfc(bbox))
  diag_len <- sqrt((bbox["xmax"] - bbox["xmin"])^2 + (bbox["ymax"] - bbox["ymin"])^2) * 1.5
  y_coords <- seq(bbox["ymin"], bbox["ymax"], length.out = n_transects + 2)[2:(n_transects + 1)]

  math_angle_deg <- 90 - angle
  math_angle_rad <- math_angle_deg * pi / 180
  rot_matrix <- matrix(c(cos(math_angle_rad), sin(math_angle_rad), -sin(math_angle_rad), cos(math_angle_rad)), 2, 2)

  horizontal_lines <- sf::st_sfc(lapply(y_coords, function(y) {
    sf::st_linestring(matrix(c(
      sf::st_coordinates(center)[1] - diag_len / 2,
      sf::st_coordinates(center)[1] + diag_len / 2,
      y, y
    ), ncol = 2))
  }), crs = sf::st_crs(domain))
  rotated_lines <- (horizontal_lines - center) * rot_matrix + center
  sf::st_crs(rotated_lines) <- sf::st_crs(domain)

  points_list <- vector("list", n_transects)
  for (i in seq_len(n_transects)) {
    line <- rotated_lines[i]
    segment <- sf::st_intersection(line, safe_domain)
    if (sf::st_is_empty(segment) || sf::st_length(segment) == 0) next
    segment_coords <- sf::st_coordinates(segment)
    start_pt <- segment_coords[1, c("X", "Y")]
    end_pt <- segment_coords[nrow(segment_coords), c("X", "Y")]
    fractions <- if (n_quadrats_per_transect > 1) seq(0, 1, length.out = n_quadrats_per_transect) else 0.5
    x_coords <- start_pt[1] + fractions * (end_pt[1] - start_pt[1])
    y_coords <- start_pt[2] + fractions * (end_pt[2] - start_pt[2])
    points_list[[i]] <- sf::st_as_sf(data.frame(x = x_coords, y = y_coords), coords = c("x", "y"), crs = sf::st_crs(domain))
  }
  points_list <- Filter(Negate(is.null), points_list)
  if (length(points_list) == 0) {
    warning("No transects intersected the safe sampling area.")
    return(sf::st_sf(quadrat_id = integer(0), geometry = sf::st_sfc(crs = sf::st_crs(domain))))
  }
  candidate_centers <- do.call(rbind, points_list)
  quadrat_geometries <- sf::st_sfc(lapply(sf::st_geometry(candidate_centers), function(pt) create_quadrat_from_center(pt, quadrat_size)), crs = sf::st_crs(domain))
  sf::st_sf(quadrat_id = 1:length(quadrat_geometries), geometry = quadrat_geometries)
}
