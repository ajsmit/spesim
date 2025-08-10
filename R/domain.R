#' Create an irregular sampling domain polygon
#'
#' @description
#' Generates a single, closed, irregular polygon intended to represent a
#' realistic, non-rectangular study area for spatial simulation. The polygon
#' outline is derived from a parameterized radius-angle curve with added
#' sinusoidal perturbations and random noise to produce an "organic" shape.
#'
#' The output is returned as an \code{sf} object containing one polygon in its
#' geometry column. This can be used directly in plotting or as the spatial
#' domain within which individuals, quadrats, or transects are placed.
#'
#' @details
#' The shape is generated in polar coordinates \eqn{(r, \theta)} using:
#' \itemize{
#'   \item a base radius of 10 units,
#'   \item sinusoidal radial modulations with periods of 3 and 5 cycles per
#'         revolution,
#'   \item independent uniform noise on both radius and Cartesian coordinates
#'         to break symmetry,
#'   \item vertical compression by a factor of 0.8 to induce aspect variation.
#' }
#' The 20 vertices are connected in sequence and closed by repeating the first
#' point. The result is wrapped in \code{sf::st_polygon()} and returned as a
#' single-row \code{sf} object.
#'
#' @return An \code{sf} object with one row and a \code{POLYGON} geometry.
#'
#' @examples
#' \dontrun{
#' domain <- create_sampling_domain()
#' plot(domain$geometry)
#' }
#'
#' @seealso \code{\link[sf]{st_polygon}}, \code{\link[sf]{st_sf}}
#' @export
create_sampling_domain <- function() {
  theta <- seq(0, 2 * pi, length.out = 20)
  r <- 10 + 3 * sin(3 * theta) + 2 * cos(5 * theta) + runif(20, -1, 1)
  x <- r * cos(theta) + runif(20, -0.5, 0.5)
  y <- r * sin(theta) * 0.8 + runif(20, -0.5, 0.5)
  coords <- cbind(x, y)
  coords <- rbind(coords, coords[1, ])
  polygon <- sf::st_polygon(list(coords))
  sf::st_sf(geometry = sf::st_sfc(polygon))
}
