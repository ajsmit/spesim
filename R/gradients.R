#' Create synthetic environmental gradients over a spatial domain
#'
#' @description
#' Generates three continuous, spatially autocorrelated environmental
#' gradients — temperature, elevation, and rainfall — over a regular grid
#' covering the specified polygonal study area. The gradients are expressed
#' in both normalised \eqn{[0, 1]} form and scaled to physical units.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Derives the bounding box of \code{domain} and constructs a
#'         regular \code{resolution} × \code{resolution} grid of points.
#'   \item Normalises \eqn{x} and \eqn{y} coordinates to \eqn{[0, 1]}.
#'   \item Computes:
#'         \itemize{
#'           \item \strong{Temperature}: weighted linear combination of
#'                 normalised \eqn{x} (0.7) and \eqn{y} (0.3) coordinates.
#'           \item \strong{Elevation}: decreases radially from domain centre.
#'           \item \strong{Rainfall}: weighted difference of \eqn{x} and \eqn{y}
#'                 components, rescaled to \eqn{[0, 1]}.
#'         }
#'   \item Adds Gaussian noise with standard deviation \code{noise_level}
#'         to each raw gradient before clipping to \eqn{[0, 1]}.
#'   \item Produces scaled versions in common physical units:
#'         \code{temperature_C} (°C), \code{elevation_m} (metres),
#'         and \code{rainfall_mm} (millimetres/year).
#' }
#'
#' Note that gradients are computed for all grid points within the
#' \emph{bounding box} of \code{domain}; no masking to the polygon interior
#' is applied here. Masking can be done later with \code{\link[sf]{st_intersects}}
#' or similar.
#'
#' @param domain An \code{sf} polygon object defining the extent over which
#'   gradients will be computed. Only its bounding box is used in this step.
#' @param resolution Integer ≥ 2; number of grid steps per axis. The output
#'   will contain \code{resolution^2} rows.
#' @param noise_level Non-negative numeric; standard deviation of Gaussian
#'   noise added independently to each gradient before clipping.
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{\code{x, y}}{Grid point coordinates (same CRS as \code{domain}).}
#'     \item{\code{temperature, elevation, rainfall}}{Normalised gradients in \eqn{[0, 1]}.}
#'     \item{\code{temperature_C}}{Temperature in degrees Celsius (approx. range: -2 to 28).}
#'     \item{\code{elevation_m}}{Elevation in metres (0–2000).}
#'     \item{\code{rainfall_mm}}{Annual rainfall in millimetres (200–900).}
#'   }
#'
#' @examples
#' \dontrun{
#' domain <- create_sampling_domain()
#' env <- create_environmental_gradients(domain, resolution = 50, noise_level = 0.05)
#' head(env)
#' }
#'
#' @seealso \code{\link[sf]{st_bbox}}, \code{\link[sf]{st_intersects}}
#' @export
create_environmental_gradients <- function(domain, resolution, noise_level) {
  bbox <- sf::st_bbox(domain)
  x_seq <- seq(bbox["xmin"], bbox["xmax"], length.out = resolution)
  y_seq <- seq(bbox["ymin"], bbox["ymax"], length.out = resolution)
  grid <- expand.grid(x = x_seq, y = y_seq)

  x_norm <- (grid$x - bbox["xmin"]) / (bbox["xmax"] - bbox["xmin"])
  y_norm <- (grid$y - bbox["ymin"]) / (bbox["ymax"] - bbox["ymin"])

  temp_vals <- x_norm * 0.7 + y_norm * 0.3 + rnorm(nrow(grid), 0, noise_level)
  dist_center <- sqrt((x_norm - 0.5)^2 + (y_norm - 0.5)^2)
  elev_vals <- 1 - (dist_center / sqrt(0.5^2 + 0.5^2)) + rnorm(nrow(grid), 0, noise_level)
  rain_vals <- (-x_norm * 0.6 + y_norm * 0.8)
  rain_vals <- (rain_vals - min(rain_vals)) / (max(rain_vals) - min(rain_vals)) + rnorm(nrow(grid), 0, noise_level)

  grid$temperature <- pmax(0, pmin(1, temp_vals))
  grid$elevation <- pmax(0, pmin(1, elev_vals))
  grid$rainfall <- pmax(0, pmin(1, rain_vals))

  grid$temperature_C <- grid$temperature * 30 - 2
  grid$elevation_m <- grid$elevation * 2000
  grid$rainfall_mm <- grid$rainfall * 700 + 200
  grid
}
