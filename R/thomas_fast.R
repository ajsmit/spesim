#' Fast Thomas Process (Rcpp-backed) in an arbitrary polygon
#'
#' @description
#' Simulate \eqn{n\_target} points from a Thomas (Gaussian Neyman–Scott) process,
#' using a fast C++ generator in the domain's bounding box and then filtering to
#' the polygon with \pkg{sf}. This is a drop-in replacement for the spatstat-based
#' Thomas simulator and is substantially faster for large simulations.
#'
#' @param domain An \pkg{sf} polygon/multipolygon defining the study area (projected CRS recommended).
#' @param n_target Integer, desired number of retained points inside \code{domain}.
#' @param kappa Optional parent intensity (parents per unit area). If \code{NULL},
#'   it is derived from \code{n_target}, \code{mu}, and domain area.
#' @param mu Mean offspring per parent (Poisson).
#' @param sigma Gaussian cluster scale (standard deviation of offspring displacement; map units).
#' @param oversample Scalar \eqn{> 1}. Multiplier used to request enough raw
#'   offspring in the bbox so that, after polygon filtering, you still retain
#'   about \code{n_target}. Default \code{1.3}.
#'
#' @return An \pkg{sf} POINT layer with exactly \code{n_target} rows (if feasible),
#'   or fewer if \code{kappa}/\code{mu} are too small to generate enough points.
#'
#' @details
#' Let \eqn{A_D} be the polygon area and \eqn{A_B} its bbox area. The expected fraction
#' of bbox points that survive the polygon filter is \eqn{f \approx A_D / A_B}. We request roughly
#' \code{n_target / f} points from C++ (times \code{oversample}), then filter, and finally
#' thin or (if short) resample with replacement to return exactly \code{n_target}. If you prefer
#' strict “no replacement”, set \code{oversample} higher or supply a larger \code{kappa}.
#'
#' @examples
#' \dontrun{
#' dom <- create_sampling_domain()
#' pts <- rthomas_fast(dom, n_target = 2000, mu = 10, sigma = 1)
#' plot(sf::st_geometry(pts), pch = 16, cex = 0.4)
#' }
#' @export
rthomas_fast <- function(domain,
                         n_target,
                         kappa = NULL,
                         mu = 10,
                         sigma = 1,
                         oversample = 1.3) {
  stopifnot(inherits(domain, "sf"))
  if (n_target <= 0) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = sf::st_crs(domain))))
  }

  bb <- sf::st_bbox(domain)
  A_B <- (bb["xmax"] - bb["xmin"]) * (bb["ymax"] - bb["ymin"])
  A_D <- as.numeric(sf::st_area(sf::st_union(domain)))
  frac <- as.numeric(A_D / A_B)
  frac <- if (is.finite(frac) && frac > 0) frac else 0.5

  # If kappa not supplied, choose one to *target* n_target inside the polygon
  if (is.null(kappa)) {
    # Expected inside count ≈ kappa * A_D * mu
    # So set kappa ≈ n_target / (A_D * mu)
    kappa <- max(1e-8, n_target / (A_D * mu))
  }

  # Ask C++ for enough to cover polygon filtering
  max_points <- ceiling((n_target / frac) * oversample)

  xy <- rthomas_bbox_cpp(
    kappa = kappa,
    mu = mu,
    sigma = sigma,
    xmin = bb["xmin"],
    ymin = bb["ymin"],
    xmax = bb["xmax"],
    ymax = bb["ymax"],
    max_points = max_points
  )

  if (nrow(xy) == 0) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = sf::st_crs(domain))))
  }

  pts_bbox <- sf::st_as_sf(as.data.frame(xy), coords = c("X1", "X2"), crs = sf::st_crs(domain))
  # Keep only those inside polygon
  inside <- sf::st_within(pts_bbox, sf::st_union(domain), sparse = TRUE)
  keep <- lengths(inside) > 0
  pts_in <- pts_bbox[keep, , drop = FALSE]

  if (nrow(pts_in) == 0) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = sf::st_crs(domain))))
  }

  # Thin or top-up to EXACT n_target
  if (nrow(pts_in) >= n_target) {
    pts_in[sample.int(nrow(pts_in), n_target), , drop = FALSE]
  } else {
    # Not enough — resample with replacement to reach n_target
    idx <- c(seq_len(nrow(pts_in)), sample.int(nrow(pts_in), n_target - nrow(pts_in), replace = TRUE))
    pts_in[idx, , drop = FALSE]
  }
}
