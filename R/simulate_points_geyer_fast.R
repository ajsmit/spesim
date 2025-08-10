# R/simulate_points_geyer_fast.R

#' Fast Geyer saturation simulator (bbox Metropolis–Hastings)
#'
#' @param domain sf polygon/multipolygon defining the window (only bbox is used).
#' @param n_target integer, number of points to return.
#' @param r interaction radius.
#' @param gamma interaction parameter (<1 inhibition, >1 clustering).
#' @param sat saturation count (positive integer).
#' @param sweeps MH sweeps (per point).
#' @param burnin burn-in sweeps.
#' @param thin unused (reserved).
#'
#' @return \code{sf} POINT layer with \code{n_target} points.
#' @export
simulate_points_geyer_fast <- function(domain, n_target,
                                       r, gamma, sat,
                                       sweeps = 2000, burnin = 200, thin = 1) {
  bb <- sf::st_bbox(domain)
  out <- rgeyer_bbox_cpp(
    n_target = as.integer(n_target),
    xmin = unname(bb["xmin"]), xmax = unname(bb["xmax"]),
    ymin = unname(bb["ymin"]), ymax = unname(bb["ymax"]),
    r = as.numeric(r), gamma = as.numeric(gamma), sat = as.integer(sat),
    sweeps = as.integer(sweeps), burnin = as.integer(burnin), thin = as.integer(thin)
  )
  sfc <- sf::st_sfc(lapply(seq_along(out$x), function(i) sf::st_point(c(out$x[i], out$y[i]))),
    crs = sf::st_crs(domain)
  )
  sf::st_sf(geometry = sfc)
}
