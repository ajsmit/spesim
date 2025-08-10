#' Fast Strauss simulator on a polygon (fixed-n, MH relocation)
#'
#' @description
#' Simulate \eqn{n} points from a Strauss process *conditioned on \eqn{n}*
#' inside an \pkg{sf} polygon using a fast C++ Metropolis–Hastings kernel:
#' at each step, randomly pick a point and propose a new uniform location;
#' accept with probability \eqn{\min(1, \gamma^{\Delta s})}, where
#' \eqn{\Delta s} is the change in the number of neighbor pairs within radius \eqn{r}.
#' The implementation uses a spatial hash for \eqn{O(1)} expected neighbor lookups.
#'
#' @param domain An \pkg{sf} polygon/multipolygon defining the window.
#' @param n_target Integer; the fixed number of points.
#' @param r Interaction radius (map units).
#' @param gamma Inhibition in \eqn{(0,1]}; \eqn{<1} yields inhibition, \eqn{=1} reduces to CSR with relocation.
#' @param sweeps Total MH sweeps; each sweep proposes \eqn{\approx n} moves.
#' @param burnin Burn-in sweeps (currently not retained; final state is returned).
#' @param thin Thinning interval (not used to average—final state is returned).
#'
#' @return An \pkg{sf} POINT layer with columns \code{x}, \code{y} as geometry.
#' @export
simulate_points_strauss_fast <- function(
    domain, n_target, r, gamma, sweeps = 2000, burnin = 200, thin = 1) {
  bb <- sf::st_bbox(domain)
  mat <- rstrauss_bbox_cpp(
    n = as.integer(n_target),
    xmin = unname(bb["xmin"]), xmax = unname(bb["xmax"]),
    ymin = unname(bb["ymin"]), ymax = unname(bb["ymax"]),
    r = r, gamma = gamma,
    sweeps = sweeps, burnin = burnin, thin = thin
  )
  sf::st_as_sf(as.data.frame(mat), coords = c("x", "y"), crs = sf::st_crs(domain))
}
