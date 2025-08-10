#' Fast Strauss simulator (Rcpp backend) returning an sf point layer
#'
#' @description
#' Fixed-\eqn{n} Metropolis–Hastings sampler for the Strauss process
#' with interaction radius \code{r} and parameter \code{gamma} (\eqn{0<\gamma\le 1}).
#' Uses a cell grid, cached neighbor counts, squared-distance checks, and
#' light auto-tuning for the proposal radius. This is a fast, self-contained
#' alternative to spatstat for baseline point generation.
#'
#' @param domain An \pkg{sf} polygon/multipolygon defining the window.
#' @param n_target Integer; number of points to return.
#' @param r Interaction radius (map units).
#' @param gamma Interaction parameter (\eqn{\le 1} for inhibition).
#' @param sweeps Number of MH sweeps (proposals per point). Default \code{2000}.
#' @param burnin Number of initial sweeps reserved for adaptation. Default \code{200}.
#' @param thin Currently unused (kept for API symmetry). Default \code{1}.
#' @param step0 Optional initial proposal radius (defaults to \code{0.5 * r}).
#' @param target_acc Target acceptance for auto-tuning. Default \code{0.35}.
#' @param tune_every Sweeps between tuning updates. Default \code{200}.
#'
#' @return An \pkg{sf} POINT layer with \code{n_target} rows.
#' @export
simulate_points_strauss_fast <- function(
    domain, n_target, r, gamma,
    sweeps = 2000, burnin = 200, thin = 1,
    step0 = NA_real_, target_acc = 0.35, tune_every = 200) {
  bb <- sf::st_bbox(domain)
  mat <- rstrauss_bbox_cpp(
    n = as.integer(n_target),
    xmin = unname(bb["xmin"]), xmax = unname(bb["xmax"]),
    ymin = unname(bb["ymin"]), ymax = unname(bb["ymax"]),
    r = as.numeric(r), gamma = as.numeric(gamma),
    sweeps = as.integer(sweeps),
    burnin = as.integer(burnin),
    thin = as.integer(thin),
    step0 = step0,
    target_acc = target_acc,
    tune_every = as.integer(tune_every)
  )
  sf::st_as_sf(as.data.frame(mat), coords = c("x", "y"), crs = sf::st_crs(domain))
}
