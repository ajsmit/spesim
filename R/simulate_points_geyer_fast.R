# R/simulate_points_geyer_fast.R

#' Fast Geyer saturation simulator (Rcpp backend; polygon-clipped)
#'
#' @description
#' Runs a fixed-n Metropolis–Hastings sampler in the domain bounding box and
#' then clips points to the polygon.
#'
#' This implementation is **heuristic but internally consistent** (see C++
#' header): it targets a saturation potential based on per-point neighbour counts.
#' It is **not** a drop-in reproduction of the exact 
#' \pkg{spatstat} Geyer saturation process likelihood.
#'
#' The result is thinned to return exactly `n_target` points. If clipping leaves
#' too few points, we generate additional points (never by duplicating existing ones).
#'
#' @param domain sf polygon/multipolygon defining the window.
#' @param n_target integer, number of points to return.
#' @param r interaction radius.
#' @param gamma interaction parameter (<1 inhibition, >1 clustering).
#' @param sat saturation count (positive integer).
#' @param sweeps MH sweeps (per point).
#' @param burnin burn-in sweeps.
#' @param thin unused (reserved).
#'
#' @return `sf` POINT layer with `n_target` points.
#' @export
simulate_points_geyer_fast <- function(domain, n_target,
                                       r, gamma, sat,
                                       sweeps = 2000, burnin = 200, thin = 1) {
  stopifnot(inherits(domain, "sf"))
  n_target <- as.integer(n_target)
  if (!is.finite(n_target) || n_target < 0L) stop("n_target must be a non-negative integer")
  if (n_target == 0L) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = sf::st_crs(domain))))
  }

  bb <- sf::st_bbox(domain)
  A_B <- as.numeric((bb["xmax"] - bb["xmin"]) * (bb["ymax"] - bb["ymin"]))
  A_D <- as.numeric(sf::st_area(sf::st_union(domain)))
  frac <- as.numeric(A_D / A_B)
  frac <- if (is.finite(frac) && frac > 0) frac else 0.5

  gen_bbox <- function(n_bbox) {
    out <- rgeyer_bbox_cpp(
      n_target = as.integer(n_bbox),
      xmin = unname(bb["xmin"]), xmax = unname(bb["xmax"]),
      ymin = unname(bb["ymin"]), ymax = unname(bb["ymax"]),
      r = as.numeric(r), gamma = as.numeric(gamma), sat = as.integer(sat),
      sweeps = as.integer(sweeps), burnin = as.integer(burnin), thin = as.integer(thin)
    )

    sfc <- sf::st_sfc(
      lapply(seq_along(out$x), function(i) sf::st_point(c(out$x[i], out$y[i]))),
      crs = sf::st_crs(domain)
    )
    sf::st_sf(geometry = sfc)
  }

  .spesim_iterative_bbox_then_clip(
    domain = domain,
    n_target = n_target,
    frac = frac,
    gen_bbox = gen_bbox,
    oversample = 1.5
  )
}
