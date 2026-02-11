#' Optional spatstat-based spatial diagnostics (edge-corrected)
#'
#' @description
#' These helpers are used by [spesim_audit()] when `diagnostics` includes
#' `"spatstat"`. They are kept optional to avoid a heavy dependency chain.
#'
#' The aim is not to replace full point-process analysis, but to provide
#' *edge-corrected* summaries (in the correct irregular window) that help users
#' verify that the realised spatial regime matches the named process they think
#' they asked for.
#'
#' @keywords internal

# Internal dependency check used by the optional audit layer.
#' @keywords internal
.spesim_require_spatstat <- function() {
  ok <- requireNamespace("spatstat.geom", quietly = TRUE) &&
    requireNamespace("spatstat.explore", quietly = TRUE)
  if (!ok) {
    stop(
      "spatstat diagnostics requested, but spatstat.geom/spatstat.explore are not installed. ",
      "Install with install.packages(c('spatstat.geom','spatstat.explore'))."
    )
  }
  invisible(TRUE)
}

# Convert an sf polygon/multipolygon to a spatstat owin.
# This uses spatstat.geom::as.owin.sf when available.
.sf_domain_to_owin <- function(domain) {
  .spesim_require_spatstat()

  # domain may be sf with multiple rows; union into one geometry
  dom_u <- sf::st_union(domain)
  dom_sf <- sf::st_as_sf(dom_u)

  # as.owin.sf is provided by spatstat.geom
  win <- spatstat.geom::as.owin(dom_sf)
  win
}

# Compute K and pcf summaries for a ppp, returning compact scalars.
.ppp_summaries <- function(X, rmax = NULL, nrval = 64) {
  .spesim_require_spatstat()

  if (is.null(rmax)) {
    # heuristic: 1/4 of min window side
    bb <- X$window$xrange
    ww <- diff(bb)
    hh <- diff(X$window$yrange)
    rmax <- 0.25 * min(ww, hh)
  }

  rvals <- seq(0, rmax, length.out = nrval)

  # Ripley's K (border edge correction is stable and fast)
  K <- spatstat.explore::Kest(X, correction = "border", r = rvals)

  # Summarise departure from CSR: L(r)-r
  L <- spatstat.explore::Lest(X, correction = "border", r = rvals)
  Lmr <- L$border - L$r

  # Pair correlation function (Ripley edge correction)
  g <- spatstat.explore::pcf(X, correction = "Ripley", r = rvals)

  # Compact scalars:
  # - mean(L(r)-r) over r>0 indicates clustering (>0) or inhibition (<0)
  # - max(g(r)) indicates peak clustering intensity; min(g(r)) indicates inhibition
  rr_ok <- is.finite(L$r) & L$r > 0
  g_ok <- is.finite(g$r) & g$r > 0

  out <- list(
    rmax = rmax,
    Lmr_mean = mean(Lmr[rr_ok], na.rm = TRUE),
    Lmr_max = max(Lmr[rr_ok], na.rm = TRUE),
    g_max = max(g$iso[g_ok], na.rm = TRUE),
    g_min = min(g$iso[g_ok], na.rm = TRUE)
  )

  out
}

#' Spatstat diagnostics for one or more species
#'
#' @description
#' Computes edge-corrected summary diagnostics (L(r)-r and pair-correlation
#' extrema) inside the true irregular domain window.
#'
#' This is primarily intended for method testing and teaching: it helps confirm
#' whether the realised pattern is clustered/inhibited/CSR-like, beyond a simple
#' nearest-neighbour heuristic.
#'
#' @param species_dist sf points with a `species` column.
#' @param domain sf polygon window.
#' @param species vector of species labels.
#' @param rmax maximum distance for summaries (optional).
#'
#' @return data.frame with one row per species.
#' @export
spesim_spatstat_diagnostics <- function(species_dist, domain, species = "all", rmax = NULL) {
  .spesim_require_spatstat()

  if (!inherits(species_dist, "sf")) stop("`species_dist` must be an sf object")
  if (!"species" %in% names(species_dist)) stop("`species_dist` must have a `species` column")

  spp_all <- sort(unique(as.character(species_dist$species)))
  spp <- if (identical(species, "all")) spp_all else intersect(as.character(species), spp_all)

  if (length(spp) == 0) return(data.frame())

  win <- .sf_domain_to_owin(domain)

  rows <- lapply(spp, function(s) {
    pts <- species_dist[species_dist$species == s, ]
    n <- nrow(pts)
    if (n < 5) {
      return(data.frame(
        species = s, n = n,
        rmax = NA_real_,
        Lmr_mean = NA_real_, Lmr_max = NA_real_,
        g_max = NA_real_, g_min = NA_real_,
        qualitative = "too_few_points"
      ))
    }

    xy <- sf::st_coordinates(pts)
    X <- spatstat.geom::ppp(x = xy[, 1], y = xy[, 2], window = win)

    sm <- .ppp_summaries(X, rmax = rmax)

    qual <- "CSR_like"
    if (is.finite(sm$Lmr_mean)) {
      if (sm$Lmr_mean > 0.05) qual <- "clustered" else if (sm$Lmr_mean < -0.05) qual <- "inhibited"
    }

    data.frame(
      species = s,
      n = n,
      rmax = sm$rmax,
      Lmr_mean = sm$Lmr_mean,
      Lmr_max = sm$Lmr_max,
      g_max = sm$g_max,
      g_min = sm$g_min,
      qualitative = qual
    )
  })

  do.call(rbind, rows)
}
