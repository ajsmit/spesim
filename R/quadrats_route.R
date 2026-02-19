#' Route-Based Quadrat Placement for Linearized Domains
#'
#' @description
#' Places quadrats along a linear route embedded in a polygonal domain. This is
#' intended for `DOMAIN_TYPE = "network"` and `"coastline"` workflows where
#' sampling should occur at known positions or equidistantly along a 1D path.
#'
#' @param domain An `sf` polygon/multipolygon sampling domain.
#' @param n_quadrats Integer target number of quadrats.
#' @param quadrat_size Numeric length-2 vector `c(width, height)`.
#' @param axis Character: `"x"` or `"y"`; axis used as the route coordinate.
#' @param mode Character: `"equidistant"` or `"specified"`.
#' @param positions Optional numeric vector in `[0,1]` used when
#'   `mode = "specified"`.
#' @param jitter_sd Numeric >= 0; perpendicular jitter around route axis.
#'
#' @return `sf` polygons with `quadrat_id`.
#' @export
place_quadrats_route <- function(domain,
                                 n_quadrats,
                                 quadrat_size,
                                 axis = "x",
                                 mode = c("equidistant", "specified"),
                                 positions = NULL,
                                 jitter_sd = 0) {
  mode <- match.arg(mode)
  axis <- tolower(as.character(axis %||% "x"))
  if (!axis %in% c("x", "y")) stop("axis must be 'x' or 'y'.")

  n_quadrats <- as.integer(n_quadrats)
  if (!is.finite(n_quadrats) || n_quadrats < 1L) stop("n_quadrats must be >= 1")

  jitter_sd <- as.numeric(jitter_sd)
  if (!is.finite(jitter_sd) || jitter_sd < 0) stop("jitter_sd must be >= 0")

  bb <- sf::st_bbox(domain)
  if (mode == "equidistant") {
    pos <- if (n_quadrats == 1L) 0.5 else seq(0, 1, length.out = n_quadrats)
  } else {
    pos <- as.numeric(positions)
    pos <- pos[is.finite(pos)]
    if (!length(pos)) stop("mode='specified' requires numeric positions in [0,1].")
    if (any(pos < 0 | pos > 1)) stop("positions must lie in [0,1].")
  }

  xmid <- (as.numeric(bb["xmin"]) + as.numeric(bb["xmax"])) / 2
  ymid <- (as.numeric(bb["ymin"]) + as.numeric(bb["ymax"])) / 2

  xy <- if (axis == "x") {
    cbind(
      as.numeric(bb["xmin"]) + pos * (as.numeric(bb["xmax"]) - as.numeric(bb["xmin"])),
      ymid + stats::rnorm(length(pos), 0, jitter_sd)
    )
  } else {
    cbind(
      xmid + stats::rnorm(length(pos), 0, jitter_sd),
      as.numeric(bb["ymin"]) + pos * (as.numeric(bb["ymax"]) - as.numeric(bb["ymin"]))
    )
  }

  centers <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]), coords = c("x", "y"), crs = sf::st_crs(domain))

  # Keep only centers whose resulting quadrats fully fit the domain.
  quads <- sf::st_sfc(lapply(sf::st_geometry(centers), function(pt) create_quadrat_from_center(pt, quadrat_size)),
                      crs = sf::st_crs(domain))
  within <- sf::st_within(quads, domain, sparse = FALSE)
  keep <- if (is.matrix(within)) as.logical(within[, 1]) else as.logical(within)

  quads <- quads[keep]
  out <- sf::st_sf(quadrat_id = seq_along(quads), geometry = quads)
  out
}
