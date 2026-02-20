# Internal helpers for the fast (bbox) point-process wrappers

# Clip an sf POINT layer to a (unioned) domain polygon.
.spesim_clip_points_to_domain <- function(pts_bbox, domain_union) {
  if (nrow(pts_bbox) == 0) return(pts_bbox)
  inside <- sf::st_within(pts_bbox, domain_union, sparse = TRUE)
  keep <- lengths(inside) > 0
  pts_bbox[keep, , drop = FALSE]
}

# Iterate bbox simulation + polygon clipping until we have >= n_target points.
#
# - gen_bbox(n_bbox) must return an sf POINT layer in the bbox CRS.
# - We never top-up by duplicating existing points.
# - If we still fall short, we fall back to uniform sampling inside the polygon.
.spesim_iterative_bbox_then_clip <- function(domain, n_target, frac,
                                            gen_bbox,
                                            oversample = 1.5,
                                            max_rounds = 6,
                                            max_bbox = 200000L) {
  stopifnot(inherits(domain, "sf"))
  n_target <- as.integer(n_target)
  if (n_target <= 0L) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = sf::st_crs(domain))))
  }

  frac <- as.numeric(frac)
  if (!is.finite(frac) || frac <= 0) frac <- 0.5

  domain_union <- sf::st_union(domain)
  pts_all <- sf::st_sf(geometry = sf::st_sfc(crs = sf::st_crs(domain)))

  dedup <- function(x) {
    if (nrow(x) <= 1) return(x)
    xy <- sf::st_coordinates(x)
    x[!duplicated(xy), , drop = FALSE]
  }

  need <- n_target
  for (k in seq_len(max_rounds)) {
    if (need <= 0L) break

    # Request enough bbox points so that, in expectation, we cover the polygon filter.
    n_bbox <- as.integer(ceiling((need / frac) * oversample))
    n_bbox <- max(1L, min(as.integer(max_bbox), n_bbox))

    pts_bbox <- gen_bbox(n_bbox)
    if (!inherits(pts_bbox, "sf")) stop("gen_bbox() must return an sf object")

    pts_in <- .spesim_clip_points_to_domain(pts_bbox, domain_union)

    if (nrow(pts_in) > 0) {
      pts_all <- dedup(rbind(pts_all, pts_in))
      need <- n_target - nrow(pts_all)
    }

    # modestly increase oversampling after a failure to reach target
    if (need > 0L) oversample <- oversample * 1.25
  }

  if (nrow(pts_all) < n_target) {
    # Fallback: uniform sampling inside polygon (no replacement).
    short <- n_target - nrow(pts_all)
    fill <- sf::st_sample(domain_union, size = short, type = "random")
    pts_all <- dedup(rbind(pts_all, sf::st_sf(geometry = fill)))
  }

  # Thin to exact n_target
  if (nrow(pts_all) > n_target) {
    pts_all <- pts_all[sample.int(nrow(pts_all), n_target), , drop = FALSE]
  }

  pts_all
}
