#' Build network edges from ordered sites
#'
#' @description
#' Convenience helper to create a simple chain network from site coordinates by
#' ordering sites along a variable (e.g., river-km, coastline order).
#'
#' @param site_coords Data frame with at least `site`.
#' @param order_by Optional numeric vector or column name in `site_coords` used
#'   to order sites.
#' @param wrap Logical; connect last site back to first.
#' @param directed Logical; whether to treat edges as directed.
#'
#' @return Data frame with columns `from`, `to`, `weight`.
#' @export
build_network_edges <- function(site_coords,
                                order_by = NULL,
                                wrap = FALSE,
                                directed = FALSE) {
  if (is.null(site_coords) || !is.data.frame(site_coords)) stop("site_coords must be a data.frame")
  if (!"site" %in% names(site_coords)) stop("site_coords must contain column 'site'")
  if (nrow(site_coords) < 2) stop("Need at least two sites to build edges")

  ord <- seq_len(nrow(site_coords))
  if (!is.null(order_by)) {
    ov <- if (length(order_by) == 1 && is.character(order_by) && order_by %in% names(site_coords)) {
      site_coords[[order_by]]
    } else {
      as.numeric(order_by)
    }
    if (length(ov) != nrow(site_coords)) stop("order_by must have length nrow(site_coords) or name a column")
    ord <- order(as.numeric(ov), na.last = TRUE)
  }
  sc <- site_coords[ord, , drop = FALSE]
  if (!all(c("x", "y") %in% names(sc))) stop("site_coords must contain x and y")

  from <- as.character(sc$site[-nrow(sc)])
  to <- as.character(sc$site[-1])
  w <- sqrt(diff(as.numeric(sc$x))^2 + diff(as.numeric(sc$y))^2)

  if (isTRUE(wrap)) {
    from <- c(from, as.character(sc$site[nrow(sc)]))
    to <- c(to, as.character(sc$site[1]))
    w <- c(w, sqrt((as.numeric(sc$x[nrow(sc)]) - as.numeric(sc$x[1]))^2 +
                   (as.numeric(sc$y[nrow(sc)]) - as.numeric(sc$y[1]))^2))
  }

  out <- data.frame(from = from, to = to, weight = as.numeric(w), stringsAsFactors = FALSE)
  if (!isTRUE(directed)) {
    out <- unique(rbind(out, data.frame(from = out$to, to = out$from, weight = out$weight, stringsAsFactors = FALSE)))
  }
  out
}


# Internal: all-pairs shortest paths on a weighted graph (base-R Dijkstra).
.compute_graph_distances <- function(site_ids, edges, directed = FALSE) {
  ids <- as.character(site_ids)
  n <- length(ids)
  if (n == 0) return(matrix(numeric(0), 0, 0))
  idx <- stats::setNames(seq_len(n), ids)
  D <- matrix(Inf, n, n, dimnames = list(ids, ids))
  diag(D) <- 0

  if (!is.null(edges) && nrow(edges) > 0) {
    e <- edges
    if (!all(c("from", "to") %in% names(e))) stop("edges must contain from, to")
    if (!"weight" %in% names(e)) e$weight <- 1
    e$from <- as.character(e$from)
    e$to <- as.character(e$to)
    e$weight <- as.numeric(e$weight)
    keep <- is.finite(e$weight) & e$weight >= 0 & e$from %in% ids & e$to %in% ids
    e <- e[keep, , drop = FALSE]
    for (i in seq_len(nrow(e))) {
      a <- idx[[e$from[i]]]
      b <- idx[[e$to[i]]]
      D[a, b] <- min(D[a, b], e$weight[i])
      if (!isTRUE(directed)) D[b, a] <- min(D[b, a], e$weight[i])
    }
  }

  dijkstra <- function(src) {
    dist <- rep(Inf, n)
    vis <- rep(FALSE, n)
    dist[src] <- 0
    for (k in seq_len(n)) {
      cand <- which(!vis)
      if (!length(cand)) break
      u <- cand[which.min(dist[cand])]
      if (!is.finite(dist[u])) break
      vis[u] <- TRUE
      nbr <- which(is.finite(D[u, ]) & D[u, ] > 0)
      if (!length(nbr)) next
      alt <- dist[u] + D[u, nbr]
      better <- alt < dist[nbr]
      if (any(better)) dist[nbr[better]] <- alt[better]
    }
    dist
  }

  out <- matrix(Inf, n, n, dimnames = list(ids, ids))
  for (i in seq_len(n)) out[i, ] <- dijkstra(i)
  diag(out) <- 0
  out
}


#' Build a `spesim_result` directly from observed data
#'
#' @description
#' Creates a full `spesim_result` object using real observed site x species
#' data, site coordinates, and optional site-level environment. No simulation is
#' performed.
#'
#' @param abund_matrix Data frame with first column `site` and remaining species
#'   abundance columns.
#' @param site_coords Data frame with columns `site`, `x`, `y`.
#' @param site_env Optional data frame with `site` and environmental columns.
#' @param domain Optional `sf` polygon domain. If `NULL`, a convex hull + small
#'   buffer is built from `site_coords`.
#' @param env_gradients Optional environmental grid. If `NULL`, a point grid is
#'   created from `site_env` joined to `site_coords`.
#' @param edges Optional edge table (`from`, `to`, optional `weight`) defining
#'   connectivity among sites.
#' @param order_by Optional ordering used to build chain edges when `edges` is
#'   `NULL`.
#' @param directed Logical; treat connectivity as directed.
#' @param wrap Logical; only used when chain edges are auto-built.
#' @param point_jitter Numeric jitter for reconstructed `species_dist` points.
#'
#' @return A `spesim_result`.
#' @export
spesim_from_observed <- function(abund_matrix,
                                 site_coords,
                                 site_env = NULL,
                                 domain = NULL,
                                 env_gradients = NULL,
                                 edges = NULL,
                                 order_by = NULL,
                                 directed = FALSE,
                                 wrap = FALSE,
                                 point_jitter = 0) {
  if (!is.data.frame(abund_matrix) || !"site" %in% names(abund_matrix)) {
    stop("abund_matrix must be a data.frame with first column 'site'")
  }
  if (!is.data.frame(site_coords) || !all(c("site", "x", "y") %in% names(site_coords))) {
    stop("site_coords must contain columns: site, x, y")
  }
  if (!all(as.character(abund_matrix$site) %in% as.character(site_coords$site))) {
    stop("All abund_matrix$site values must exist in site_coords$site")
  }

  site_coords <- site_coords[match(as.character(abund_matrix$site), as.character(site_coords$site)), , drop = FALSE]
  if (!is.null(site_env)) {
    if (!("site" %in% names(site_env))) stop("site_env must contain column 'site'")
    site_env <- site_env[match(as.character(abund_matrix$site), as.character(site_env$site)), , drop = FALSE]
  } else {
    site_env <- data.frame(site = abund_matrix$site, stringsAsFactors = FALSE)
  }

  if (is.null(domain)) {
    pts <- sf::st_as_sf(site_coords, coords = c("x", "y"), crs = NA)
    hull <- sf::st_convex_hull(sf::st_union(pts))
    bb <- sf::st_bbox(hull)
    d <- sqrt((bb["xmax"] - bb["xmin"])^2 + (bb["ymax"] - bb["ymin"])^2)
    domain <- sf::st_as_sf(sf::st_buffer(hull, d * 0.05))
  }

  if (is.null(edges) && (!is.null(order_by) || "linear_pos" %in% names(site_coords))) {
    ob <- if (!is.null(order_by)) order_by else "linear_pos"
    edges <- build_network_edges(site_coords, order_by = ob, wrap = wrap, directed = directed)
  }
  if (!is.null(edges)) {
    gdist <- .compute_graph_distances(site_coords$site, edges, directed = directed)
    attr(site_coords, "graph_dist_matrix") <- gdist
    attr(site_coords, "graph_directed") <- isTRUE(directed)
  }

  if (is.null(env_gradients)) {
    env_gradients <- dplyr::left_join(
      site_coords[, c("x", "y", "site"), drop = FALSE],
      site_env,
      by = "site"
    )
  }

  # Build tiny square quadrats around each site point.
  dx <- stats::sd(as.numeric(site_coords$x), na.rm = TRUE)
  dy <- stats::sd(as.numeric(site_coords$y), na.rm = TRUE)
  w <- if (is.finite(dx) && dx > 0) dx * 0.03 else 0.01
  h <- if (is.finite(dy) && dy > 0) dy * 0.03 else 0.01
  mk_quad <- function(x, y) {
    sf::st_polygon(list(matrix(c(
      x - w, y - h,
      x + w, y - h,
      x + w, y + h,
      x - w, y + h,
      x - w, y - h
    ), ncol = 2, byrow = TRUE)))
  }
  quads <- sf::st_sf(
    quadrat_id = site_coords$site,
    geometry = sf::st_sfc(lapply(seq_len(nrow(site_coords)), function(i) mk_quad(site_coords$x[i], site_coords$y[i])),
      crs = sf::st_crs(domain))
  )

  # Reconstruct species points from abundances at site locations.
  spp <- setdiff(names(abund_matrix), "site")
  rows <- list()
  k <- 0L
  for (i in seq_len(nrow(abund_matrix))) {
    for (sp in spp) {
      n <- as.integer(abund_matrix[[sp]][i])
      if (!is.finite(n) || n <= 0) next
      x0 <- rep(as.numeric(site_coords$x[i]), n)
      y0 <- rep(as.numeric(site_coords$y[i]), n)
      if (point_jitter > 0) {
        x0 <- x0 + stats::rnorm(n, 0, point_jitter)
        y0 <- y0 + stats::rnorm(n, 0, point_jitter)
      }
      k <- k + 1L
      rows[[k]] <- data.frame(species = rep(sp, n), x = x0, y = y0, stringsAsFactors = FALSE)
    }
  }
  spp_df <- if (length(rows)) do.call(rbind, rows) else data.frame(species = character(0), x = numeric(0), y = numeric(0))
  species_dist <- sf::st_as_sf(spp_df, coords = c("x", "y"), crs = sf::st_crs(domain))

  P <- list(
    SEED = NA_integer_,
    N_SPECIES = length(spp),
    N_INDIVIDUALS = sum(as.matrix(abund_matrix[, spp, drop = FALSE]), na.rm = TRUE),
    SAMPLING_SCHEME = "observed",
    N_QUADRATS = nrow(site_coords),
    DOMAIN_TYPE = if (!is.null(edges)) "network" else "polygon",
    LINEAR_WRAP = isTRUE(wrap)
  )

  res <- list(
    P = P,
    domain = domain,
    species_dist = species_dist,
    quadrats = quads,
    env_gradients = env_gradients,
    abund_matrix = abund_matrix,
    site_env = site_env,
    site_coords = site_coords
  )
  new_spesim_result(res)
}


#' Simulate unsampled sites between observed sites
#'
#' @description
#' Generates pseudo-sites between observed sites by interpolating local
#' abundance structure and optional environment along `linear_pos`.
#'
#' @param res A `spesim_result` (typically from [spesim_from_observed()]).
#' @param n_new_sites Integer number of pseudo-sites to generate.
#' @param distance_power IDW power for weighting neighbouring observed sites.
#' @param direction_bias Numeric in \eqn{[-1,1]}; positive favors downstream/higher
#'   `linear_pos`, negative favors upstream/lower `linear_pos`.
#'
#' @return A list with `site_coords`, `abund_matrix`, and `site_env` for the new
#'   pseudo-sites.
#' @export
simulate_between_observed <- function(res,
                                      n_new_sites = 20,
                                      distance_power = 2,
                                      direction_bias = 0) {
  if (is.null(res$abund_matrix) || is.null(res$site_coords)) {
    stop("res must contain abund_matrix and site_coords")
  }
  A <- res$abund_matrix
  C <- res$site_coords
  if (!("linear_pos" %in% names(C))) {
    stop("site_coords must contain linear_pos for interpolation")
  }
  n_new_sites <- as.integer(n_new_sites)
  if (!is.finite(n_new_sites) || n_new_sites < 1) stop("n_new_sites must be >= 1")
  pwr <- as.numeric(distance_power)
  if (!is.finite(pwr) || pwr <= 0) pwr <- 2
  b <- as.numeric(direction_bias)
  if (!is.finite(b)) b <- 0
  b <- max(-1, min(1, b))

  spp <- setdiff(names(A), "site")
  A_num <- as.matrix(A[, spp, drop = FALSE])
  lp_obs <- as.numeric(C$linear_pos)
  ord <- order(lp_obs)
  lp_obs <- lp_obs[ord]
  C <- C[ord, , drop = FALSE]
  A_num <- A_num[ord, , drop = FALSE]
  tot_obs <- rowSums(A_num)

  lp_new <- seq(min(lp_obs), max(lp_obs), length.out = n_new_sites + 2)
  lp_new <- lp_new[-c(1, length(lp_new))]
  x_new <- stats::approx(lp_obs, as.numeric(C$x), xout = lp_new, rule = 2)$y
  y_new <- stats::approx(lp_obs, as.numeric(C$y), xout = lp_new, rule = 2)$y
  t_new <- pmax(1, round(stats::approx(lp_obs, as.numeric(tot_obs), xout = lp_new, rule = 2)$y))

  prob_for_pos <- function(p) {
    d <- abs(lp_obs - p)
    w <- 1 / ((d + 1e-6)^pwr)
    if (b != 0) {
      down <- lp_obs >= p
      up <- !down
      w[down] <- w[down] * (1 + max(0, b))
      w[up] <- w[up] * (1 + max(0, -b))
    }
    w <- w / sum(w)
    p_sp <- colSums(A_num * w)
    if (!all(is.finite(p_sp)) || sum(p_sp) <= 0) p_sp <- rep(1, length(spp))
    p_sp / sum(p_sp)
  }

  A_new <- matrix(0, nrow = n_new_sites, ncol = length(spp), dimnames = list(NULL, spp))
  for (i in seq_len(n_new_sites)) {
    p_sp <- prob_for_pos(lp_new[i])
    A_new[i, ] <- as.integer(stats::rmultinom(1, size = t_new[i], prob = p_sp)[, 1])
  }

  site_new <- paste0("interp_", seq_len(n_new_sites))
  out_coords <- data.frame(site = site_new, x = x_new, y = y_new, linear_pos = lp_new, stringsAsFactors = FALSE)
  out_abund <- data.frame(site = site_new, as.data.frame(A_new), check.names = FALSE)

  out_env <- data.frame(site = site_new, stringsAsFactors = FALSE)
  if (!is.null(res$site_env) && is.data.frame(res$site_env) && "site" %in% names(res$site_env)) {
    E <- res$site_env
    E <- E[match(C$site, E$site), , drop = FALSE]
    num_cols <- names(E)[vapply(E, is.numeric, logical(1))]
    num_cols <- setdiff(num_cols, "site")
    for (nm in num_cols) {
      out_env[[nm]] <- stats::approx(lp_obs, as.numeric(E[[nm]]), xout = lp_new, rule = 2)$y
    }
  }

  list(site_coords = out_coords, abund_matrix = out_abund, site_env = out_env)
}
