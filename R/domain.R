#' Create a random sampling domain polygon
#'
#' @description
#' Generates a single polygon intended to represent a realistic study area for
#' spatial simulation.
#'
#' The original domain generator in **spesim** produced an "organic" outline
#' derived from a radius–angle curve; while random, this tended to yield
#' variations on the same motif. This function now supports multiple *shape
#' families* so that default domains can look genuinely different between runs:
#'
#' - **Concave / irregular** (star-like, anisotropic, with indented bays)
#' - **Convex / irregular** (random convex hulls)
#' - **Regular convex** (circles/ellipses and boxy rectangles/squares)
#'
#' **Units:** The default domain uses an arbitrary planar coordinate system
#' (i.e., it is *unitless*). If you need real-world units (m, km) or a known CRS,
#' supply your own `sf` polygon domain to the simulation/plotting functions.
#'
#' @param shape Character. Which family of shapes to generate.
#'   - `"mixed"` (default): randomly chooses among concave, convex, ellipse, and
#'     rectangle-like shapes.
#'   - `"concave"`: irregular, often concave, optionally anisotropic.
#'   - `"convex"`: irregular but convex (convex hull of random points).
#'   - `"ellipse"`: oval/circular (superellipse with exponent 2).
#'   - `"rectangle"`: boxy rectangle/square (superellipse with large exponent).
#'
#' @param n_vertices Integer >= 12. Number of vertices used to trace the outline
#'   (more vertices => smoother shapes).
#' @param size Positive numeric. Controls the overall scale (roughly the
#'   semi-axis length for regular shapes, and base radius for irregular shapes).
#' @param aspect Optional numeric. Controls anisotropy as a *y/x* scale factor.
#'   - If `NULL` (default), an aspect ratio is drawn at random.
#'   - If a single number, it is used directly.
#'   - If length-2, a random value is drawn uniformly from that range.
#' @param rotation Optional numeric (radians). If `NULL`, a random rotation is
#'   used.
#' @param probs Numeric named vector of probabilities used when
#'   `shape = "mixed"`. Names must be a subset of `c("concave","convex","ellipse","rectangle")`.
#' @param attempts Integer. If a generated polygon is invalid (e.g.
#'   self-intersecting), retry up to this many times before falling back to a
#'   safe convex shape.
#'
#' @return An `sf` object with one row and a `POLYGON` geometry.
#'
#' @examples
#' \dontrun{
#' # default: a genuinely varied random domain
#' dom1 <- create_sampling_domain()
#' dom2 <- create_sampling_domain()
#'
#' # force a family
#' dom_c <- create_sampling_domain(shape = "concave")
#' dom_e <- create_sampling_domain(shape = "ellipse", aspect = 1) # circle
#' dom_r <- create_sampling_domain(shape = "rectangle", aspect = c(0.6, 1.8))
#' }
#'
#' @seealso [sf::st_polygon()], [sf::st_sf()]
#' @export
create_sampling_domain <- function(shape = c("mixed", "concave", "convex", "ellipse", "rectangle"),
                                  n_vertices = 80L,
                                  size = 10,
                                  aspect = NULL,
                                  rotation = NULL,
                                  probs = c(concave = 0.35, convex = 0.25, ellipse = 0.20, rectangle = 0.20),
                                  attempts = 20L) {
  shape <- match.arg(shape)

  n_vertices <- as.integer(n_vertices)
  if (!is.finite(n_vertices) || n_vertices < 12L) {
    stop("`n_vertices` must be an integer >= 12.")
  }
  if (!is.numeric(size) || length(size) != 1L || !is.finite(size) || size <= 0) {
    stop("`size` must be a single positive number.")
  }
  attempts <- as.integer(attempts)
  if (!is.finite(attempts) || attempts < 1L) {
    stop("`attempts` must be an integer >= 1.")
  }

  .sample_aspect <- function(aspect) {
    if (is.null(aspect)) {
      # mixture: sometimes near-isotropic, sometimes anisotropic
      if (runif(1) < 0.55) {
        return(runif(1, 0.9, 1.1))
      }
      return(exp(runif(1, log(0.5), log(2.0))))
    }
    if (!is.numeric(aspect) || any(!is.finite(aspect))) {
      stop("`aspect` must be NULL or finite numeric.")
    }
    if (length(aspect) == 1L) {
      if (aspect <= 0) stop("`aspect` must be > 0.")
      return(aspect)
    }
    if (length(aspect) == 2L) {
      lo <- min(aspect)
      hi <- max(aspect)
      if (lo <= 0) stop("`aspect` range must be > 0.")
      return(runif(1, lo, hi))
    }
    stop("`aspect` must be NULL, length-1, or length-2.")
  }

  .sample_rotation <- function(rotation) {
    if (is.null(rotation)) return(runif(1, 0, 2 * pi))
    if (!is.numeric(rotation) || length(rotation) != 1L || !is.finite(rotation)) {
      stop("`rotation` must be NULL or a single finite number (radians).")
    }
    rotation
  }

  .rotate <- function(xy, ang) {
    ca <- cos(ang)
    sa <- sin(ang)
    cbind(xy[, 1] * ca - xy[, 2] * sa, xy[, 1] * sa + xy[, 2] * ca)
  }

  .close_ring <- function(xy) {
    # sf requires the ring to be *exactly* closed (first == last).
    # (all.equal() is not sufficient here.)
    if (!(nrow(xy) >= 2L && all(xy[1, ] == xy[nrow(xy), ]))) {
      xy <- rbind(xy, xy[1, , drop = FALSE])
    }
    xy
  }

  .as_sf_polygon <- function(coords) {
    poly <- sf::st_polygon(list(coords))
    sf::st_sf(geometry = sf::st_sfc(poly))
  }

  .make_superellipse <- function(n, a, b, p, rot) {
    theta <- seq(0, 2 * pi, length.out = n)
    ct <- cos(theta)
    st <- sin(theta)
    # superellipse parameterization
    x <- a * sign(ct) * abs(ct)^(2 / p)
    y <- b * sign(st) * abs(st)^(2 / p)
    xy <- cbind(x, y)
    xy <- .rotate(xy, rot)
    .close_ring(xy)
  }

  .make_concave <- function(n, base_r, asp, rot) {
    theta <- seq(0, 2 * pi, length.out = n + 1L)
    theta <- theta[-length(theta)]

    # Random Fourier-like radial curve with occasional deep inlets.
    k <- sample(4:8, 1)
    freqs <- sample(2:14, k, replace = FALSE)
    phases <- runif(k, 0, 2 * pi)

    # amplitudes scaled to base radius; keep moderate to reduce self-intersection
    amps <- runif(k, 0.10, 0.35) * base_r

    r <- rep(base_r * runif(1, 0.9, 1.2), n)
    for (i in seq_len(k)) {
      r <- r + amps[i] * sin(freqs[i] * theta + phases[i])
    }

    # Add short-scale roughness
    r <- r + rnorm(n, sd = 0.06 * base_r)

    # Introduce a few "bays" (local reductions) to encourage concavity
    nbays <- sample(1:4, 1)
    bay_centres <- sample(seq_len(n), nbays)
    for (bc in bay_centres) {
      width <- sample(3:9, 1)
      idx <- ((bc - width):(bc + width) - 1L) %% n + 1L
      depth <- runif(1, 0.35, 0.70)
      taper <- exp(-((seq_along(idx) - (length(idx) + 1) / 2)^2) / (0.35 * length(idx))^2)
      r[idx] <- r[idx] * (1 - depth * (taper / max(taper)))
    }

    r <- pmax(r, 0.15 * base_r)

    x <- r * cos(theta)
    y <- r * sin(theta) * asp

    # light Cartesian jitter to break any remaining symmetry
    x <- x + rnorm(n, sd = 0.03 * base_r)
    y <- y + rnorm(n, sd = 0.03 * base_r)

    xy <- cbind(x, y)
    xy <- .rotate(xy, rot)
    .close_ring(xy)
  }

  .make_convex_hull <- function(base_r, asp, rot) {
    m <- max(150L, 3L * n_vertices)
    pts <- cbind(rnorm(m, sd = base_r * 0.65), rnorm(m, sd = base_r * 0.65) * asp)
    pts <- .rotate(pts, rot)

    # convex hull of a MULTIPOINT yields a polygon
    hull <- sf::st_convex_hull(sf::st_sfc(sf::st_multipoint(pts)))

    # st_convex_hull returns an sfc; wrap as sf
    sf::st_sf(geometry = hull)
  }

  # Choose a family if mixed
  if (shape == "mixed") {
    allowed <- c("concave", "convex", "ellipse", "rectangle")
    if (is.null(names(probs)) || any(!names(probs) %in% allowed)) {
      stop("When `shape = \"mixed\"`, `probs` must be a named vector with names in: ",
           paste(allowed, collapse = ", "))
    }
    pp <- probs / sum(probs)
    shape <- sample(names(pp), size = 1, prob = pp)
  }

  asp <- .sample_aspect(aspect)
  rot <- .sample_rotation(rotation)

  # Try multiple times to ensure a valid polygon for concave shapes.
  for (i in seq_len(attempts)) {
    if (shape == "concave") {
      coords <- .make_concave(n_vertices, base_r = size, asp = asp, rot = rot)
      dom <- .as_sf_polygon(coords)
    } else if (shape == "ellipse") {
      a <- size * runif(1, 0.85, 1.20)
      b <- a * asp
      coords <- .make_superellipse(n_vertices, a = a, b = b, p = 2, rot = rot)
      dom <- .as_sf_polygon(coords)
    } else if (shape == "rectangle") {
      a <- size * runif(1, 0.85, 1.20)
      b <- a * asp
      p <- runif(1, 8, 24) # higher => boxier
      coords <- .make_superellipse(n_vertices, a = a, b = b, p = p, rot = rot)
      dom <- .as_sf_polygon(coords)
    } else if (shape == "convex") {
      dom <- .make_convex_hull(base_r = size, asp = asp, rot = rot)
    } else {
      stop("Unknown `shape`: ", shape)
    }

    # Validate; concave shapes can occasionally self-intersect
    ok <- TRUE
    v <- try(sf::st_is_valid(dom), silent = TRUE)
    if (!inherits(v, "try-error")) {
      ok <- isTRUE(as.logical(v[1]))
    }

    if (ok) return(dom)

    # If invalid, randomise rotation/aspect slightly and retry.
    rot <- rot + runif(1, -0.5, 0.5)
    asp <- max(0.2, asp * exp(runif(1, -0.15, 0.15)))
  }

  # Fallback: always return a safe convex shape.
  a <- size
  b <- a * asp
  coords <- .make_superellipse(n_vertices, a = a, b = b, p = 10, rot = rot)
  .as_sf_polygon(coords)
}
