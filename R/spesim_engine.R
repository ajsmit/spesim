#' Internal simulation engine (no I/O)
#'
#' @description
#' `spesim_engine()` runs the core simulation given a fully-resolved parameter
#' list `P`, a study-area `domain`, and (optionally) pre-resolved interaction
#' fields inside `P`.
#'
#' This function is intentionally **pure** with respect to I/O:
#' it does not read configuration files, does not create timestamped prefixes,
#' and does not write CSVs/figures/reports. It returns a plain list with the
#' standard simulation components.
#'
#' The public orchestrator [spesim_run()] is responsible for:
#' - reading init files (via [load_config()])
#' - creating a default domain (possibly varying across runs)
#' - resolving interspecific interactions (file pointer / inline rules)
#' - optional writing of outputs
#'
#' @param P Fully-resolved parameter list. Must include at least:
#'   `N_SPECIES`, `N_INDIVIDUALS`, `SAMPLING_SCHEME`, `N_QUADRATS`,
#'   `QUADRAT_SIZE`, `SAMPLING_RESOLUTION`, and `ENVIRONMENTAL_NOISE`.
#'   If interactions are used, `INTERACTION_RADIUS` and `INTERACTION_MATRIX`
#'   should be present.
#' @param domain An `sf` polygon study area.
#'
#' @return A named list with at least:
#' - `P`
#' - `domain`
#' - `species_dist`
#' - `quadrats`
#' - `env_gradients`
#' - `abund_matrix`
#' - `site_env`
#' - `site_coords`
#'
#' @keywords internal
#' @noRd
spesim_engine <- function(P, domain) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (is.null(P) || !is.list(P)) stop("`P` must be a list.")
  if (is.null(domain) || !inherits(domain, "sf")) stop("`domain` must be an sf object.")

  # Reproducibility: always seed before any randomness (best-effort)
  if (!is.null(P$SEED) && is.finite(P$SEED)) {
    set.seed(as.integer(P$SEED))
  }

  # Domain is provided by caller (spesim_run); do NOT generate one here.

  # --- Core simulation -------------------------------------------------------
  if (!is.null(P$SEED) && is.finite(P$SEED)) {
    set.seed(as.integer(P$SEED) + 1L)
  }
  species_dist <- generate_heterogeneous_distribution(domain, P)

  if (!is.null(P$SEED) && is.finite(P$SEED)) {
    set.seed(as.integer(P$SEED) + 2L)
  }
  quadrats <- switch(P$SAMPLING_SCHEME,
    "random"     = place_quadrats(domain, P$N_QUADRATS, P$QUADRAT_SIZE),
    "tiled"      = place_quadrats_tiled(domain, P$N_QUADRATS, P$QUADRAT_SIZE),
    "systematic" = place_quadrats_systematic(domain, P$N_QUADRATS, P$QUADRAT_SIZE),
    "transect"   = place_quadrats_transect(domain, P$N_TRANSECTS, P$N_QUADRATS_PER_TRANSECT, P$QUADRAT_SIZE, P$TRANSECT_ANGLE),
    "voronoi"    = place_quadrats_voronoi(domain, P$N_QUADRATS, P$QUADRAT_SIZE, P$VORONOI_SEED_FACTOR),
    stop("Invalid 'SAMPLING_SCHEME'.")
  )
  if (nrow(quadrats) == 0) stop("Failed to place any quadrats.")

  all_species <- LETTERS[1:P$N_SPECIES]
  env_gradients <- create_environmental_gradients(
    domain,
    P$SAMPLING_RESOLUTION,
    P$ENVIRONMENTAL_NOISE,
    covariates = P$ENV_COVARIATES %||% NULL
  )
  abund_matrix <- create_abundance_matrix(species_dist, quadrats, all_species)
  site_env <- calculate_quadrat_environment(env_gradients, quadrats, sf::st_crs(domain))

  site_coords <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(quadrats))) |>
    as.data.frame() |>
    dplyr::mutate(site = quadrats$quadrat_id) |>
    dplyr::select(site, x = X, y = Y)

  if (tolower(as.character(P$DOMAIN_TYPE %||% "polygon")) %in% c("network", "coastline")) {
    bb <- sf::st_bbox(domain)
    if (tolower(as.character(P$LINEAR_AXIS %||% "x")) == "x") {
      den <- as.numeric(bb["xmax"] - bb["xmin"])
      if (!is.finite(den) || den == 0) den <- 1
      site_coords$linear_pos <- (site_coords$x - as.numeric(bb["xmin"])) / den
    } else {
      den <- as.numeric(bb["ymax"] - bb["ymin"])
      if (!is.finite(den) || den == 0) den <- 1
      site_coords$linear_pos <- (site_coords$y - as.numeric(bb["ymin"])) / den
    }
    site_coords$linear_pos <- pmax(0, pmin(1, site_coords$linear_pos))
    attr(site_coords, "distance_metric") <- "along_path"
    attr(site_coords, "linear_wrap") <- isTRUE(P$LINEAR_WRAP)
  }

  list(
    P = P,
    domain = domain,
    species_dist = species_dist,
    quadrats = quadrats,
    env_gradients = env_gradients,
    abund_matrix = abund_matrix,
    site_env = site_env,
    site_coords = site_coords
  )
}

# ---- Parameter materialisation (programmatic workflows) ----------------------
#
# Users often call load_config() and then tweak canonical fields like
# GRADIENT_* or QUADRAT_SIZE_OPTION. Those tweaks must be re-materialised into
# derived fields (P$GRADIENT and P$QUADRAT_SIZE) before the engine runs.
#
# This helper is internal and designed to be idempotent.
#
# @keywords internal
# @noRd
spesim_materialize_P <- function(P) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (is.null(P) || !is.list(P)) stop("`P` must be a list.")

  # Point-process strings
  P$SPATIAL_PROCESS_A <- tolower(as.character(P$SPATIAL_PROCESS_A %||% "poisson"))
  P$SPATIAL_PROCESS_OTHERS <- tolower(as.character(P$SPATIAL_PROCESS_OTHERS %||% "poisson"))
  P$MODEL_FAMILY <- tolower(as.character(P$MODEL_FAMILY %||% "manual"))
  P$MODEL_FAMILY <- gsub("-", "_", P$MODEL_FAMILY)
  P$DOMAIN_TYPE <- tolower(as.character(P$DOMAIN_TYPE %||% "polygon"))
  P$DOMAIN_TYPE <- gsub("-", "_", P$DOMAIN_TYPE)
  P$LINEAR_AXIS <- tolower(as.character(P$LINEAR_AXIS %||% "x"))
  P$DISTANCE_METRIC <- tolower(as.character(P$DISTANCE_METRIC %||% "auto"))
  P$DISTANCE_METRIC <- gsub("-", "_", P$DISTANCE_METRIC)
  P$DISPERSAL_KERNEL <- tolower(as.character(P$DISPERSAL_KERNEL %||% "gaussian"))
  P$DISPERSAL_KERNEL <- gsub("-", "_", P$DISPERSAL_KERNEL)
  P$NEUTRAL_META_MODEL <- tolower(as.character(P$NEUTRAL_META_MODEL %||% "zsm"))
  P$NEUTRAL_META_MODEL <- gsub("_", "-", P$NEUTRAL_META_MODEL)
  if (!P$SPATIAL_PROCESS_A %in% c("poisson", "thomas")) {
    stop("SPATIAL_PROCESS_A must be 'poisson' or 'thomas'.")
  }
  if (!P$SPATIAL_PROCESS_OTHERS %in% c("poisson", "strauss", "geyer", "thomas")) {
    stop("SPATIAL_PROCESS_OTHERS must be 'poisson', 'thomas', 'strauss', or 'geyer'.")
  }
  if (!P$MODEL_FAMILY %in% c("manual", "niche_filtering", "neutral_csr", "neutral_hubbell_like", "hybrid")) {
    stop("MODEL_FAMILY must be one of: manual, niche_filtering, neutral_csr, neutral_hubbell_like, hybrid.")
  }
  if (!P$DISPERSAL_KERNEL %in% c("gaussian", "exponential", "power_law")) {
    stop("DISPERSAL_KERNEL must be one of: gaussian, exponential, power_law.")
  }
  if (!P$DOMAIN_TYPE %in% c("polygon", "network", "coastline")) {
    stop("DOMAIN_TYPE must be one of: polygon, network, coastline.")
  }
  if (!P$LINEAR_AXIS %in% c("x", "y")) {
    stop("LINEAR_AXIS must be 'x' or 'y'.")
  }
  if (!P$DISTANCE_METRIC %in% c("auto", "euclidean", "along_path")) {
    stop("DISTANCE_METRIC must be one of: auto, euclidean, along_path.")
  }
  if (!is.finite(as.numeric(P$DISPERSAL_DIRECTION_BIAS %||% 0)) ||
      as.numeric(P$DISPERSAL_DIRECTION_BIAS %||% 0) < -1 ||
      as.numeric(P$DISPERSAL_DIRECTION_BIAS %||% 0) > 1) {
    stop("DISPERSAL_DIRECTION_BIAS must be in [-1, 1].")
  }
  ok_sad <- c(
    "fisher", "geometric", "brokenstick", "zipf", "zipf-mandelbrot",
    "lognormal", "poisson-lognormal", "poisson-gamma", "zsm", "custom"
  )
  if (P$NEUTRAL_META_MODEL %in% c("broken-stick", "brokenstick")) P$NEUTRAL_META_MODEL <- "brokenstick"
  if (P$NEUTRAL_META_MODEL %in% c("zipfmandelbrot", "zipf-mandelbrot")) P$NEUTRAL_META_MODEL <- "zipf-mandelbrot"
  if (P$NEUTRAL_META_MODEL %in% c("poissonlognormal", "poisson-lognormal")) P$NEUTRAL_META_MODEL <- "poisson-lognormal"
  if (P$NEUTRAL_META_MODEL %in% c("poissongamma", "poisson-gamma")) P$NEUTRAL_META_MODEL <- "poisson-gamma"
  if (!P$NEUTRAL_META_MODEL %in% ok_sad) {
    stop("NEUTRAL_META_MODEL must match a supported SAD model.")
  }

  # Gradients: rebuild P$GRADIENT from the canonical GRADIENT_* fields
  gs <- as.character(P$GRADIENT_SPECIES %||% character(0))
  ga <- as.character(P$GRADIENT_ASSIGNMENTS %||% character(0))
  if (length(gs) != length(ga)) {
    stop("GRADIENT_SPECIES and GRADIENT_ASSIGNMENTS must have the same length.")
  }
  if (length(gs) && !all(ga %in% c("temperature", "elevation", "rainfall"))) {
    stop("Unknown gradient(s): ", paste(setdiff(ga, c("temperature", "elevation", "rainfall")), collapse = ", "))
  }

  opt_raw <- .parse_named_pairs_numeric(P$GRADIENT_OPTIMA %||% 0.5)
  tol_raw <- .parse_named_pairs_numeric(P$GRADIENT_TOLERANCE %||% 0.1)

  .resolve_param_vector <- function(values, key_species, key_gradients) {
    uq_grad <- unique(key_gradients)
    valnames <- names(values)

    if (length(values) == 1L && is.numeric(values)) {
      return(rep(as.numeric(values), length(key_species)))
    }
    if (!is.null(valnames) && all(key_species %in% valnames)) {
      return(as.numeric(values[key_species]))
    }
    if (!is.null(valnames) && all(uq_grad %in% valnames)) {
      return(as.numeric(values[key_gradients]))
    }
    if (length(values) == length(key_species)) {
      return(as.numeric(values))
    }
    if (length(values) == length(uq_grad)) {
      map <- stats::setNames(as.numeric(values), uq_grad)
      return(as.numeric(map[key_gradients]))
    }
    stop(
      "Cannot resolve gradient parameters; supply scalar, length |species|, ",
      "named by species, named by gradients, or length |unique(gradients)|."
    )
  }

  if (length(gs)) {
    opt <- .resolve_param_vector(opt_raw, gs, ga)
    tol <- .resolve_param_vector(tol_raw, gs, ga)
    clamp01 <- function(x) pmax(0, pmin(1, x))
    opt <- clamp01(opt)
    if (any(!is.finite(tol) | tol <= 0)) {
      stop("GRADIENT_TOLERANCE must be positive and finite.")
    }
    P$GRADIENT <- tibble::tibble(
      species  = as.character(gs),
      gradient = as.character(ga),
      optimum  = as.numeric(opt),
      tol      = as.numeric(tol)
    )
  } else {
    P$GRADIENT <- NULL
  }

  # Quadrat size
  QUADRAT_SIZES <- list(
    small  = c(1, 1),
    medium = c(1.5, 1.5),
    large  = c(2, 2)
  )
  qs <- QUADRAT_SIZES[[as.character(P$QUADRAT_SIZE_OPTION %||% "medium")]]
  if (is.null(qs)) stop("Invalid QUADRAT_SIZE_OPTION. Use small|medium|large.")
  P$QUADRAT_SIZE <- qs

  P
}

# ---- Interaction resolution --------------------------------------------------
#
# @keywords internal
# @noRd
spesim_resolve_interactions <- function(P,
                                       interactions_file = NULL,
                                       interactions_validate = TRUE,
                                       interactions_strict = TRUE,
                                       interactions_print = TRUE,
                                       interactions_top_n = 20) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # 1) explicit function argument wins
  src_file <- interactions_file

  # 2) else, check for a pointer inside P
  if (is.null(src_file) && !is.null(P$INTERACTIONS_FILE)) {
    src_file <- as.character(P$INTERACTIONS_FILE)
  }

  # 3) load from file | inline rules | neutral fallback
  if (!is.null(src_file) && nzchar(src_file)) {
    I <- load_interactions(src_file, P$N_SPECIES)
  } else if (!is.null(P$INTERACTIONS_EDGELIST)) {
    I <- load_interactions_inline(
      rules = as.character(P$INTERACTIONS_EDGELIST),
      n_species = P$N_SPECIES,
      radius = as.numeric(P$INTERACTION_RADIUS %||% 0)
    )
  } else {
    I <- list(
      radius = as.numeric(P$INTERACTION_RADIUS %||% 0),
      matrix = matrix(1, P$N_SPECIES, P$N_SPECIES,
        dimnames = list(LETTERS[1:P$N_SPECIES], LETTERS[1:P$N_SPECIES])
      )
    )
  }

  if (isTRUE(interactions_validate)) {
    validate_interactions(
      I,
      spp_names = LETTERS[1:P$N_SPECIES],
      stop_on_error = isTRUE(interactions_strict)
    )
  }
  if (isTRUE(interactions_print)) {
    print_interactions(I, digits = 3, top_n = interactions_top_n)
  }

  P$INTERACTION_RADIUS <- I$radius
  P$INTERACTION_MATRIX <- I$matrix

  list(P = P, I = I)
}
