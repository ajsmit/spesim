# -----------------------------------------------------------------------------
# SPATIAL SAMPLING SIMULATION (Package-ready, refactored)
# -----------------------------------------------------------------------------
# Author: Gemini (refactor)
# Date: 2025-08-08
# -----------------------------------------------------------------------------

# ---- 1. SETUP & LIBRARIES ---------------------------------------------------

#' Install Packages If Missing
#'
#' @description Checks for, and installs, any missing packages.
#' @param pkgs Character vector of package names.
#' @return Invisibly returns \code{NULL}. Installs packages for their side effects.
#' @keywords internal
install_if_missing <- function(pkgs) {
  new_pkgs <- pkgs[!(pkgs %in% utils::installed.packages()[, "Package"])]
  if (length(new_pkgs)) {
    message("Installing missing packages: ", paste(new_pkgs, collapse = ", "))
    utils::install.packages(new_pkgs, dependencies = TRUE)
  }
  invisible(NULL)
}

.required_packages <- c(
  "sf", "ggplot2", "dplyr", "tidyr", "viridis", "stringr", "patchwork",
  "FNN", "RANN", "lwgeom", "vegan", "colorspace", "metR", "tibble", "scales"
)

install_if_missing(.required_packages)

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(vegan)
  library(viridis)
  library(stringr)
  library(patchwork)
  library(FNN)
  library(RANN)
  library(lwgeom)
  library(metR)
  library(colorspace)
  library(tibble)
  library(scales)
})

# Internal helpers
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

#' @keywords internal
.strip_quotes <- function(x) gsub('^\\s*["\\\']?|["\\\']?\\s*$', "", x)

#' @keywords internal
.split_vector <- function(val) {
  val <- trimws(val)
  if (grepl("^c\\s*\\(", val) && grepl("\\)\\s*$", val)) {
    inner <- sub("^c\\s*\\((.*)\\)\\s*$", "\\1", val)
    items <- strsplit(inner, ",")[[1]]
  } else if (grepl(",", val)) {
    items <- strsplit(val, ",")[[1]]
  } else {
    return(.strip_quotes(val))
  }
  trimws(.strip_quotes(items))
}

#' @keywords internal
.parse_named_pairs <- function(items) {
  if (length(items) == 0) {
    return(NULL)
  }
  m <- regexec("^([^:]+):(.+)$", items)
  parts <- regmatches(items, m)
  if (any(vapply(parts, length, integer(1)) != 3)) {
    return(NULL)
  }
  keys <- trimws(vapply(parts, `[`, character(1), 2))
  vals <- trimws(vapply(parts, `[`, character(1), 3))
  nums <- suppressWarnings(as.numeric(vals))
  out <- if (!anyNA(nums)) nums else vals
  names(out) <- keys
  out
}

#' @keywords internal
.parse_named_numeric <- function(x) {
  parts <- strsplit(as.character(x), "\\s*,\\s*")[[1]]
  out <- numeric(length(parts))
  nms <- character(length(parts))
  for (i in seq_along(parts)) {
    kv <- strsplit(parts[i], "\\s*:\\s*")[[1]]
    if (length(kv) == 2) {
      nms[i] <- kv[1]
      out[i] <- as.numeric(kv[2])
    } else {
      nms[i] <- NA_character_
      out[i] <- as.numeric(kv[1])
    }
  }
  if (any(!is.na(nms))) names(out) <- nms
  out
}

#' @keywords internal
.parse_matrix_spec <- function(x, S, spp_names) {
  row_strs <- if (length(x) == 1L) unlist(strsplit(x, "\\s*[;|]\\s*")) else x
  rows <- lapply(row_strs, function(r) as.numeric(strsplit(trimws(r), "\\s*,\\s*")[[1]]))
  if (length(rows) != S || any(vapply(rows, length, 1L) != S)) {
    stop("Row-delimited INTERACTION_MATRIX must be S rows of S numeric values (S×S).")
  }
  mat <- do.call(rbind, rows)
  dimnames(mat) <- list(spp_names, spp_names)
  mat
}

#' @keywords internal
.parse_named_pairs_numeric <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.character(x)) {
    return(x)
  }
  items <- unlist(strsplit(x, "\\s*,\\s*"))
  has_colon <- grepl(":", items, fixed = TRUE)
  if (!any(has_colon)) {
    return(x)
  }
  nms <- sub(":.*$", "", items)
  vals_chr <- sub("^[^:]*:", "", items)
  vals_num <- suppressWarnings(as.numeric(vals_chr))
  if (any(is.na(vals_num))) stop("Could not parse named pairs: ", paste(items[is.na(vals_num)], collapse = ", "))
  stats::setNames(vals_num, trimws(nms))
}

# ---- 2. CONFIGURATION -------------------------------------------------------

#' Parse a Parameter File (Simple Key = Value)
#'
#' @description Reads a configuration file with \code{key = value} pairs.
#'   Trailing comments beginning with \code{#} are allowed and ignored; hex colors like \code{#22223b} are preserved.
#'   Supports scalars, comma-separated vectors, \code{c(...)} vectors, and \code{name:value} entries.
#'
#' @param filename Path to the initialization file (text).
#'
#' @return A named list of parsed parameters.
#' @examples
#' \dontrun{
#' params <- parse_init_file("simul_init.txt")
#' }
#' @export
#' @importFrom utils read.csv
parse_init_file <- function(filename) {
  if (!file.exists(filename)) stop("Parameter file not found: ", filename)
  lines <- readLines(filename, warn = FALSE)
  # Only strip a trailing comment starting with whitespace + '#'
  lines <- sub("\\s+#.*$", "", lines)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  params <- list()
  for (ln in lines) {
    if (!grepl("=", ln, fixed = TRUE)) next
    kv <- strsplit(ln, "=", fixed = TRUE)[[1]]
    key <- toupper(trimws(kv[1]))
    val <- trimws(paste(kv[-1], collapse = "="))
    val <- gsub('^"|"$', "", val)
    # vector?
    if (grepl(",", val) || grepl("^c\\s*\\(", val)) {
      items <- .split_vector(val)
      named <- .parse_named_pairs(items)
      if (!is.null(named)) {
        params[[key]] <- named
        next
      }
      nums <- suppressWarnings(as.numeric(items))
      params[[key]] <- if (!anyNA(nums)) nums else items
      next
    }
    # scalar coercions
    if (val %in% c("TRUE", "FALSE")) {
      params[[key]] <- as.logical(val)
    } else if (grepl("^[+-]?[0-9]+$", val)) {
      params[[key]] <- as.integer(val)
    } else if (grepl("^[+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?$", val)) {
      params[[key]] <- as.numeric(val)
    } else if (grepl("^[^:]+:.+$", val)) {
      parts <- strsplit(val, ":", fixed = TRUE)[[1]]
      nm <- trimws(parts[1])
      vv <- trimws(parts[2])
      num <- suppressWarnings(as.numeric(vv))
      out <- if (!is.na(num)) c(num) else c(vv)
      names(out) <- nm
      params[[key]] <- out
    } else {
      params[[key]] <- val
    }
  }
  params
}

#' Load Simulation Configuration (With Defaults)
#'
#' @description Loads parameters from a file and merges with defaults.
#'   Validates inputs, materializes gradient specs into a tidy table,
#'   resolves quadrat sizes, and validates colors.
#'
#' @param init_file Path to the configuration file.
#'
#' @return A fully-populated parameter list used by the simulator.
#' @export
#' @importFrom tibble tibble
#' @importFrom grDevices col2rgb
#' @importFrom stats setNames
load_config <- function(init_file) {
  message("========== INITIALISING SPATIAL SAMPLING SIMULATION ==========")

  config <- list(
    SEED = 77,
    OUTPUT_PREFIX = "output",
    N_INDIVIDUALS = 2000,
    N_SPECIES = 10,
    DOMINANT_FRACTION = 0.30,
    FISHER_ALPHA = 3.0,
    FISHER_X = 0.95,
    GRADIENT_SPECIES = character(0),
    GRADIENT_ASSIGNMENTS = character(0),
    GRADIENT_OPTIMA = NULL,
    GRADIENT_TOLERANCE = NULL,
    SAMPLING_RESOLUTION = 50,
    ENVIRONMENTAL_NOISE = 0.05,
    MAX_CLUSTERS_DOMINANT = 5,
    CLUSTER_SPREAD_DOMINANT = 3.0,
    INTERACTION_RADIUS = 0,
    INTERACTION_MATRIX = NULL,
    SAMPLING_SCHEME = "random",
    N_QUADRATS = 20,
    QUADRAT_SIZE_OPTION = "medium",
    N_TRANSECTS = 1,
    N_QUADRATS_PER_TRANSECT = 8,
    TRANSECT_ANGLE = 90,
    VORONOI_SEED_FACTOR = 10,
    POINT_SIZE = 0.2,
    POINT_ALPHA = 1.0,
    QUADRAT_ALPHA = 0.05,
    BACKGROUND_COLOUR = "white",
    FOREGROUND_COLOUR = "#22223b",
    QUADRAT_COLOUR = "black",
    ADVANCED_ANALYSIS = FALSE
  )

  if (!file.exists(init_file)) stop("init_file not found: ", init_file)
  raw <- readLines(init_file, warn = FALSE)
  strip_comments <- function(x) sub("\\s+#.*$", "", x) # keeps hex colors
  lines <- trimws(vapply(raw, strip_comments, character(1)))
  lines <- lines[nzchar(lines)]

  kvs <- list()
  i <- 1L
  while (i <= length(lines)) {
    ln <- lines[i]
    if (!grepl("=", ln, fixed = TRUE)) {
      i <- i + 1
      next
    }
    parts <- strsplit(ln, "=", fixed = TRUE)[[1]]
    key <- trimws(parts[1])
    val <- trimws(paste(parts[-1], collapse = "="))
    if (grepl("^c\\s*\\(", val)) {
      open <- stringr::str_count(val, "\\(")
      close <- stringr::str_count(val, "\\)")
      while (open > close && i < length(lines)) {
        i <- i + 1L
        val <- paste0(val, " ", lines[i])
        open <- stringr::str_count(val, "\\(")
        close <- stringr::str_count(val, "\\)")
      }
    }
    kvs[[length(kvs) + 1L]] <- list(key = key, val = val)
    i <- i + 1L
  }

  parse_value <- function(v) {
    v <- trimws(v)
    v <- gsub('^"(.*)"$', "\\1", v)
    v <- gsub("^'(.*)'$", "\\1", v)
    if (grepl("^c\\s*\\(", v)) {
      out <- tryCatch(eval(parse(text = v)), error = function(e) e)
      if (inherits(out, "error")) stop("Could not parse vector literal: ", v, "\n", out$message)
      return(out)
    }
    if (grepl(",", v, fixed = TRUE)) {
      items <- trimws(strsplit(v, ",", fixed = TRUE)[[1]])
      nums <- suppressWarnings(as.numeric(items))
      return(if (!anyNA(nums)) nums else items)
    }
    if (v %in% c("TRUE", "FALSE")) {
      return(as.logical(v))
    }
    if (grepl("^[+-]?[0-9]+$", v)) {
      return(as.integer(v))
    }
    if (grepl("^[+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?$", v)) {
      return(as.numeric(v))
    }
    if (grepl("^[^:]+:.+$", v)) {
      kv <- strsplit(v, ":", fixed = TRUE)[[1]]
      nm <- trimws(kv[1])
      vv <- trimws(kv[2])
      num <- suppressWarnings(as.numeric(vv))
      out <- if (!is.na(num)) c(num) else c(vv)
      names(out) <- nm
      return(out)
    }
    v
  }

  for (kv in kvs) config[[kv$key]] <- parse_value(kv$val)

  .validate_colour <- function(x, default) {
    x <- as.character(x)
    ok <- tryCatch(
      {
        grDevices::col2rgb(x)
        TRUE
      },
      error = function(...) FALSE
    )
    if (ok) x else default
  }
  config$BACKGROUND_COLOUR <- .validate_colour(config$BACKGROUND_COLOUR, "white")
  config$FOREGROUND_COLOUR <- .validate_colour(config$FOREGROUND_COLOUR, "#22223b")
  config$QUADRAT_COLOUR <- .validate_colour(config$QUADRAT_COLOUR, "black")

  num_fields <- c(
    "SEED", "N_INDIVIDUALS", "N_SPECIES", "DOMINANT_FRACTION", "FISHER_ALPHA", "FISHER_X",
    "SAMPLING_RESOLUTION", "ENVIRONMENTAL_NOISE", "MAX_CLUSTERS_DOMINANT", "CLUSTER_SPREAD_DOMINANT",
    "INTERACTION_RADIUS", "N_QUADRATS", "N_TRANSECTS", "N_QUADRATS_PER_TRANSECT", "TRANSECT_ANGLE",
    "VORONOI_SEED_FACTOR", "POINT_SIZE", "POINT_ALPHA", "QUADRAT_ALPHA"
  )
  for (f in intersect(names(config), num_fields)) config[[f]] <- as.numeric(config[[f]])
  for (f in intersect(names(config), "ADVANCED_ANALYSIS")) config[[f]] <- as.logical(as.character(config[[f]]))

  gs <- as.character(config$GRADIENT_SPECIES)
  ga <- as.character(config$GRADIENT_ASSIGNMENTS)
  if (length(gs) != length(ga)) stop("GRADIENT_SPECIES and GRADIENT_ASSIGNMENTS must be same length.")
  allowed_grads <- c("temperature", "elevation", "rainfall")
  if (!all(ga %in% allowed_grads)) {
    stop(
      "Unknown gradient(s): ", paste(setdiff(ga, allowed_grads), collapse = ", "),
      ". Allowed: temperature, elevation, rainfall."
    )
  }

  opt_raw <- .parse_named_pairs_numeric(config$GRADIENT_OPTIMA %||% 0.5)
  tol_raw <- .parse_named_pairs_numeric(config$GRADIENT_TOLERANCE %||% 0.1)

  .resolve_param_vector <- function(values, key_species, key_gradients) {
    uq_grad <- unique(key_gradients)
    val_names <- names(values)
    if (length(values) == 1L && is.numeric(values)) {
      return(rep(as.numeric(values), length(key_species)))
    }
    if (!is.null(val_names) && all(key_species %in% val_names)) {
      return(as.numeric(values[key_species]))
    }
    if (!is.null(val_names) && all(uq_grad %in% val_names)) {
      return(as.numeric(values[key_gradients]))
    }
    if (length(values) == length(key_species)) {
      return(as.numeric(values))
    }
    if (length(values) == length(uq_grad)) {
      map <- stats::setNames(as.numeric(values), uq_grad)
      return(as.numeric(map[key_gradients]))
    }
    stop("Cannot resolve parameter vector for gradients; supply scalar, length S, named by species, named by gradients, or length |unique(gradients)|.")
  }

  opt <- .resolve_param_vector(opt_raw, gs, ga)
  tol <- .resolve_param_vector(tol_raw, gs, ga)
  clamp01 <- function(x) pmax(0, pmin(1, x))
  opt <- clamp01(opt)
  if (any(!is.finite(tol) | tol <= 0)) stop("GRADIENT_TOLERANCE must be positive and finite.")

  config$GRADIENT <- tibble::tibble(
    species  = gs,
    gradient = ga,
    optimum  = as.numeric(opt),
    tol      = as.numeric(tol)
  )
  config$GRADIENT$species <- as.character(config$GRADIENT$species)

  QUADRAT_SIZES <- list(small = c(1, 1), medium = c(1.5, 1.5), large = c(2, 2))
  qs <- QUADRAT_SIZES[[as.character(config$QUADRAT_SIZE_OPTION)]]
  if (is.null(qs)) stop("Invalid QUADRAT_SIZE_OPTION. Use small|medium|large.")
  config$QUADRAT_SIZE <- qs

  set.seed(config$SEED)
  config
}

#' Load Interspecific Interactions
#'
#' @description Loads a radius and interaction matrix from a separate config file.
#'   Supports a full \code{MATRIX_CSV} or an \code{EDGELIST_CSV} with columns
#'   \code{focal, neighbor, value}. If missing, returns a neutral (all-ones) matrix.
#'
#' @param interactions_file Path to an interactions config file (may be \code{NULL}).
#' @param n_species Integer; number of species (defines matrix size).
#'
#' @return A list with \code{radius} (numeric) and \code{matrix} (S×S numeric).
#' @export
#' @importFrom utils read.csv
load_interactions <- function(interactions_file, n_species) {
  spp_names <- LETTERS[1:n_species]
  out_radius <- 0
  IM <- matrix(1.0, nrow = n_species, ncol = n_species, dimnames = list(spp_names, spp_names))

  if (is.null(interactions_file) || !file.exists(interactions_file)) {
    if (!is.null(interactions_file) && !file.exists(interactions_file)) {
      warning("interactions_file not found: ", interactions_file, " — interactions disabled.")
    }
    return(list(radius = out_radius, matrix = IM))
  }

  raw <- readLines(interactions_file, warn = FALSE)
  raw <- trimws(sub("\\s+#.*$", "", raw))
  raw <- raw[nzchar(raw)]
  kv <- strsplit(raw, "=", fixed = TRUE)
  K <- toupper(trimws(vapply(kv, `[`, "", 1)))
  V <- trimws(vapply(kv, function(x) paste(x[-1], collapse = "="), ""))
  params <- as.list(stats::setNames(V, K))

  getp <- function(nm, default = NULL) {
    val <- params[[nm]]
    if (is.null(val) || !nzchar(val)) default else val
  }

  rstr <- getp("INTERACTION_RADIUS")
  if (!is.null(rstr)) {
    rnum <- suppressWarnings(as.numeric(rstr))
    if (length(rnum) == 1 && is.finite(rnum) && rnum >= 0) {
      out_radius <- rnum
    } else {
      warning("Bad INTERACTION_RADIUS; using 0 (disabled).")
      out_radius <- 0
    }
  }

  matrix_csv <- getp("MATRIX_CSV")
  edgelist_csv <- getp("EDGELIST_CSV")
  auto <- tolower(getp("AUTO", "false")) %in% c("true", "1", "yes")

  if (!is.null(matrix_csv)) {
    if (!file.exists(matrix_csv)) stop("MATRIX_CSV not found: ", matrix_csv)
    M <- as.matrix(utils::read.csv(matrix_csv, row.names = 1, check.names = FALSE))
    fullM <- matrix(1.0, nrow = n_species, ncol = n_species, dimnames = list(spp_names, spp_names))
    rn <- intersect(rownames(M), spp_names)
    cn <- intersect(colnames(M), spp_names)
    if (length(rn) && length(cn)) {
      suppressWarnings({
        fullM[rn, cn] <- as.numeric(M[rn, cn, drop = FALSE])
      })
    }
    IM <- fullM
  } else if (!is.null(edgelist_csv)) {
    if (!file.exists(edgelist_csv)) stop("EDGELIST_CSV not found: ", edgelist_csv)
    E <- utils::read.csv(edgelist_csv, stringsAsFactors = FALSE, check.names = FALSE)
    need_cols <- c("focal", "neighbor", "value")
    if (!all(need_cols %in% names(E))) stop("EDGELIST_CSV must have columns: focal, neighbor, value")
    for (i in seq_len(nrow(E))) {
      f <- as.character(E$focal[i])
      n <- as.character(E$neighbor[i])
      v <- suppressWarnings(as.numeric(E$value[i]))
      if (f %in% spp_names && n %in% spp_names && is.finite(v)) IM[f, n] <- v
    }
  } else if (auto) {
    if (all(c("D", "E") %in% spp_names)) IM["E", "D"] <- 3.0
  }

  if (any(!is.finite(IM))) stop("INTERACTION_MATRIX contains non-finite values.")
  list(radius = out_radius, matrix = IM)
}

# ---- 3. CORE SIMULATION -----------------------------------------------------

#' Create an Irregular Sampling Domain
#'
#' @description Generates a single organic polygon to define the study area.
#' @return An \code{sf} polygon (\code{GEOMETRYCOLLECTION} of one polygon).
#' @export
#' @importFrom sf st_polygon st_sfc st_sf
create_sampling_domain <- function() {
  theta <- seq(0, 2 * pi, length.out = 20)
  r <- 10 + 3 * sin(3 * theta) + 2 * cos(5 * theta) + runif(20, -1, 1)
  x <- r * cos(theta) + runif(20, -0.5, 0.5)
  y <- r * sin(theta) * 0.8 + runif(20, -0.5, 0.5)
  coords <- cbind(x, y)
  coords <- rbind(coords, coords[1, ])
  polygon <- sf::st_polygon(list(coords))
  sf::st_sf(geometry = sf::st_sfc(polygon))
}

#' Generate Fisher's Log-Series Abundances
#'
#' @param n_species Total number of species (integer).
#' @param n_individuals Total number of individuals (integer).
#' @param dominant_fraction Fraction of individuals assigned to species "A" (numeric in [0,1]).
#' @param alpha Fisher's alpha parameter (numeric).
#' @param x Log-series scaling parameter near 1 (numeric).
#'
#' @return Named numeric vector of abundances for species A, B, C, ...
#' @export
generate_fisher_log_series <- function(n_species, n_individuals, dominant_fraction, alpha, x) {
  n_dominant <- round(n_individuals * dominant_fraction)
  n_remaining <- n_individuals - n_dominant
  ranks <- 2:n_species
  rel_abundances <- alpha * (x^ranks) / ranks
  abundances <- round(rel_abundances / sum(rel_abundances) * n_remaining)
  all_abundances <- c(n_dominant, abundances)
  adjustment <- n_individuals - sum(all_abundances)
  if (length(all_abundances) >= 2) {
    all_abundances[2] <- all_abundances[2] + adjustment
  } else {
    all_abundances[1] <- all_abundances[1] + adjustment
  }
  species_names <- LETTERS[1:n_species]
  names(all_abundances) <- species_names
  all_abundances[all_abundances > 0]
}

#' Create Environmental Gradients
#'
#' @param domain \code{sf} polygon defining the domain.
#' @param resolution Grid resolution (number of steps per axis).
#' @param noise_level Standard deviation of Gaussian noise added to gradients.
#'
#' @return Data frame with columns \code{x, y, temperature, elevation, rainfall, temperature_C, elevation_m, rainfall_mm}.
#' @export
#' @importFrom sf st_bbox
create_environmental_gradients <- function(domain, resolution, noise_level) {
  bbox <- sf::st_bbox(domain)
  x_seq <- seq(bbox["xmin"], bbox["xmax"], length.out = resolution)
  y_seq <- seq(bbox["ymin"], bbox["ymax"], length.out = resolution)
  grid <- expand.grid(x = x_seq, y = y_seq)

  x_norm <- (grid$x - bbox["xmin"]) / (bbox["xmax"] - bbox["xmin"])
  y_norm <- (grid$y - bbox["ymin"]) / (bbox["ymax"] - bbox["ymin"])

  temp_vals <- x_norm * 0.7 + y_norm * 0.3 + rnorm(nrow(grid), 0, noise_level)
  dist_center <- sqrt((x_norm - 0.5)^2 + (y_norm - 0.5)^2)
  elev_vals <- 1 - (dist_center / sqrt(0.5^2 + 0.5^2)) + rnorm(nrow(grid), 0, noise_level)
  rain_vals <- (-x_norm * 0.6 + y_norm * 0.8)
  rain_vals <- (rain_vals - min(rain_vals)) / (max(rain_vals) - min(rain_vals)) + rnorm(nrow(grid), 0, noise_level)

  grid$temperature <- pmax(0, pmin(1, temp_vals))
  grid$elevation <- pmax(0, pmin(1, elev_vals))
  grid$rainfall <- pmax(0, pmin(1, rain_vals))

  grid$temperature_C <- grid$temperature * 30 - 2
  grid$elevation_m <- grid$elevation * 2000
  grid$rainfall_mm <- grid$rainfall * 700 + 200
  grid
}

#' Generate Species Distribution with Interactions
#'
#' @description Assigns species identities to random points in the domain
#'   given abundances, environmental preferences, dominant clustering for A,
#'   and a local interspecific interaction modifier within a radius.
#'
#' @param domain \code{sf} polygon domain.
#' @param P List of parameters from \code{\link{load_config}} (requires fields:
#'   \code{N_SPECIES, N_INDIVIDUALS, DOMINANT_FRACTION, FISHER_ALPHA, FISHER_X,
#'   SAMPLING_RESOLUTION, ENVIRONMENTAL_NOISE, GRADIENT, MAX_CLUSTERS_DOMINANT,
#'   CLUSTER_SPREAD_DOMINANT, INTERACTION_RADIUS, INTERACTION_MATRIX}).
#'
#' @return \code{sf} point object with columns for environmental values and \code{species}.
#' @export
#' @importFrom sf st_as_sf st_sample st_join st_nearest_feature st_coordinates
#' @importFrom FNN get.knnx
#' @importFrom RANN nn2
generate_heterogeneous_distribution <- function(domain, P) {
  abundances <- generate_fisher_log_series(P$N_SPECIES, P$N_INDIVIDUALS, P$DOMINANT_FRACTION, P$FISHER_ALPHA, P$FISHER_X)
  env_grid <- create_environmental_gradients(domain, P$SAMPLING_RESOLUTION, P$ENVIRONMENTAL_NOISE)
  env_sf <- sf::st_as_sf(env_grid, coords = c("x", "y"), crs = sf::st_crs(domain))
  all_points <- sf::st_sample(domain, size = P$N_INDIVIDUALS, type = "random")
  points_with_env <- sf::st_join(sf::st_sf(geometry = all_points), env_sf, join = sf::st_nearest_feature)
  points_with_env$species <- ""
  available_indices <- seq_len(nrow(points_with_env))

  for (sp in names(abundances)) {
    n_ind <- abundances[sp]
    if (length(available_indices) < n_ind) n_ind <- length(available_indices)
    if (n_ind == 0) next

    base_probs <- rep(1.0, length(available_indices))
    if (!is.null(P$GRADIENT) && sp %in% P$GRADIENT$species) {
      row <- P$GRADIENT[P$GRADIENT$species == sp, ][1, ]
      env_values <- points_with_env[[row$gradient]][available_indices]
      base_probs <- exp(-((env_values - row$optimum)^2) / (2 * row$tol^2))
    } else if (sp == "A") {
      n_avail <- length(available_indices)
      n_clusters <- min(max(1, P$MAX_CLUSTERS_DOMINANT %||% 5), n_avail)
      spread <- P$CLUSTER_SPREAD_DOMINANT %||% 3
      if (n_clusters >= 1 && n_avail >= 1) {
        centre_idx_rel <- sample(seq_len(n_avail), n_clusters)
        centre_idx_abs <- available_indices[centre_idx_rel]
        centre_coords <- sf::st_coordinates(points_with_env[centre_idx_abs, , drop = FALSE])
        avail_coords <- sf::st_coordinates(points_with_env[available_indices, , drop = FALSE])
        nn <- RANN::nn2(data = centre_coords, query = avail_coords, k = 1)
        d <- as.numeric(nn$nn.dists[, 1])
        base_probs <- exp(-d / spread)
        if (!all(is.finite(base_probs)) || all(base_probs <= 0) || anyNA(base_probs)) base_probs <- rep(1.0, n_avail)
      } else {
        base_probs <- rep(1.0, n_avail)
      }
    }

    interaction_modifier <- rep(1.0, length(available_indices))
    assigned_indices <- which(points_with_env$species != "")
    if (length(assigned_indices) > 0 && is.finite(P$INTERACTION_RADIUS) && P$INTERACTION_RADIUS > 0) {
      if (!is.null(dim(P$INTERACTION_MATRIX)) && !is.na(match(sp, rownames(P$INTERACTION_MATRIX)))) {
        available_coords <- as.matrix(sf::st_coordinates(points_with_env[available_indices, , drop = FALSE]))
        assigned_coords <- as.matrix(sf::st_coordinates(points_with_env[assigned_indices, , drop = FALSE]))
        if (NROW(assigned_coords) > 0 && NROW(available_coords) > 0) {
          k_use <- min(5L, NROW(assigned_coords))
          knn <- FNN::get.knnx(data = assigned_coords, query = available_coords, k = k_use)
          nn_index <- knn$nn.index
          nn_dist <- knn$nn.dist
          interaction_scores <- vapply(seq_len(NROW(nn_index)), function(row_i) {
            idxs <- nn_index[row_i, ]
            dists <- nn_dist[row_i, ]
            keep <- which(is.finite(dists) & dists <= P$INTERACTION_RADIUS & idxs > 0)
            if (length(keep) == 0L) {
              return(1.0)
            }
            neighbour_assigned_idx <- assigned_indices[idxs[keep]]
            neighbour_species <- points_with_env$species[neighbour_assigned_idx]
            iv <- P$INTERACTION_MATRIX[sp, neighbour_species, drop = TRUE]
            iv[is.na(iv)] <- 1.0
            iv <- pmax(iv, 1e-12)
            exp(mean(log(iv)))
          }, numeric(1))
          interaction_modifier <- interaction_scores
        }
      }
    }

    final_probs <- base_probs * interaction_modifier
    if (all(!is.finite(final_probs)) || all(final_probs <= 0)) final_probs <- rep(1, length(available_indices))
    selected <- sample(available_indices, size = n_ind, prob = final_probs, replace = FALSE)
    if (length(selected) > 0) {
      points_with_env$species[selected] <- sp
      available_indices <- setdiff(available_indices, selected)
    }
  }
  dplyr::filter(points_with_env, species != "")
}

#' Create a Quadrat Polygon From Center
#'
#' @param center_point \code{sf} point geometry (single point).
#' @param size Numeric vector \code{c(width, height)}.
#'
#' @return \code{sf} polygon geometry (single rectangle).
#' @keywords internal
#' @importFrom sf st_coordinates st_polygon st_sfc
create_quadrat_from_center <- function(center_point, size) {
  half_w <- size[1] / 2
  half_h <- size[2] / 2
  coords <- sf::st_coordinates(center_point)[1, ]
  sf::st_polygon(list(cbind(
    c(coords["X"] - half_w, coords["X"] + half_w, coords["X"] + half_w, coords["X"] - half_w, coords["X"] - half_w),
    c(coords["Y"] - half_h, coords["Y"] - half_h, coords["Y"] + half_h, coords["Y"] + half_h, coords["Y"] - half_h)
  )))
}

#' Randomly Place Non-Overlapping Quadrats
#'
#' @param domain \code{sf} polygon domain.
#' @param n_quadrats Target number of quadrats.
#' @param quadrat_size Numeric vector \code{c(width, height)}.
#'
#' @return \code{sf} polygon layer of quadrats with \code{quadrat_id}.
#' @export
#' @importFrom sf st_bbox st_sfc st_within st_intersects st_sf st_crs
place_quadrats <- function(domain, n_quadrats, quadrat_size) {
  bbox <- sf::st_bbox(domain)
  quadrats <- list()
  attempts <- 0
  max_attempts <- n_quadrats * 100
  while (length(quadrats) < n_quadrats && attempts < max_attempts) {
    attempts <- attempts + 1
    x_center <- runif(1, bbox["xmin"] + quadrat_size[1] / 2, bbox["xmax"] - quadrat_size[1] / 2)
    y_center <- runif(1, bbox["ymin"] + quadrat_size[2] / 2, bbox["ymax"] - quadrat_size[2] / 2)
    center_pt_sfc <- sf::st_sfc(sf::st_point(c(x_center, y_center)), crs = sf::st_crs(domain))
    new_quadrat_poly <- create_quadrat_from_center(center_pt_sfc, quadrat_size) |> sf::st_sfc(crs = sf::st_crs(domain))
    within_mat <- sf::st_within(new_quadrat_poly, domain, sparse = FALSE)
    is_within <- isTRUE(within_mat[, 1, drop = TRUE][1])
    if (is_within) {
      is_overlapping <- FALSE
      if (length(quadrats) > 0) {
        existing_quadrats_sfc <- do.call(c, quadrats)
        overlap_mat <- sf::st_intersects(new_quadrat_poly, existing_quadrats_sfc, sparse = FALSE)
        is_overlapping <- any(overlap_mat[1, , drop = TRUE])
      }
      if (!is_overlapping) quadrats[[length(quadrats) + 1]] <- new_quadrat_poly
    }
  }
  if (length(quadrats) == 0) {
    return(sf::st_sf(quadrat_id = integer(0), geometry = sf::st_sfc(crs = sf::st_crs(domain))))
  }
  final_quadrats <- do.call(c, quadrats) |> sf::st_sf(quadrat_id = 1:length(.), geometry = .)
  final_quadrats
}

#' Tiled (Systematic Cell) Quadrat Placement
#'
#' @param domain \code{sf} polygon domain.
#' @param n_quadrats Target number of quadrats.
#' @param quadrat_size Numeric vector \code{c(width, height)}.
#'
#' @return \code{sf} polygons of selected quadrats.
#' @export
#' @importFrom sf st_make_grid st_within st_sf st_crs
place_quadrats_tiled <- function(domain, n_quadrats, quadrat_size) {
  candidate_grid <- sf::st_make_grid(domain, cellsize = quadrat_size, what = "polygons")
  within_mat <- sf::st_within(candidate_grid, domain, sparse = FALSE)
  inside <- if (is.matrix(within_mat)) drop(within_mat[, 1, drop = TRUE]) else as.logical(within_mat)
  valid_locations <- candidate_grid[inside]
  num_possible <- length(valid_locations)
  if (num_possible == 0) {
    warning("Systematic placement failed: no quadrats fit inside the domain.")
    return(sf::st_sf(quadrat_id = integer(0), geometry = sf::st_sfc(crs = sf::st_crs(domain))))
  }
  if (num_possible < n_quadrats) {
    warning(sprintf("Could only place %d of %d requested quadrats.", num_possible, n_quadrats))
    n_quadrats <- num_possible
  }
  sampled_indices <- sample.int(num_possible, size = n_quadrats)
  final_quadrats_sfc <- valid_locations[sampled_indices]
  sf::st_sf(quadrat_id = seq_len(n_quadrats), geometry = final_quadrats_sfc)
}

#' Voronoi-based Quadrat Placement
#'
#' @param domain \code{sf} polygon domain.
#' @param n_quadrats Target number of quadrats.
#' @param quadrat_size Numeric vector \code{c(width, height)}.
#' @param voronoi_seed_factor Multiplier controlling number of initial seeds.
#'
#' @return \code{sf} polygons of quadrats.
#' @export
#' @importFrom sf st_union st_sample st_voronoi st_cast st_intersection st_inscribed_circle st_centroid st_sfc st_sf st_area
place_quadrats_voronoi <- function(domain, n_quadrats, quadrat_size, voronoi_seed_factor) {
  n_seeds <- n_quadrats * voronoi_seed_factor
  domain_union <- sf::st_union(domain)
  seed_points <- sf::st_sample(domain_union, size = n_seeds, type = "random")
  voronoi_polys <- sf::st_voronoi(sf::st_union(seed_points))
  voronoi_clipped <- sf::st_intersection(sf::st_cast(voronoi_polys), domain_union)
  inscribed_circles <- suppressWarnings(sf::st_inscribed_circle(voronoi_clipped))
  quadrat_half_diag <- sqrt(quadrat_size[1]^2 + quadrat_size[2]^2) / 2
  radii <- sqrt(sf::st_area(inscribed_circles) / pi)
  suitable_indices <- which(radii >= quadrat_half_diag)
  num_possible <- length(suitable_indices)
  if (num_possible == 0) {
    warning("Voronoi placement failed: no suitable cells.")
    return(sf::st_sf(quadrat_id = integer(0), geometry = sf::st_sfc(crs = sf::st_crs(domain))))
  }
  if (num_possible < n_quadrats) {
    warning(sprintf("Voronoi placement found %d suitable locations; using all.", num_possible))
    n_quadrats <- num_possible
  }
  sampled_indices <- sample(suitable_indices, size = n_quadrats)
  final_centers <- sf::st_centroid(inscribed_circles[sampled_indices])
  quadrat_list <- lapply(sf::st_geometry(final_centers), function(pt) create_quadrat_from_center(pt, quadrat_size))
  final_quadrats_sfc <- sf::st_sfc(quadrat_list, crs = sf::st_crs(domain))
  sf::st_sf(quadrat_id = 1:length(final_quadrats_sfc), geometry = final_quadrats_sfc)
}

#' Systematic Grid Quadrat Placement
#'
#' @param domain \code{sf} polygon domain.
#' @param n_quadrats Approximate number of quadrats.
#' @param quadrat_size Numeric vector \code{c(width, height)}.
#'
#' @return \code{sf} polygons of quadrats.
#' @export
#' @importFrom sf st_bbox st_make_grid st_sfc st_within st_sf st_crs
place_quadrats_systematic <- function(domain, n_quadrats, quadrat_size) {
  bbox <- sf::st_bbox(domain)
  aspect_ratio <- (bbox["ymax"] - bbox["ymin"]) / (bbox["xmax"] - bbox["xmin"])
  nx <- round(sqrt(n_quadrats / aspect_ratio))
  ny <- round(aspect_ratio * nx)
  candidate_centers <- sf::st_make_grid(domain, n = c(nx, ny), what = "centers")
  candidate_quadrats <- sf::st_sfc(lapply(candidate_centers, function(pt) create_quadrat_from_center(pt, quadrat_size)), crs = sf::st_crs(domain))
  within_mat <- sf::st_within(candidate_quadrats, domain, sparse = FALSE)
  is_within <- within_mat[, 1, drop = TRUE]
  valid_quadrats <- candidate_quadrats[is_within]
  if (length(valid_quadrats) == 0) {
    warning("Systematic sampling failed to place any quadrats.")
    return(sf::st_sf(quadrat_id = integer(0), geometry = sf::st_sfc(crs = sf::st_crs(domain))))
  }
  sf::st_sf(quadrat_id = seq_len(length(valid_quadrats)), geometry = valid_quadrats)
}

#' Parallel-Transect Quadrat Placement
#'
#' @param domain \code{sf} polygon domain.
#' @param n_transects Number of parallel transects.
#' @param n_quadrats_per_transect Quadrats per transect.
#' @param quadrat_size Numeric vector \code{c(width, height)}.
#' @param angle Compass angle in degrees (0 = North, 90 = East).
#'
#' @return \code{sf} polygons of quadrats.
#' @export
#' @importFrom sf st_buffer st_is_empty st_area st_bbox st_centroid st_linestring st_sfc st_crs st_intersection st_coordinates st_sf
place_quadrats_transect <- function(domain, n_transects, n_quadrats_per_transect, quadrat_size, angle) {
  buffer_dist <- sqrt(quadrat_size[1]^2 + quadrat_size[2]^2) / 2
  safe_domain <- sf::st_buffer(domain, -buffer_dist)
  if (sf::st_is_empty(safe_domain) || sf::st_area(safe_domain) == 0) {
    warning("Domain too small for given quadrat size.")
    return(sf::st_sf(quadrat_id = integer(0), geometry = sf::st_sfc(crs = sf::st_crs(domain))))
  }

  bbox <- sf::st_bbox(domain)
  center <- sf::st_centroid(sf::st_as_sfc(bbox))
  diag_len <- sqrt((bbox["xmax"] - bbox["xmin"])^2 + (bbox["ymax"] - bbox["ymin"])^2) * 1.5
  y_coords <- seq(bbox["ymin"], bbox["ymax"], length.out = n_transects + 2)[2:(n_transects + 1)]

  math_angle_deg <- 90 - angle
  math_angle_rad <- math_angle_deg * pi / 180
  rot_matrix <- matrix(c(cos(math_angle_rad), sin(math_angle_rad), -sin(math_angle_rad), cos(math_angle_rad)), 2, 2)

  horizontal_lines <- sf::st_sfc(lapply(y_coords, function(y) {
    sf::st_linestring(matrix(c(
      sf::st_coordinates(center)[1] - diag_len / 2,
      sf::st_coordinates(center)[1] + diag_len / 2,
      y, y
    ), ncol = 2))
  }), crs = sf::st_crs(domain))
  rotated_lines <- (horizontal_lines - center) * rot_matrix + center
  sf::st_crs(rotated_lines) <- sf::st_crs(domain)

  points_list <- vector("list", n_transects)
  for (i in seq_len(n_transects)) {
    line <- rotated_lines[i]
    segment <- sf::st_intersection(line, safe_domain)
    if (sf::st_is_empty(segment) || sf::st_length(segment) == 0) next
    segment_coords <- sf::st_coordinates(segment)
    start_pt <- segment_coords[1, c("X", "Y")]
    end_pt <- segment_coords[nrow(segment_coords), c("X", "Y")]
    fractions <- if (n_quadrats_per_transect > 1) seq(0, 1, length.out = n_quadrats_per_transect) else 0.5
    x_coords <- start_pt[1] + fractions * (end_pt[1] - start_pt[1])
    y_coords <- start_pt[2] + fractions * (end_pt[2] - start_pt[2])
    points_list[[i]] <- sf::st_as_sf(data.frame(x = x_coords, y = y_coords), coords = c("x", "y"), crs = sf::st_crs(domain))
  }
  points_list <- Filter(Negate(is.null), points_list)
  if (length(points_list) == 0) {
    warning("No transects intersected the safe sampling area.")
    return(sf::st_sf(quadrat_id = integer(0), geometry = sf::st_sfc(crs = sf::st_crs(domain))))
  }
  candidate_centers <- do.call(rbind, points_list)
  quadrat_geometries <- sf::st_sfc(lapply(sf::st_geometry(candidate_centers), function(pt) create_quadrat_from_center(pt, quadrat_size)), crs = sf::st_crs(domain))
  sf::st_sf(quadrat_id = 1:length(quadrat_geometries), geometry = quadrat_geometries)
}

# ---- 4. ADVANCED ANALYSES ---------------------------------------------------

#' Rank-Abundance Data
#'
#' @param species_dist \code{sf} points with column \code{species}.
#' @param P Parameter list from \code{\link{load_config}}.
#'
#' @return Data frame with columns \code{Rank, Abundance, Source}.
#' @export
calculate_rank_abundance <- function(species_dist, P) {
  observed_counts <- table(species_dist$species)
  observed_data <- as.data.frame(observed_counts, stringsAsFactors = FALSE)
  names(observed_data) <- c("Species", "Abundance")
  observed_data <- observed_data %>%
    arrange(desc(Abundance)) %>%
    mutate(Rank = dplyr::row_number(), Source = "Observed") %>%
    select(Rank, Abundance, Source)

  theoretical_abundances <- generate_fisher_log_series(P$N_SPECIES, P$N_INDIVIDUALS, P$DOMINANT_FRACTION, P$FISHER_ALPHA, P$FISHER_X)
  theoretical_data <- tibble::tibble(
    Abundance = sort(as.numeric(theoretical_abundances), decreasing = TRUE),
    Rank = seq_along(theoretical_abundances),
    Source = "Theoretical"
  )
  dplyr::bind_rows(observed_data, theoretical_data)
}

#' Plot Rank-Abundance Curve
#'
#' @param rank_abundance_data Data frame from \code{\link{calculate_rank_abundance}}.
#' @return A \code{ggplot} object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_y_log10 scale_color_manual scale_linetype_manual scale_shape_manual labs theme_bw theme element_text
plot_rank_abundance <- function(rank_abundance_data) {
  ggplot(rank_abundance_data, aes(x = Rank, y = Abundance, color = Source)) +
    geom_line(aes(linetype = Source), linewidth = 1.1) +
    geom_point(aes(shape = Source), size = 3, fill = "white", stroke = 1.2) +
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_color_manual(values = c("Observed" = "black", "Theoretical" = "#e41a1c")) +
    scale_linetype_manual(values = c("Observed" = "solid", "Theoretical" = "dashed")) +
    scale_shape_manual(values = c("Observed" = 21, "Theoretical" = 22)) +
    labs(
      title = "Species-Abundance Distribution (SAD)",
      subtitle = "Observed vs. Fisher's log-series",
      x = "Species Rank",
      y = "Abundance (log10)",
      color = "Distribution",
      linetype = "Distribution",
      shape = "Distribution"
    ) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom", legend.title = element_text(face = "bold"), plot.title = element_text(face = "bold"))
}

#' Occupancy-Abundance Table
#'
#' @param abund_matrix Data frame: first column \code{site}, then species columns of non-negative integers.
#' @return Data frame with \code{Species, TotalAbundance, Occupancy}.
#' @export
calculate_occupancy_abundance <- function(abund_matrix) {
  abund_numeric <- abund_matrix[, -which(names(abund_matrix) == "site"), drop = FALSE]
  data.frame(
    Species = names(abund_numeric),
    TotalAbundance = colSums(abund_numeric),
    Occupancy = colSums(abund_numeric > 0),
    row.names = NULL
  )
}

#' Plot Occupancy-Abundance Relationship
#'
#' @param oa_data Data frame from \code{\link{calculate_occupancy_abundance}}.
#' @return A \code{ggplot} object.
#' @export
plot_occupancy_abundance <- function(oa_data) {
  if (nrow(oa_data) == 0 || all(oa_data$TotalAbundance == 0)) {
    return(ggplot2::ggplot() +
      ggplot2::labs(title = "Occupancy-Abundance Relationship", subtitle = "No data to plot.") +
      ggplot2::theme_void())
  }
  ggplot(oa_data, aes(x = TotalAbundance, y = Occupancy)) +
    geom_point(alpha = 0.6, size = 3) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.8) +
    scale_x_log10(labels = scales::label_log()) +
    scale_y_log10(labels = scales::label_log()) +
    labs(
      title = "Occupancy-Abundance Relationship",
      x = "Total Abundance (log10)",
      y = "Sites Occupied (log10)"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
}

#' Species-Area (Accumulation) Data
#'
#' @param abund_matrix Data frame: first column \code{site}, then species columns.
#' @return Data frame with \code{Sites, Richness, SD}.
#' @export
#' @importFrom vegan specaccum
calculate_species_area <- function(abund_matrix) {
  abund_numeric <- abund_matrix[, -which(names(abund_matrix) == "site"), drop = FALSE]
  sar_curve <- vegan::specaccum(abund_numeric, method = "random", permutations = 100)
  data.frame(Sites = sar_curve$sites, Richness = sar_curve$richness, SD = sar_curve$sd)
}

#' Plot Species-Area Relationship
#'
#' @param sar_data Data frame from \code{\link{calculate_species_area}}.
#' @return A \code{ggplot} object.
#' @export
plot_species_area <- function(sar_data) {
  ggplot(sar_data, aes(x = Sites, y = Richness)) +
    geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.4) +
    geom_line(linewidth = 1.2) +
    labs(
      title = "Species-Area Relationship (SAR)",
      x = "Number of Quadrats",
      y = "Cumulative Species Richness"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
}

#' Distance-Decay Data
#'
#' @param abund_matrix Data frame: first column \code{site}, then species columns.
#' @param site_coords Data frame with columns \code{site, x, y} (same order as \code{abund_matrix} rows).
#'
#' @return Data frame with \code{Distance, Dissimilarity} (Sørensen as binary Bray-Curtis).
#' @export
#' @importFrom stats dist
#' @importFrom vegan vegdist
calculate_distance_decay <- function(abund_matrix, site_coords) {
  coords <- site_coords[, c("x", "y")]
  abund_numeric <- abund_matrix[, -which(names(abund_matrix) == "site"), drop = FALSE]
  geo_dist <- stats::dist(coords, method = "euclidean")
  comm_dissim <- vegan::vegdist(abund_numeric, method = "bray", binary = TRUE)
  data.frame(Distance = as.vector(geo_dist), Dissimilarity = as.vector(comm_dissim))
}

#' Plot Distance-Decay Relationship
#'
#' @param decay_data Data frame from \code{\link{calculate_distance_decay}}.
#' @return A \code{ggplot} object.
#' @export
plot_distance_decay <- function(decay_data) {
  ggplot(decay_data, aes(x = Distance, y = Dissimilarity)) +
    geom_point(alpha = 0.3, shape = 16) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.1) +
    ylim(0, 1) +
    labs(
      title = "Distance-Decay of Community Similarity",
      x = "Geographic Distance",
      y = "Community Dissimilarity (Sørensen)"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
}

#' Rarefaction Curves Data
#'
#' @param abund_matrix Data frame: first column \code{site}, then species columns.
#' @return Data frame with \code{SiteID, SampleSize, RarefiedRichness}.
#' @export
#' @importFrom vegan rarecurve
calculate_rarefaction <- function(abund_matrix) {
  abund_numeric <- abund_matrix[, -which(names(abund_matrix) == "site"), drop = FALSE]
  site_ids <- abund_matrix$site
  rarefaction_list <- vegan::rarecurve(abund_numeric, step = 1)
  output_list <- vector("list", length(rarefaction_list))
  for (i in seq_along(rarefaction_list)) {
    richness_values <- rarefaction_list[[i]]
    sample_sizes <- attr(richness_values, "Subsample")
    output_list[[i]] <- data.frame(
      SiteID = as.factor(site_ids[i]),
      SampleSize = sample_sizes,
      RarefiedRichness = richness_values
    )
  }
  do.call(rbind, output_list)
}

#' Plot Rarefaction Curves
#'
#' @param rarefaction_data Data frame from \code{\link{calculate_rarefaction}}.
#' @return A \code{ggplot} object.
#' @export
plot_rarefaction <- function(rarefaction_data) {
  ggplot(rarefaction_data, aes(x = SampleSize, y = RarefiedRichness, group = SiteID, color = SiteID)) +
    geom_line(linewidth = 0.8) +
    scale_color_viridis_d(option = "plasma") +
    labs(
      title = "Rarefaction Curves by Quadrat",
      x = "Individuals Sampled",
      y = "Expected Species (Rarefied)",
      color = "Quadrat ID"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
}

#' Advanced Analysis Panel (Patchwork)
#'
#' @param res List produced by the main simulator (\code{P, species_dist, abund_matrix, site_coords}).
#' @return A \code{patchwork} object combining multiple plots.
#' @export
#' @importFrom patchwork plot_spacer plot_annotation
generate_advanced_panel <- function(res) {
  theme_panel <- theme(text = element_text(size = 11), legend.title = element_text(size = 10), legend.text = element_text(size = 9))
  p_rank <- plot_rank_abundance(calculate_rank_abundance(res$species_dist, res$P)) + labs(subtitle = NULL) + theme_panel + theme(legend.position = "bottom")
  p_oa <- plot_occupancy_abundance(calculate_occupancy_abundance(res$abund_matrix)) + labs(subtitle = NULL) + theme_panel
  p_sar <- plot_species_area(calculate_species_area(res$abund_matrix)) + labs(subtitle = NULL) + theme_panel
  p_decay <- plot_distance_decay(calculate_distance_decay(res$abund_matrix, res$site_coords)) + labs(subtitle = NULL) + theme_panel
  p_rare <- plot_rarefaction(calculate_rarefaction(res$abund_matrix)) + labs(subtitle = NULL) + theme_panel + guides(color = "none")
  (p_rank | p_oa) /
    (p_sar | p_decay) /
    (p_rare | patchwork::plot_spacer()) +
    plot_annotation(
      title = "Advanced Ecological Analysis Panel",
      theme = theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
    )
}

# ---- 5. DATA EXTRACTION -----------------------------------------------------

#' Site-by-Species Abundance Matrix
#'
#' @param species_dist \code{sf} points with \code{species}.
#' @param quadrats \code{sf} polygons with \code{quadrat_id}.
#' @param all_species_names Character vector of species column names to include (e.g., \code{LETTERS[1:S]}).
#'
#' @return Data frame with column \code{site} and species abundance columns.
#' @export
#' @importFrom sf st_intersection st_drop_geometry
#' @importFrom tidyr pivot_wider
create_abundance_matrix <- function(species_dist, quadrats, all_species_names) {
  intersections <- sf::st_intersection(species_dist, quadrats)
  if (nrow(intersections) == 0) {
    abund_df <- data.frame(site = quadrats$quadrat_id)
    for (sp in all_species_names) abund_df[[sp]] <- 0
    return(abund_df)
  }
  abund_df <- intersections %>%
    sf::st_drop_geometry() %>%
    dplyr::count(quadrat_id, species, name = "abundance") %>%
    tidyr::pivot_wider(names_from = species, values_from = abundance, values_fill = 0)

  missing_species <- setdiff(all_species_names, names(abund_df))
  for (sp in missing_species) abund_df[[sp]] <- 0

  data.frame(quadrat_id = quadrats$quadrat_id) %>%
    dplyr::left_join(abund_df, by = "quadrat_id") %>%
    dplyr::mutate(dplyr::across(-quadrat_id, ~ tidyr::replace_na(., 0))) %>%
    dplyr::select(site = quadrat_id, dplyr::all_of(all_species_names)) %>%
    dplyr::arrange(site)
}

#' Mean Environmental Conditions per Quadrat
#'
#' @param env_grid Data frame with columns \code{x, y, temperature_C, elevation_m, rainfall_mm}.
#' @param quadrats \code{sf} polygons with \code{quadrat_id}.
#' @param domain_crs CRS object or integer EPSG matching the domain/quadrats.
#'
#' @return Data frame with columns \code{site, temperature_C, elevation_m, rainfall_mm}.
#' @export
#' @importFrom sf st_as_sf st_join st_drop_geometry
calculate_quadrat_environment <- function(env_grid, quadrats, domain_crs) {
  env_sf <- sf::st_as_sf(env_grid, coords = c("x", "y"), crs = domain_crs)
  joined_data <- sf::st_join(quadrats, env_sf)
  site_env <- joined_data %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(site = quadrat_id) %>%
    dplyr::summarise(
      dplyr::across(c(temperature_C, elevation_m, rainfall_mm), ~ mean(., na.rm = TRUE)),
      .groups = "drop"
    )
  dplyr::left_join(data.frame(site = quadrats$quadrat_id), site_env, by = "site")
}

# ---- 6. PLOTTING & REPORT ---------------------------------------------------

#' Master Spatial Plot
#'
#' @description Plot the domain, points, and quadrats, optionally overlaying an environmental gradient.
#'
#' @param domain \code{sf} polygon.
#' @param species \code{sf} points with column \code{species}.
#' @param quadrats \code{sf} polygons with \code{quadrat_id}.
#' @param P Parameter list (colors, sizes, etc.).
#' @param show_gradient Logical; overlay a gradient raster/contour.
#' @param env_gradients Data frame from \code{\link{create_environmental_gradients}}.
#' @param gradient_type Character; one of \code{"temperature_C"}, \code{"elevation_m"}, \code{"rainfall_mm"}.
#'
#' @return A \code{ggplot} object.
#' @export
#' @importFrom ggplot2 ggplot geom_sf geom_sf_text aes scale_color_manual labs theme_void theme element_rect element_text coord_sf
plot_spatial_sampling <- function(domain, species, quadrats, P,
                                  show_gradient = FALSE, env_gradients = NULL,
                                  gradient_type = "temperature_C") {
  all_species_names <- LETTERS[1:P$N_SPECIES]
  species_colors <- rev(colorspace::sequential_hcl(P$N_SPECIES, palette = "RdPu"))
  names(species_colors) <- all_species_names

  p <- ggplot() +
    geom_sf(data = domain, fill = "grey70", color = P$FOREGROUND_COLOUR, linewidth = 0.5, alpha = 0.4)

  if (isTRUE(show_gradient) && !is.null(env_gradients) && gradient_type %in% names(env_gradients)) {
    p <- p +
      metR::geom_contour_fill(data = env_gradients, aes(x = x, y = y, z = .data[[gradient_type]]), alpha = 0.5) +
      scale_fill_viridis(option = "viridis", name = stringr::str_to_title(gsub("_", " ", gradient_type)), guide = "colorbar")
  }

  p +
    geom_sf(data = species, aes(color = species), size = P$POINT_SIZE, alpha = P$POINT_ALPHA) +
    scale_color_manual(values = species_colors, name = "Species", drop = FALSE) +
    geom_sf(data = quadrats, fill = NA, color = P$QUADRAT_COLOUR, linewidth = 0.4) +
    geom_sf_text(data = quadrats, aes(label = quadrat_id), color = P$QUADRAT_COLOUR, size = 2.5, fontface = "bold") +
    coord_sf(expand = FALSE) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = P$BACKGROUND_COLOUR, color = NA),
      plot.title = element_text(color = P$FOREGROUND_COLOUR, size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(color = P$FOREGROUND_COLOUR, size = 10, hjust = 0.5, margin = margin(b = 10)),
      legend.position = "none"
    ) +
    labs(
      title = if (isTRUE(show_gradient)) stringr::str_replace_all(stringr::str_to_title(gradient_type), "_", " ") else "Species Distribution",
      subtitle = paste(P$N_INDIVIDUALS, "Individuals |", P$N_SPECIES, "Species |", nrow(quadrats), "Quadrats")
    )
}

#' Generate Full Text Report
#'
#' @param res List of results (\code{P, env_gradients, species_dist, quadrats, abund_matrix, site_coords}).
#' @return Character scalar containing a formatted multi-section report.
#' @export
#' @importFrom stats cor cor.test dist
#' @importFrom vegan diversity vegdist fisher.alpha
generate_full_report <- function(res) {
  .opt_to_units <- function(opt, gname) {
    switch(gname,
      "temperature" = opt * 30 - 2,
      "elevation"   = opt * 2000,
      "rainfall"    = opt * 700 + 200,
      NA_real_
    )
  }
  .tol_to_units <- function(tol, gname) {
    switch(gname,
      "temperature" = tol * 30,
      "elevation"   = tol * 2000,
      "rainfall"    = tol * 700,
      NA_real_
    )
  }
  .fmt_grad_line <- function(rows, gname, units_label) {
    if (nrow(rows) == 0) {
      return("None")
    }
    paste(sprintf(
      "%s (opt=%.2f [%.1f %s], tol=%.2f [%.1f %s])",
      rows$species, rows$optimum, .opt_to_units(rows$optimum, gname), units_label,
      rows$tol, .tol_to_units(rows$tol, gname), units_label
    ), collapse = ", ")
  }

  env_report <- c("\nEnvironmental Gradients:")
  temp_range <- range(res$env_gradients$temperature_C, na.rm = TRUE)
  elev_range <- range(res$env_gradients$elevation_m, na.rm = TRUE)
  rain_range <- range(res$env_gradients$rainfall_mm, na.rm = TRUE)

  if (!is.null(res$P$GRADIENT)) {
    G <- res$P$GRADIENT
    temp_species_str <- .fmt_grad_line(G[G$gradient == "temperature", , drop = FALSE], "temperature", "°C")
    elev_species_str <- .fmt_grad_line(G[G$gradient == "elevation", , drop = FALSE], "elevation", "m")
    rain_species_str <- .fmt_grad_line(G[G$gradient == "rainfall", , drop = FALSE], "rainfall", "mm")
  } else {
    temp_species_str <- "None"
    elev_species_str <- "None"
    rain_species_str <- "None"
  }

  env_report <- c(
    env_report,
    sprintf("  Temperature: %.1f–%.1f °C (range: %.1f °C)", temp_range[1], temp_range[2], diff(temp_range)),
    "    Pattern: Diagonal (NW cool → SE warm)",
    paste0("    Responsive species: ", temp_species_str),
    sprintf("  Elevation: %.0f–%.0f m (range: %.0f m)", elev_range[1], elev_range[2], diff(elev_range)),
    "    Pattern: Central peak (mountain-like topology)",
    paste0("    Responsive species: ", elev_species_str),
    sprintf("  Rainfall: %.0f–%.0f mm (range: %.0f mm)", rain_range[1], rain_range[2], diff(rain_range)),
    "    Pattern: Perpendicular (NE dry → SW wet)",
    paste0("    Responsive species: ", rain_species_str)
  )

  cor_mat <- cor(res$env_gradients[, c("temperature_C", "elevation_m", "rainfall_mm")], use = "complete.obs")
  interp <- if (max(abs(cor_mat[upper.tri(cor_mat)])) < 0.3) "Gradients are approximately orthogonal (low correlation)" else "Some gradient correlation detected"
  corr_report <- c(
    "\nGradient Correlations:",
    sprintf("  Temperature-Elevation: r=%.3f", cor_mat[1, 2]),
    sprintf("  Temperature-Rainfall:  r=%.3f", cor_mat[1, 3]),
    sprintf("  Elevation-Rainfall:    r=%.3f", cor_mat[2, 3]),
    paste0("  Interpretation: ", interp)
  )

  sad <- as.data.frame(table(res$species_dist$species)) %>%
    `colnames<-`(c("Species", "Count")) %>%
    arrange(desc(Count)) %>%
    mutate(
      Percent = 100 * Count / sum(Count),
      Role = dplyr::case_when(
        !is.null(res$P$GRADIENT) & Species %in% res$P$GRADIENT$species ~ {
          g <- res$P$GRADIENT$gradient[match(Species, res$P$GRADIENT$species)]
          sprintf("[%s-RESPONSIVE]", toupper(g))
        },
        Species == "A" ~ "[DOMINANT - clustered]",
        TRUE ~ "[SUBORDINATE]"
      )
    )
  sad_report <- c(
    "\nSpecies Abundance Distribution:",
    sprintf("  Species %s: %3d individuals (%5.1f%%) %s", sad$Species, sad$Count, sad$Percent, sad$Role)
  )

  alpha_report <- c("\nSpatial Distribution of Alpha Diversity:")
  for (i in 1:nrow(res$quadrats)) {
    spp_in_q <- suppressWarnings(sf::st_intersection(res$species_dist, res$quadrats[i, ]))
    if (nrow(spp_in_q) > 0) {
      abunds <- as.data.frame(table(spp_in_q$species))
      abunds_str <- paste(sprintf("%s(%s)", abunds$Var1, abunds$Freq), collapse = ", ")
      alpha_report <- c(
        alpha_report,
        sprintf("  Quadrat %2d: α = %2d species | N = %3d individuals", i, dplyr::n_distinct(spp_in_q$species), nrow(spp_in_q)),
        sprintf("    Species: %s", abunds_str)
      )
    } else {
      alpha_report <- c(
        alpha_report,
        sprintf("  Quadrat %2d: α = %2d species | N = %3d individuals", i, 0, 0),
        "    Species: None"
      )
    }
  }

  abund_data <- res$abund_matrix %>% dplyr::select(-site)
  richness_data <- abund_data %>% dplyr::mutate(richness = rowSums(. > 0), n_ind = rowSums(.))
  shannon_H <- vegan::diversity(abund_data, index = "shannon")
  simpson_D <- vegan::diversity(abund_data, index = "simpson")
  mean_alpha <- mean(richness_data$richness)
  se_alpha <- stats::sd(richness_data$richness) / sqrt(nrow(richness_data))
  gamma_div <- dplyr::n_distinct(res$species_dist$species)
  beta_whittaker <- gamma_div / mean_alpha
  beta_additive <- gamma_div - mean_alpha
  pa_matrix <- abund_data %>%
    as.matrix() %>%
    `>`(0) %>%
    `*`(1)
  mean_sorensen <- mean(vegan::vegdist(pa_matrix, method = "bray", binary = TRUE))

  div_report <- c(
    "\nDiversity Partitioning:",
    sprintf("Alpha (mean local richness): %.2f ± %.2f SE", mean_alpha, se_alpha),
    sprintf("Shannon's H' (mean): %.3f ± %.3f SE", mean(shannon_H), stats::sd(shannon_H) / sqrt(length(shannon_H))),
    sprintf("Simpson's (1-D, mean): %.3f ± %.3f SE", mean(simpson_D), stats::sd(simpson_D) / sqrt(length(simpson_D))),
    sprintf("Gamma (regional species pool): %d species", gamma_div),
    sprintf("Beta (Whittaker): %.2f", beta_whittaker),
    sprintf("Beta (additive): %.2f", beta_additive),
    sprintf("Mean pairwise beta (Sørensen): %.3f", mean_sorensen),
    sprintf("Mean quadrat abundance: %.1f ± %.1f", mean(richness_data$n_ind), stats::sd(richness_data$n_ind)),
    sprintf("Abundance variation (CV): %.3f", stats::sd(richness_data$n_ind) / mean(richness_data$n_ind))
  )

  space_dist <- stats::dist(res$site_coords[, c("x", "y")])
  richness_dist <- stats::dist(richness_data$richness)
  mantel_interp <- "  No significant spatial autocorrelation."
  if (nrow(res$site_coords) >= 4 && var(richness_data$richness) > 0) {
    mantel_test <- suppressWarnings(stats::cor.test(space_dist, richness_dist, method = "pearson"))
    mantel_r <- mantel_test$estimate
    mantel_p <- mantel_test$p.value
    if (!is.na(mantel_p) && mantel_p < 0.05) {
      mantel_interp <- if (mantel_r > 0) "  Significant positive autocorrelation." else "  Significant negative autocorrelation."
    }
  } else {
    mantel_r <- NA
    mantel_p <- NA
  }
  spat_report <- c(
    "\nSpatial Autocorrelation Analysis:",
    sprintf("Spatial autocorrelation in richness (Mantel r): %.3f", mantel_r),
    sprintf("Statistical significance: p = %.3f", mantel_p),
    mantel_interp,
    sprintf("  Mean inter-quadrat distance: %.2f", mean(space_dist)),
    sprintf("  Richness variance: %.2f", var(richness_data$richness))
  )

  obs_abund <- sort(table(res$species_dist$species), decreasing = TRUE)
  n_rem <- res$P$N_INDIVIDUALS * (1 - res$P$DOMINANT_FRACTION)
  ranks <- 2:res$P$N_SPECIES
  rel_abund_theory <- res$P$FISHER_ALPHA * (res$P$FISHER_X^ranks) / ranks
  theory_abund <- c(res$P$N_INDIVIDUALS * res$P$DOMINANT_FRACTION, rel_abund_theory / sum(rel_abund_theory) * n_rem)
  len_obs <- length(obs_abund)
  len_theory <- length(theory_abund)
  max_len <- max(len_obs, len_theory)
  length(obs_abund) <- max_len
  length(theory_abund) <- max_len
  obs_abund[is.na(obs_abund)] <- 0
  residuals <- obs_abund - theory_abund
  rmse <- sqrt(mean(residuals^2, na.rm = TRUE))
  r_squared <- if (var(obs_abund, na.rm = TRUE) > 0) stats::cor(obs_abund, theory_abund, use = "complete.obs")^2 else 1
  eff_alpha <- vegan::fisher.alpha(table(res$species_dist$species))

  fisher_report <- c(
    "\nFisher's Log Series Model Validation:",
    sprintf("  RMSE: %.2f", rmse),
    sprintf("  R-squared: %.3f", r_squared),
    sprintf("  Max residual: %.1f", max(abs(residuals), na.rm = TRUE)),
    sprintf("  Specified alpha: %.2f", res$P$FISHER_ALPHA),
    sprintf("  Effective alpha from data: %.2f", eff_alpha)
  )

  full_report <- c(
    "========== ANALYSIS REPORT ==========",
    env_report, corr_report, sad_report, alpha_report, div_report, spat_report, fisher_report,
    "\nSIMULATION COMPLETED SUCCESSFULLY."
  )
  paste(full_report, collapse = "\n")
}

# ---- 7. MAIN ORCHESTRATOR ---------------------------------------------------

#' Run the Full Spatial Sampling Simulation
#'
#' @description Top-level function that orchestrates the simulation:
#'   loads config and interactions, simulates communities, places quadrats,
#'   writes tables and figures, optionally runs advanced analyses, and
#'   appends a text report.
#'
#' @param init_file Path to the simulation config file (text).
#' @param interactions_file Optional path to interactions config (text). If \code{NULL}, interactions default to neutral.
#' @param output_prefix String used as the base for all output filenames (timestamp appended automatically).
#'
#' @return Invisibly returns \code{NULL}. Files are written to disk.
#' @export
#' @importFrom sf st_crs st_centroid st_coordinates
#' @importFrom patchwork plot_layout
#' @importFrom ggplot2 ggsave theme
run_spatial_simulation <- function(init_file = "simul_init.txt",
                                   interactions_file = NULL,
                                   output_prefix = "simulation_output") {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_prefix <- paste0(output_prefix, "_", timestamp)
  output_dir <- dirname(output_prefix)
  report_path <- paste0(output_prefix, "_report.txt")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  con <- NULL
  on.exit(
    {
      num_sinks <- sink.number()
      if (num_sinks > 0) for (i in seq_len(num_sinks)) sink()
      try(close(con), silent = TRUE)
    },
    add = TRUE
  )

  con <- file(report_path, open = "wt")
  sink(con, type = "output")
  sink(con, type = "message")

  results_list <- NULL
  tryCatch(
    {
      P <- load_config(init_file)

      I <- load_interactions(interactions_file, P$N_SPECIES)
      P$INTERACTION_RADIUS <- I$radius
      P$INTERACTION_MATRIX <- I$matrix

      cat("========== SIMULATION PARAMETERS ==========\n")
      cat(paste(capture.output(dput(P[!names(P) %in% "QUADRAT_SIZE"])), collapse = "\n"), "\n\n")

      cat("========== RUNNING SIMULATION ==========\n")
      domain <- create_sampling_domain()
      species_dist <- generate_heterogeneous_distribution(domain, P)

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
      env_gradients <- create_environmental_gradients(domain, P$SAMPLING_RESOLUTION, P$ENVIRONMENTAL_NOISE)
      abund_matrix <- create_abundance_matrix(species_dist, quadrats, all_species)
      site_env <- calculate_quadrat_environment(env_gradients, quadrats, sf::st_crs(domain))
      site_coords <- suppressWarnings(sf::st_coordinates(sf::st_centroid(quadrats))) |>
        as.data.frame() |>
        dplyr::mutate(site = quadrats$quadrat_id) |>
        dplyr::select(site, x = X, y = Y)

      cat("\n========== SAVING OUTPUT FILES ==========\n")
      f_abund <- paste0(output_prefix, "_abundances.csv")
      write.csv(abund_matrix, f_abund, row.names = FALSE)
      f_env <- paste0(output_prefix, "_environments.csv")
      write.csv(site_env, f_env, row.names = FALSE)
      f_coord <- paste0(output_prefix, "_quadrat_centroids.csv")
      write.csv(site_coords, f_coord, row.names = FALSE)

      f_png <- paste0(output_prefix, "_fig_panel.png")
      panel_plot <- (plot_spatial_sampling(domain, species_dist, quadrats, P) |
        plot_spatial_sampling(domain, species_dist, quadrats, P, TRUE, env_gradients, "temperature_C")) /
        (plot_spatial_sampling(domain, species_dist, quadrats, P, TRUE, env_gradients, "elevation_m") |
          plot_spatial_sampling(domain, species_dist, quadrats, P, TRUE, env_gradients, "rainfall_mm")) +
        patchwork::plot_layout(guides = "collect") & theme(legend.position = "right")
      ggsave(f_png, panel_plot, width = 14, height = 10, dpi = 300, bg = "white")
      cat(sprintf("Main panel figure saved as %s\n", f_png))
      cat("Data tables saved.\n")

      results_list <- list(
        P = P, domain = domain, species_dist = species_dist, quadrats = quadrats,
        env_gradients = env_gradients, abund_matrix = abund_matrix, site_coords = site_coords
      )
    },
    error = function(e) {
      cat("\nERROR in core simulation:\n", e$message, "\n", sep = "")
    }
  )

  if (!is.null(results_list) && isTRUE(results_list$P$ADVANCED_ANALYSIS)) {
    tryCatch(
      {
        cat("\n========== ADVANCED ANALYSIS ==========\n")
        f_adv_panel <- paste0(output_prefix, "_fig_advanced_panel.png")
        advanced_panel_plot <- generate_advanced_panel(results_list)
        ggsave(f_adv_panel, advanced_panel_plot, width = 12, height = 14, dpi = 300, bg = "white")
        cat(sprintf("Advanced analysis panel saved as %s\n", f_adv_panel))
      },
      error = function(e) {
        cat("\nERROR during advanced analysis:\n", e$message, "\n", sep = "")
      }
    )
  }

  if (!is.null(results_list)) {
    cat("\n", generate_full_report(results_list), "\n", sep = "")
    cat("Outputs saved to: ", normalizePath(output_dir), "\n", sep = "")
  }
  invisible(NULL)
}
