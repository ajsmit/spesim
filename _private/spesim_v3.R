# -----------------------------------------------------------------------------
#
# SPATIAL SAMPLING SIMULATION
#
# Description
#   End-to-end simulator of a spatially heterogeneous ecological community in an
#   irregular domain. Individuals are placed according to:
#     • Fisher’s log-series species abundance distribution (with a dominant spp.).
#     • Environmental responses (per-species Gaussian along named gradients).
#     • Conspecific clustering for the dominant species.
#     • Local interspecific interactions (competition/facilitation) via a
#       radius-limited neighbourhood modifier.
#
#   The script then samples the landscape with user-selectable quadrat schemes,
#   assembles site×species tables, computes α/β/γ diversity and other metrics,
#   and writes a full text report plus publication-quality figures.
#
# Author: Gemini (Refactored from original script)
# Date: 2025-08-07
#
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# 1. SETUP & LIBRARIES
# -----------------------------------------------------------------------------

# Helper to install packages if missing
#' @title Install Packages
#' @description Checks for, and installs, any missing packages.
#' @param pkgs A character vector of package names.
install_if_missing <- function(pkgs) {
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new_pkgs)) {
    message(paste("Installing missing packages:", paste(new_pkgs, collapse = ", ")))
    install.packages(new_pkgs, dependencies = TRUE)
  }
}

# Load required libraries
required_packages <- c(
  "sf", "ggplot2", "dplyr", "tidyr", "viridis", "stringr", "patchwork", "FNN",
  "lwgeom", "vegan", "colorspace", "metR", "tibble", "scales", "RANN"
)

install_if_missing(required_packages)

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


# -----------------------------------------------------------------------------
# 2. CONFIGURATION & PARAMETER FUNCTIONS
# -----------------------------------------------------------------------------

# --- helpers for parse_init_file() -------------------------------------------
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
}

.strip_quotes <- function(x) gsub('^\\s*["\\\']?|["\\\']?\\s*$', "", x)

# Split a value that may be either "a,b,c" or c(a,b,c) or c("a","b","c")
# -> returns a trimmed character vector of items with quotes removed
.split_vector <- function(val) {
  val <- trimws(val)
  if (grepl("^c\\s*\\(", val) && grepl("\\)\\s*$", val)) {
    inner <- sub("^c\\s*\\((.*)\\)\\s*$", "\\1", val)
    # split on commas not caring about quotes (we'll strip them anyway)
    items <- strsplit(inner, ",")[[1]]
  } else if (grepl(",", val)) {
    items <- strsplit(val, ",")[[1]]
  } else {
    return(.strip_quotes(val))
  }
  trimws(.strip_quotes(items))
}

# If every item is "name:value", return a named vector (numeric when possible)
# Otherwise return NULL (caller will handle other coercions)
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

  if (!anyNA(nums)) {
    out <- nums
  } else {
    out <- vals
  }
  names(out) <- keys
  out
}

#' @title Parse Simulation Parameters
#' @description Reads a configuration file, parsing key-value pairs.
#'   It ignores comments, trims whitespace, and attempts to coerce values to
#'   appropriate types (numeric, integer, logical, vectors).
#' @param filename Path to the initialization file.
#' @return A named list of parameters.
parse_init_file <- function(filename) {
  if (!file.exists(filename)) stop(paste("Parameter file not found:", filename))

  lines <- readLines(filename, warn = FALSE)
  lines <- gsub("#.*", "", lines) # strip comments
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)] # drop empty

  params <- list()
  for (ln in lines) {
    if (!grepl("=", ln)) next
    kv <- strsplit(ln, "=", fixed = TRUE)[[1]]
    key <- toupper(trimws(kv[1]))
    val <- trimws(kv[2])

    # remove surrounding quotes around the whole value, if present
    val <- gsub('^"|"$', "", val)

    # 1) Vector forms ----------------------------------------------------------
    if (grepl(",", val) || grepl("^c\\s*\\(", val)) {
      items <- .split_vector(val)

      # a) key:value vectors -> named numeric/character vector
      named <- .parse_named_pairs(items)
      if (!is.null(named)) {
        params[[key]] <- named
        next
      }

      # b) plain numeric vector
      nums <- suppressWarnings(as.numeric(items))
      if (!anyNA(nums)) {
        params[[key]] <- nums
      } else {
        params[[key]] <- items
      }
      next
    }

    # 2) Scalars ---------------------------------------------------------------
    if (val %in% c("TRUE", "FALSE")) {
      params[[key]] <- as.logical(val)
    } else if (grepl("^[+-]?[0-9]+$", val)) {
      params[[key]] <- as.integer(val)
    } else if (grepl("^[+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?$", val)) {
      params[[key]] <- as.numeric(val)
    } else if (grepl("^[^:]+:.+$", val)) {
      # single "name:value" scalar -> named length-1 (numeric if possible)
      parts <- sub("^([^:]+):(.+)$", "\\1,\\2", val)
      pvec <- strsplit(parts, ",", fixed = TRUE)[[1]]
      nm <- trimws(pvec[1])
      vv <- trimws(pvec[2])
      num <- suppressWarnings(as.numeric(vv))
      if (!is.na(num)) {
        out <- c(num)
        names(out) <- nm
      } else {
        out <- c(vv)
        names(out) <- nm
      }
      params[[key]] <- out
    } else {
      params[[key]] <- val
    }
  }
  params
}


# ---- helpers for gradient params -------------------------------------------
.resolve_param_vector <- function(values, key_species, key_gradients) {
  # values: can be length 1; length S; length G; or named by species or gradient.
  # key_species: character vector (e.g., c("D","E",...))
  # key_gradients: character vector aligned with key_species (e.g., c("temperature","elevation",...))
  uq_grad <- unique(key_gradients)

  # Coerce to character names if present
  val_names <- names(values)

  pick_by_species <- function() {
    # out <- rep(NA_real_, length(key_species))
    out <- as.numeric(values[key_species])
    out
  }
  pick_by_gradient <- function() {
    # out <- rep(NA_real_, length(key_species))
    # map named gradient values onto species by their assignment
    out <- as.numeric(values[key_gradients])
    out
  }

  if (length(values) == 1L && is.numeric(values)) {
    rep(as.numeric(values), length(key_species))
  } else if (!is.null(val_names) && all(key_species %in% val_names)) {
    pick_by_species()
  } else if (!is.null(val_names) && all(uq_grad %in% val_names)) {
    pick_by_gradient()
  } else if (length(values) == length(key_species)) {
    as.numeric(values)
  } else if (length(values) == length(uq_grad)) {
    # assume ordering corresponds to unique(key_gradients) sequence
    map <- setNames(as.numeric(values), uq_grad)
    as.numeric(map[key_gradients])
  } else {
    stop(
      "Cannot resolve parameter vector. Supply either:\n",
      "  • a single number (recycled),\n",
      "  • a vector length S (same order as GRADIENT_SPECIES),\n",
      "  • a vector named by species (A,B,C,...),\n",
      "  • a vector named by gradients (temperature,elevation,rainfall), or\n",
      "  • a vector length G (unique gradients, order matching unique(GRADIENT_ASSIGNMENTS))."
    )
  }
}

.clamp01 <- function(x) pmax(0, pmin(1, x))


# helper: parse "key:value" form into named numeric vector
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


# helper: parse row-delimited string or vector into S×S numeric matrix
.parse_matrix_spec <- function(x, S, spp_names) {
  if (length(x) == 1L) {
    row_strs <- unlist(strsplit(x, "\\s*[;|]\\s*"))
  } else {
    row_strs <- x
  }
  rows <- lapply(row_strs, function(r) as.numeric(strsplit(trimws(r), "\\s*,\\s*")[[1]]))
  if (length(rows) != S || any(vapply(rows, length, 1L) != S)) {
    stop("Row-delimited INTERACTION_MATRIX must be S rows of S numeric values (S×S).")
  }
  mat <- do.call(rbind, rows)
  dimnames(mat) <- list(spp_names, spp_names)
  mat
}


.parse_named_pairs_numeric <- function(x) {
  # Accepts a character scalar or vector like c("temperature:0.5","elevation:0.3")
  # Returns a named numeric vector c(temperature=0.5, elevation=0.3)
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
  } # nothing to do
  nms <- sub(":.*$", "", items)
  vals_chr <- sub("^[^:]*:", "", items)
  vals_num <- suppressWarnings(as.numeric(vals_chr))
  if (any(is.na(vals_num))) stop("Could not parse named pairs: ", paste(items[is.na(vals_num)], collapse = ", "))
  stats::setNames(vals_num, trimws(nms))
}


#' @title Load Configuration with Defaults
#' @description Loads parameters from a file and merges them with a list of
#'    default values. User-defined values override defaults. It now dynamically
#'    generates the interaction matrix based on the final N_SPECIES.
#' @param init_file Path to the user's configuration file.
#' @return A list containing all configuration parameters for the simulation.
load_config <- function(init_file) {
  message("========== INITIALISING SPATIAL SAMPLING SIMULATION ==========")

  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

  # -------------------- defaults (complete) --------------------
  config <- list(
    SEED = 77,
    OUTPUT_PREFIX = "output",
    # community
    N_INDIVIDUALS = 2000,
    N_SPECIES = 10,
    DOMINANT_FRACTION = 0.30,
    FISHER_ALPHA = 3.0,
    FISHER_X = 0.95,
    # gradients
    GRADIENT_SPECIES = character(0),
    GRADIENT_ASSIGNMENTS = character(0),
    GRADIENT_OPTIMA = NULL,
    GRADIENT_TOLERANCE = NULL,
    SAMPLING_RESOLUTION = 50,
    ENVIRONMENTAL_NOISE = 0.05,
    # dominant clustering
    MAX_CLUSTERS_DOMINANT = 5,
    CLUSTER_SPREAD_DOMINANT = 3.0,
    # interactions
    INTERACTION_RADIUS = 0, # 0 disables interactions
    INTERACTION_MATRIX = NULL,
    # sampling/quadrats
    SAMPLING_SCHEME = "random", # random|systematic|transect|tiled|voronoi
    N_QUADRATS = 20,
    QUADRAT_SIZE_OPTION = "medium", # small|medium|large
    N_TRANSECTS = 1,
    N_QUADRATS_PER_TRANSECT = 8,
    TRANSECT_ANGLE = 90,
    VORONOI_SEED_FACTOR = 10,
    # plotting
    POINT_SIZE = 0.2,
    POINT_ALPHA = 1.0,
    QUADRAT_ALPHA = 0.05,
    BACKGROUND_COLOUR = "white",
    FOREGROUND_COLOUR = "#22223b",
    QUADRAT_COLOUR = "black",
    # extras
    ADVANCED_ANALYSIS = FALSE
  )

  # -------------------- read + coalesce lines --------------------
  if (!file.exists(init_file)) stop("init_file not found: ", init_file)
  raw <- readLines(init_file, warn = FALSE)

  # strip comments and blanks
  strip_comments <- function(x) sub("\\s+#.*$", "", x)
  lines <- trimws(vapply(raw, strip_comments, character(1)))
  lines <- lines[nzchar(lines)]

  # join multi-line values when RHS starts a c(...) and continues
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

    # if RHS starts a c( and parentheses aren't balanced, keep appending rows
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

  # -------------------- parse scalar/vector values --------------------
  parse_value <- function(v) {
    v <- trimws(v)
    # quoted scalars
    v <- gsub('^"(.*)"$', "\\1", v)
    v <- gsub("^'(.*)'$", "\\1", v)

    # c(...) vectors (possibly multi-line)
    if (grepl("^c\\s*\\(", v)) {
      # safe-ish eval of a vector literal
      out <- tryCatch(eval(parse(text = v)), error = function(e) e)
      if (inherits(out, "error")) stop("Could not parse vector literal: ", v, "\n", out$message)
      return(out)
    }

    # simple comma vector "a,b,c"
    if (grepl(",", v, fixed = TRUE)) {
      items <- trimws(strsplit(v, ",", fixed = TRUE)[[1]])
      # try numeric first, else character
      nums <- suppressWarnings(as.numeric(items))
      if (!anyNA(nums)) {
        return(nums)
      }
      return(items)
    }

    # logical
    if (v %in% c("TRUE", "FALSE")) {
      return(as.logical(v))
    }
    # integer / numeric
    if (grepl("^[+-]?[0-9]+$", v)) {
      return(as.integer(v))
    }
    if (grepl("^[+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?$", v)) {
      return(as.numeric(v))
    }
    # single name:value -> named scalar (numeric if possible)
    if (grepl("^[^:]+:.+$", v)) {
      kv <- strsplit(v, ":", fixed = TRUE)[[1]]
      nm <- trimws(kv[1])
      vv <- trimws(kv[2])
      num <- suppressWarnings(as.numeric(vv))
      out <- if (!is.na(num)) c(num) else c(vv)
      names(out) <- nm
      return(out)
    }
    # plain string
    v
  }

  for (kv in kvs) {
    config[[kv$key]] <- parse_value(kv$val)
  }

  # --- validate colors ----------------------------------------------------------
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

  # -------------------- coerce numerics/logicals --------------------
  num_fields <- c(
    "SEED", "N_INDIVIDUALS", "N_SPECIES", "DOMINANT_FRACTION",
    "FISHER_ALPHA", "FISHER_X",
    "SAMPLING_RESOLUTION", "ENVIRONMENTAL_NOISE",
    "MAX_CLUSTERS_DOMINANT", "CLUSTER_SPREAD_DOMINANT",
    "INTERACTION_RADIUS",
    "N_QUADRATS", "N_TRANSECTS", "N_QUADRATS_PER_TRANSECT", "TRANSECT_ANGLE",
    "VORONOI_SEED_FACTOR",
    "POINT_SIZE", "POINT_ALPHA", "QUADRAT_ALPHA"
  )
  for (f in intersect(names(config), num_fields)) {
    config[[f]] <- as.numeric(config[[f]])
  }
  log_fields <- c("ADVANCED_ANALYSIS")
  for (f in intersect(names(config), log_fields)) {
    config[[f]] <- as.logical(as.character(config[[f]]))
  }

  # -------------------- gradient table --------------------
  gs <- as.character(config$GRADIENT_SPECIES)
  ga <- as.character(config$GRADIENT_ASSIGNMENTS)
  if (length(gs) != length(ga)) {
    stop("GRADIENT_SPECIES and GRADIENT_ASSIGNMENTS must be same length.")
  }

  allowed_grads <- c("temperature", "elevation", "rainfall")
  if (!all(ga %in% allowed_grads)) {
    stop(
      "Unknown gradient(s): ", paste(setdiff(ga, allowed_grads), collapse = ", "),
      ". Allowed: temperature,elevation,rainfall."
    )
  }

  # parse possible name:value lists
  opt_raw <- .parse_named_pairs_numeric(config$GRADIENT_OPTIMA %||% 0.5)
  tol_raw <- .parse_named_pairs_numeric(config$GRADIENT_TOLERANCE %||% 0.1)

  # resolve to per-species vectors
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
    stop(
      "Cannot resolve parameter vector for gradients; supply a scalar, ",
      "a vector length S, named by species, named by gradients, or length |unique(gradients)|."
    )
  }

  opt <- .resolve_param_vector(opt_raw, gs, ga)
  tol <- .resolve_param_vector(tol_raw, gs, ga)

  clamp01 <- function(x) pmax(0, pmin(1, x))
  opt <- clamp01(opt)
  if (any(!is.finite(tol) | tol <= 0)) {
    stop("GRADIENT_TOLERANCE must be positive and finite.")
  }

  config$GRADIENT <- tibble::tibble(
    species  = gs,
    gradient = ga,
    optimum  = as.numeric(opt),
    tol      = as.numeric(tol)
  )

  config$GRADIENT$species <- as.character(config$GRADIENT$species)

  # -------------------- quadrat size materialization --------------------
  QUADRAT_SIZES <- list(
    small = c(1, 1),
    medium = c(1.5, 1.5),
    large = c(2, 2)
  )
  qs <- QUADRAT_SIZES[[as.character(config$QUADRAT_SIZE_OPTION)]]
  if (is.null(qs)) stop("Invalid QUADRAT_SIZE_OPTION. Use small|medium|large.")
  config$QUADRAT_SIZE <- qs

  set.seed(config$SEED)
  config
}


#' @title Load interactions from a separate config
#' @param interactions_file Path to interactions_init.txt (may be NULL)
#' @param n_species Integer S
#' @return list(radius=numeric, matrix=SxS numeric with dimnames)
#' @title Load interactions from a separate config
#' @param interactions_file Path to interactions_init.txt (may be NULL)
#' @param n_species Integer S
#' @return list(radius=numeric, matrix=SxS numeric with dimnames)
load_interactions <- function(interactions_file, n_species) {
  spp_names <- LETTERS[1:n_species]
  # defaults = disabled
  out_radius <- 0
  IM <- matrix(1.0,
    nrow = n_species, ncol = n_species,
    dimnames = list(spp_names, spp_names)
  )

  if (is.null(interactions_file)) {
    return(list(radius = out_radius, matrix = IM))
  }
  if (!file.exists(interactions_file)) {
    warning("interactions_file not found: ", interactions_file, " — interactions disabled.")
    return(list(radius = out_radius, matrix = IM))
  }

  # parse simple key=value lines into a **named list**
  raw <- readLines(interactions_file, warn = FALSE)
  raw <- trimws(sub("#.*$", "", raw))
  raw <- raw[nzchar(raw)]
  kv <- strsplit(raw, "=", fixed = TRUE)
  K <- toupper(trimws(vapply(kv, `[`, "", 1)))
  V <- trimws(vapply(kv, function(x) paste(x[-1], collapse = "="), ""))
  params <- as.list(stats::setNames(V, K)) # <— list, not atomic vector!

  # small helper for safe gets
  getp <- function(nm, default = NULL) {
    val <- params[[nm]]
    if (is.null(val) || !nzchar(val)) default else val
  }

  # radius
  rstr <- getp("INTERACTION_RADIUS")
  if (!is.null(rstr)) {
    rnum <- suppressWarnings(as.numeric(rstr))
    if (length(rnum) == 1 && is.finite(rnum) && rnum >= 0) {
      out_radius <- rnum
    } else {
      warning("Bad INTERACTION_RADIUS in interactions file; using 0 (disabled).")
      out_radius <- 0
    }
  }

  # choose source: MATRIX_CSV or EDGELIST_CSV or AUTO
  matrix_csv <- getp("MATRIX_CSV")
  edgelist_csv <- getp("EDGELIST_CSV")
  auto <- tolower(getp("AUTO", "false")) %in% c("true", "1", "yes")

  if (!is.null(matrix_csv)) {
    if (!file.exists(matrix_csv)) stop("MATRIX_CSV not found: ", matrix_csv)
    M <- as.matrix(utils::read.csv(matrix_csv, row.names = 1, check.names = FALSE))
    # coerce order and fill missing with 1
    fullM <- matrix(1.0,
      nrow = n_species, ncol = n_species,
      dimnames = list(spp_names, spp_names)
    )
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
      if (f %in% spp_names && n %in% spp_names && is.finite(v)) {
        IM[f, n] <- v
      }
    }
  } else if (auto) {
    if (all(c("D", "E") %in% spp_names)) IM["E", "D"] <- 3.0
  } # else: leave as all 1s

  if (any(!is.finite(IM))) stop("INTERACTION_MATRIX contains non-finite values.")
  list(radius = out_radius, matrix = IM)
}



# -----------------------------------------------------------------------------
# 3. CORE SIMULATION FUNCTIONS
# -----------------------------------------------------------------------------

#' @title Assign Individuals by Clustering (Efficiently)
#' @description Helper function to assign individuals to locations based on proximity
#'   to random cluster centers. Uses a fast nearest-neighbour search (k-d tree)
#'   which is significantly faster than calculating a full distance matrix.
#' @param points_with_env The full sf data frame of all potential point locations.
#' @param available_indices Numeric vector of indices in `points_with_env` that are available for assignment.
#' @param n_ind_to_assign The number of individuals of the current species to assign.
#' @param n_clusters The number of clusters to form.
#' @param spread The spread parameter (sigma) for the exponential decay function.
#' @return A numeric vector of indices selected from `available_indices`.
assign_by_clustering_fast <- function(points_with_env, available_indices, n_ind_to_assign, n_clusters, spread) {
  # 1. Select cluster centers from available locations
  # Ensure we don't try to pick more cluster centers than available points
  n_clusters <- min(n_clusters, length(available_indices))
  if (n_clusters == 0) {
    return(integer(0))
  } # Return empty if no clusters can be formed
  cluster_centers_idx <- sample(available_indices, n_clusters)

  # 2. Extract coordinates for fast nearest-neighbour search
  # Data points are the cluster centers
  center_coords <- st_coordinates(points_with_env[cluster_centers_idx, ])
  # Query points are all available locations
  available_coords <- st_coordinates(points_with_env[available_indices, ])

  # 3. Perform fast nearest-neighbour search using RANN::nn2
  # nn2 returns a list with distances (nn.dists) and indices (nn.idx)
  # We only need the distance to the single nearest neighbour (k=1)
  nn_results <- RANN::nn2(data = center_coords, query = available_coords, k = 1)

  # The result is a matrix, we extract the distance column as a vector
  min_dists <- nn_results$nn.dists[, 1]

  # 4. Calculate sampling probabilities based on distance
  # Points closer to a cluster center have a higher probability of being selected
  probs <- exp(-min_dists / spread)

  # Handle cases where all probabilities might be zero (e.g., very large distances and small spread)
  if (all(probs == 0)) {
    probs <- rep(1, length(available_indices))
  }

  # 5. Sample from the available indices based on calculated probabilities
  selected_indices <- sample(available_indices, size = n_ind_to_assign, prob = probs)

  return(selected_indices)
}


#' @title Create Irregular Sampling Domain
#' @description Generates a vaguely organic, amoeba-like polygon to serve as the study area.
#' @return An `sf` object containing a single polygon.
create_sampling_domain <- function() {
  theta <- seq(0, 2 * pi, length.out = 20)
  r <- 10 + 3 * sin(3 * theta) + 2 * cos(5 * theta) + runif(20, -1, 1)
  x <- r * cos(theta) + runif(20, -0.5, 0.5)
  y <- r * sin(theta) * 0.8 + runif(20, -0.5, 0.5)
  coords <- cbind(x, y)
  coords <- rbind(coords, coords[1, ]) # Close the polygon
  polygon <- st_polygon(list(coords))
  return(st_sf(geometry = st_sfc(polygon)))
}


#' @title Generate Fisher's Log-Series Abundances
#' @description Creates a species abundance distribution (SAD) following Fisher's log-series.
#'   A single dominant species is set manually, and the rest follow the series.
#' @param n_species Total number of species.
#' @param n_individuals Total number of individuals.
#' @param dominant_fraction The proportion of total individuals belonging to the dominant species.
#' @param alpha Fisher's alpha, a diversity parameter.
#' @param x A parameter of the log-series, typically close to 1.
#' @return A named numeric vector of species abundances.
generate_fisher_log_series <- function(n_species, n_individuals, dominant_fraction, alpha, x) {
  n_dominant <- round(n_individuals * dominant_fraction)
  n_remaining <- n_individuals - n_dominant

  # Abundances for non-dominant species
  ranks <- 2:n_species
  rel_abundances <- alpha * (x^ranks) / ranks
  abundances <- round(rel_abundances / sum(rel_abundances) * n_remaining)

  all_abundances <- c(n_dominant, abundances)

  # Adjust to ensure the total number of individuals is exact
  adjustment <- n_individuals - sum(all_abundances)
  all_abundances[2] <- all_abundances[2] + adjustment

  species_names <- LETTERS[1:n_species]
  names(all_abundances) <- species_names
  return(all_abundances[all_abundances > 0]) # Return only species with >0 individuals
}


#' @title Create Environmental Gradient Fields
#' @description Generates spatial environmental gradients (temperature, elevation, rainfall)
#'   across a grid within the domain.
#' @param domain The `sf` polygon object of the study area.
#' @param resolution The number of grid cells along each axis.
#' @param noise_level The standard deviation of random noise to add to gradients.
#' @return A data frame representing a grid with raw and rescaled environmental values.
create_environmental_gradients <- function(domain, resolution, noise_level) {
  bbox <- st_bbox(domain)
  x_seq <- seq(bbox["xmin"], bbox["xmax"], length.out = resolution)
  y_seq <- seq(bbox["ymin"], bbox["ymax"], length.out = resolution)
  grid <- expand.grid(x = x_seq, y = y_seq)

  # Normalize coordinates
  x_norm <- (grid$x - bbox["xmin"]) / (bbox["xmax"] - bbox["xmin"])
  y_norm <- (grid$y - bbox["ymin"]) / (bbox["ymax"] - bbox["ymin"])

  # Define gradients (normalized 0-1)
  temp_vals <- x_norm * 0.7 + y_norm * 0.3 + rnorm(nrow(grid), 0, noise_level)

  dist_center <- sqrt((x_norm - 0.5)^2 + (y_norm - 0.5)^2)
  elev_vals <- 1 - (dist_center / sqrt(0.5^2 + 0.5^2)) + rnorm(nrow(grid), 0, noise_level)

  rain_vals <- (-x_norm * 0.6 + y_norm * 0.8)
  rain_vals <- (rain_vals - min(rain_vals)) / (max(rain_vals) - min(rain_vals)) + rnorm(nrow(grid), 0, noise_level)

  # Clamp values to [0, 1] range
  grid$temperature <- pmax(0, pmin(1, temp_vals))
  grid$elevation <- pmax(0, pmin(1, elev_vals))
  grid$rainfall <- pmax(0, pmin(1, rain_vals))

  # Rescale to more realistic units for reporting
  grid$temperature_C <- grid$temperature * 30 - 2 # Range: -2 to 28 C
  grid$elevation_m <- grid$elevation * 2000 # Range: 0 to 2000 m
  grid$rainfall_mm <- grid$rainfall * 700 + 200 # Range: 200 to 900 mm

  return(grid)
}


#' @title Generate Heterogeneous Species Distribution with Interspecific Interactions
#' @description Places individuals based on abundance, clustering, environmental
#'   gradients, and interspecific interactions (competition/facilitation).
#' @param domain The `sf` polygon of the study area.
#' @param P The configuration list, including INTERACTION_MATRIX and INTERACTION_RADIUS.
#' @return An `sf` point object with each point assigned a species.
generate_heterogeneous_distribution <- function(domain, P) {
  # --- Initial Setup ---
  abundances <- generate_fisher_log_series(
    P$N_SPECIES, P$N_INDIVIDUALS, P$DOMINANT_FRACTION, P$FISHER_ALPHA, P$FISHER_X
  )
  env_grid <- create_environmental_gradients(domain, P$SAMPLING_RESOLUTION, P$ENVIRONMENTAL_NOISE)
  env_sf <- st_as_sf(env_grid, coords = c("x", "y"), crs = st_crs(domain))
  all_points <- st_sample(domain, size = P$N_INDIVIDUALS, type = "random")
  points_with_env <- st_join(st_sf(geometry = all_points), env_sf, join = st_nearest_feature)
  points_with_env$species <- ""
  available_indices <- 1:nrow(points_with_env)

  # Main placement loop
  for (sp in names(abundances)) {
    n_ind <- abundances[sp]
    if (length(available_indices) < n_ind) n_ind <- length(available_indices)
    if (n_ind == 0) next

    # 1) Base probability (environmental response OR dominant clustering)
    base_probs <- rep(1.0, length(available_indices))

    if (sp %in% P$GRADIENT$species) {
      row <- P$GRADIENT[P$GRADIENT$species == sp, ][1, ]
      gradient_type <- row$gradient
      optimum <- row$optimum
      tol <- row$tol

      env_values <- points_with_env[[gradient_type]][available_indices]
      base_probs <- exp(-((env_values - optimum)^2) / (2 * tol^2))
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
        if (!all(is.finite(base_probs)) || all(base_probs <= 0) || anyNA(base_probs)) {
          base_probs <- rep(1.0, n_avail)
        }
      } else {
        base_probs <- rep(1.0, n_avail)
      }
    }

    # 2) Inter-specific interaction modifier
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
            iv <- pmax(iv, 1e-12) # guard: avoid log(0) while allowing near-exclusion

            exp(mean(log(iv)))
          }, numeric(1))

          interaction_modifier <- interaction_scores
        }
      }
    }

    # 3) Final probability and assignment
    final_probs <- base_probs * interaction_modifier
    if (all(!is.finite(final_probs)) || all(final_probs <= 0)) {
      final_probs <- rep(1, length(available_indices)) # guard: only if everything is degenerate
    }

    selected <- sample(available_indices, size = n_ind, prob = final_probs, replace = FALSE)
    if (length(selected) > 0) {
      points_with_env$species[selected] <- sp
      available_indices <- setdiff(available_indices, selected)
    }
  }

  points_with_env %>% dplyr::filter(species != "")
}


#' @title Create a Quadrat Polygon from its Center
#' @description Helper function to generate an axis-aligned rectangular polygon
#'   from a center point and dimensions.
#' @param center_point An `sf` point geometry object.
#' @param size A numeric vector `c(width, height)`.
#' @return An `st_polygon` geometry object.
create_quadrat_from_center <- function(center_point, size) {
  half_w <- size[1] / 2
  half_h <- size[2] / 2
  # st_coordinates returns a matrix, so we get the first row
  coords <- st_coordinates(center_point)[1, ]
  st_polygon(list(cbind(
    c(coords["X"] - half_w, coords["X"] + half_w, coords["X"] + half_w, coords["X"] - half_w, coords["X"] - half_w),
    c(coords["Y"] - half_h, coords["Y"] - half_h, coords["Y"] + half_h, coords["Y"] + half_h, coords["Y"] - half_h)
  )))
}


#' @title Place Non-Overlapping Quadrats
#' @description Randomly places rectangular quadrats within the domain, ensuring they
#'   are fully contained and do not overlap.
#' @param domain The `sf` polygon of the study area.
#' @param n_quadrats The target number of quadrats.
#' @param quadrat_size A numeric vector `c(width, height)` for the quadrats.
#' @return An `sf` object containing the quadrat polygons.
place_quadrats <- function(domain, n_quadrats, quadrat_size) {
  bbox <- st_bbox(domain)
  quadrats <- list()
  attempts <- 0
  max_attempts <- n_quadrats * 100

  while (length(quadrats) < n_quadrats && attempts < max_attempts) {
    attempts <- attempts + 1

    x_center <- runif(1, bbox["xmin"] + quadrat_size[1] / 2, bbox["xmax"] - quadrat_size[1] / 2)
    y_center <- runif(1, bbox["ymin"] + quadrat_size[2] / 2, bbox["ymax"] - quadrat_size[2] / 2)

    center_pt_sfc <- st_sfc(st_point(c(x_center, y_center)), crs = st_crs(domain))
    new_quadrat_poly <- create_quadrat_from_center(center_pt_sfc, quadrat_size) %>%
      st_sfc(crs = st_crs(domain))

    # NOTE: st_within(..., sparse = FALSE) returns a matrix; select the first (only) column
    within_mat <- st_within(new_quadrat_poly, domain, sparse = FALSE)
    is_within <- isTRUE(within_mat[, 1, drop = TRUE][1])

    if (is_within) {
      is_overlapping <- FALSE
      if (length(quadrats) > 0) {
        existing_quadrats_sfc <- do.call(c, quadrats)
        # Returns a logical matrix [1 x nExisting]; reduce it explicitly
        overlap_mat <- st_intersects(new_quadrat_poly, existing_quadrats_sfc, sparse = FALSE)
        is_overlapping <- any(overlap_mat[1, , drop = TRUE])
      }
      if (!is_overlapping) {
        quadrats[[length(quadrats) + 1]] <- new_quadrat_poly
      }
    }
  }

  if (length(quadrats) < n_quadrats) {
    warning(sprintf(
      "Placed only %d of %d requested quadrats after %d attempts.",
      length(quadrats), n_quadrats, max_attempts
    ))
  }

  if (length(quadrats) == 0) {
    return(st_sf(quadrat_id = integer(0), geometry = st_sfc(crs = st_crs(domain))))
  }

  final_quadrats <- do.call(c, quadrats) %>%
    st_sf(quadrat_id = 1:length(.), geometry = .)
  return(final_quadrats)
}


#' @title Place Non-Overlapping Quadrats (Systematic Tiling Method)
#' @description Places rectangular quadrats within the domain using an efficient
#'   tiling method. It creates a grid of all possible non-overlapping quadrat
#'   locations and randomly samples from those that are fully contained
#'   within the domain.
#' @param domain The `sf` polygon of the study area.
#' @param n_quadrats The target number of quadrats.
#' @param quadrat_size A numeric vector `c(width, height)` for the quadrats.
#' @return An `sf` object containing the quadrat polygons.
place_quadrats_tiled <- function(domain, n_quadrats, quadrat_size) {
  # 1. Create a grid of potential non-overlapping quadrats (cells)
  candidate_grid <- st_make_grid(domain, cellsize = quadrat_size, what = "polygons")

  # 2. Filter the grid to those fully within the domain
  within_mat <- st_within(candidate_grid, domain, sparse = FALSE)
  inside <- if (is.matrix(within_mat)) drop(within_mat[, 1, drop = TRUE]) else as.logical(within_mat)
  valid_locations <- candidate_grid[inside]

  num_possible <- length(valid_locations)

  # 3. Check if enough valid locations were found
  if (num_possible == 0) {
    warning("Systematic placement failed: No quadrats of the given size can fit entirely within the domain.")
    return(st_sf(quadrat_id = integer(0), geometry = st_sfc(crs = st_crs(domain))))
  }

  if (num_possible < n_quadrats) {
    warning(sprintf(
      "Could only place %d of %d requested quadrats, as this is the maximum that can fit.",
      num_possible, n_quadrats
    ))
    n_quadrats <- num_possible
  }

  # 4. Randomly sample from the valid locations
  sampled_indices <- sample.int(num_possible, size = n_quadrats)
  final_quadrats_sfc <- valid_locations[sampled_indices]

  # 5. Finalize the sf object with IDs
  st_sf(quadrat_id = seq_len(n_quadrats), geometry = final_quadrats_sfc)
}



#' @title Place Quadrats Using Voronoi Tessellation
#' @description Places quadrats by generating a Voronoi diagram from random seed
#'    points. It identifies cells large enough to contain a quadrat and places
#'    the quadrat at the cell's "pole of inaccessibility" (center of the
#'    largest inscribed circle), ensuring quadrats are well-spaced.
#' @param domain The `sf` polygon of the study area.
#' @param n_quadrats The target number of quadrats.
#' @param quadrat_size A numeric vector `c(width, height)`.
#' @return An `sf` object containing the quadrat polygons.
place_quadrats_voronoi <- function(domain, n_quadrats, quadrat_size, voronoi_seed_factor) {
  # 1. Generate more seed points than needed to ensure good coverage and choice.
  n_seeds <- n_quadrats * voronoi_seed_factor

  # Ensure the domain is a single polygon for st_voronoi to work correctly
  domain_union <- st_union(domain)
  seed_points <- st_sample(domain_union, size = n_seeds, type = "random")

  # 2. Create the Voronoi tessellation and clip it to the domain boundary.
  voronoi_polys <- st_voronoi(st_union(seed_points))
  voronoi_clipped <- st_intersection(st_cast(voronoi_polys), domain_union)

  # 3. For each Voronoi cell, find the largest possible circle that fits inside.
  # We suppress warnings as it can be noisy on complex geometries.
  inscribed_circles <- suppressWarnings(st_inscribed_circle(voronoi_clipped))

  # 4. Determine which cells are large enough to hold our quadrat.
  quadrat_half_diag <- sqrt(quadrat_size[1]^2 + quadrat_size[2]^2) / 2
  radii <- sqrt(st_area(inscribed_circles) / pi)
  suitable_indices <- which(radii >= quadrat_half_diag)

  # 5. Check if we found enough suitable locations.
  num_possible <- length(suitable_indices)
  if (num_possible == 0) {
    warning("Voronoi placement failed: No cells were large enough to contain a quadrat.")
    return(st_sf(quadrat_id = integer(0), geometry = st_sfc(crs = st_crs(domain))))
  }

  if (num_possible < n_quadrats) {
    warning(sprintf(
      "Voronoi placement could only find %d suitable locations out of %d requested.",
      num_possible, n_quadrats
    ))
    n_quadrats <- num_possible
  }

  # 6. Randomly sample from the suitable locations.
  sampled_indices <- sample(suitable_indices, size = n_quadrats)
  final_centers <- st_centroid(inscribed_circles[sampled_indices])

  # 7. Create the final axis-aligned quadrats using the helper function.
  quadrat_list <- lapply(st_geometry(final_centers), function(pt) {
    create_quadrat_from_center(pt, quadrat_size)
  })

  final_quadrats_sfc <- st_sfc(quadrat_list, crs = st_crs(domain))
  final_quadrats <- st_sf(quadrat_id = 1:length(final_quadrats_sfc), geometry = final_quadrats_sfc)

  return(final_quadrats)
}


#' @title Place Quadrats in a Systematic Grid
#' @description Generates a regular grid of points over the domain and places a
#'   quadrat at each point that falls completely within the domain.
#' @param domain The `sf` polygon of the study area.
#' @param n_quadrats The target number of quadrats (approximate).
#' @param quadrat_size A numeric vector `c(width, height)`.
#' @return An `sf` object containing the quadrat polygons.
place_quadrats_systematic <- function(domain, n_quadrats, quadrat_size) {
  bbox <- st_bbox(domain)
  aspect_ratio <- (bbox["ymax"] - bbox["ymin"]) / (bbox["xmax"] - bbox["xmin"])
  nx <- round(sqrt(n_quadrats / aspect_ratio))
  ny <- round(aspect_ratio * nx)

  candidate_centers <- st_make_grid(domain, n = c(nx, ny), what = "centers")

  candidate_quadrats <- st_sfc(lapply(candidate_centers, function(pt) {
    create_quadrat_from_center(pt, quadrat_size)
  }), crs = st_crs(domain))

  # Already correct: select first (only) domain column from logical matrix
  within_mat <- st_within(candidate_quadrats, domain, sparse = FALSE)
  is_within <- within_mat[, 1, drop = TRUE]
  valid_quadrats <- candidate_quadrats[is_within]

  if (length(valid_quadrats) == 0) {
    warning("Systematic sampling failed to place any quadrats within the domain.")
    return(st_sf(quadrat_id = integer(0), geometry = st_sfc(crs = st_crs(domain))))
  }

  st_sf(quadrat_id = seq_len(length(valid_quadrats)), geometry = valid_quadrats)
}


#' @title Place Quadrats Along Parallel Transects (First Principles Rewrite)
#' @description This is a complete rewrite using basic geometry and linear
#'    interpolation, avoiding advanced sf functions to bypass a persistent bug.
#' @param domain The `sf` polygon of the study area.
#' @param n_transects The number of parallel transect lines.
#' @param n_quadrats_per_transect The number of quadrats on each transect.
#' @param quadrat_size A numeric vector `c(width, height)`.
#' @param angle The compass direction of the transects in degrees (0=N, 90=E).
#' @return An `sf` object containing the quadrat polygons.
place_quadrats_transect <- function(domain, n_transects, n_quadrats_per_transect, quadrat_size, angle) {
  # 1. Create a "safe sampling area".
  buffer_dist <- sqrt(quadrat_size[1]^2 + quadrat_size[2]^2) / 2
  safe_domain <- st_buffer(domain, -buffer_dist)

  if (st_is_empty(safe_domain) || st_area(safe_domain) == 0) {
    warning("Domain is too small to create a safe sampling area for the given quadrat size.")
    return(st_sf(quadrat_id = integer(0), geometry = st_sfc(crs = st_crs(domain))))
  }

  # 2. Define the ideal mathematical transect lines.
  bbox <- st_bbox(domain)
  center <- st_centroid(st_as_sfc(bbox))
  diag_len <- sqrt((bbox["xmax"] - bbox["xmin"])^2 + (bbox["ymax"] - bbox["ymin"])^2) * 1.5

  y_coords <- seq(bbox["ymin"], bbox["ymax"], length.out = n_transects + 2)[2:(n_transects + 1)]

  # 3. Rotate these lines to the desired angle.
  math_angle_deg <- 90 - angle
  math_angle_rad <- math_angle_deg * pi / 180
  rot_matrix <- matrix(c(cos(math_angle_rad), sin(math_angle_rad), -sin(math_angle_rad), cos(math_angle_rad)), 2, 2)

  horizontal_lines <- st_sfc(lapply(y_coords, function(y) {
    st_linestring(matrix(c(st_coordinates(center)[1] - diag_len / 2, st_coordinates(center)[1] + diag_len / 2, y, y), ncol = 2))
  }), crs = st_crs(domain))
  rotated_lines <- (horizontal_lines - center) * rot_matrix + center
  st_crs(rotated_lines) <- st_crs(domain)

  points_list <- list()
  for (i in seq_len(n_transects)) {
    line <- rotated_lines[i]
    segment <- st_intersection(line, safe_domain)

    if (st_is_empty(segment) || st_length(segment) == 0) next

    # Manually get the start and end coordinates of the valid segment.
    segment_coords <- st_coordinates(segment)
    start_pt <- segment_coords[1, c("X", "Y")]
    end_pt <- segment_coords[nrow(segment_coords), c("X", "Y")]

    # Manually calculate N points using linear interpolation.
    if (n_quadrats_per_transect > 1) {
      fractions <- seq(0, 1, length.out = n_quadrats_per_transect)
    } else {
      fractions <- 0.5 # Place a single point in the middle
    }

    x_coords <- start_pt[1] + fractions * (end_pt[1] - start_pt[1])
    y_coords <- start_pt[2] + fractions * (end_pt[2] - start_pt[2])

    points_list[[i]] <- st_as_sf(data.frame(x = x_coords, y = y_coords), coords = c("x", "y"), crs = st_crs(domain))
  }

  if (length(points_list) == 0) {
    warning("No transects intersected the safe sampling area.")
    return(st_sf(quadrat_id = integer(0), geometry = st_sfc(crs = st_crs(domain))))
  }

  candidate_centers <- do.call(rbind, points_list)

  # Create axis-aligned quadrats around the calculated centers using the helper.
  quadrat_geometries <- st_sfc(lapply(st_geometry(candidate_centers), function(pt) {
    create_quadrat_from_center(pt, quadrat_size)
  }), crs = st_crs(domain))

  final_quadrats <- st_sf(quadrat_id = 1:length(quadrat_geometries), geometry = quadrat_geometries)

  return(final_quadrats)
}


# -----------------------------------------------------------------------------
# 4. ADVANCED ANALYSES
# -----------------------------------------------------------------------------

#' @title Calculate Rank-Abundance Data (Base R Method)
#' @description Computes both the observed rank-abundance distribution from the
#'    simulation results and the theoretical distribution from Fisher's log-series
#'    parameters. This version uses base R's `table()` to avoid dplyr conflicts.
#' @param species_dist An `sf` object of all individuals with a `species` column.
#' @param P The configuration list for the simulation.
#' @return A single data frame containing ranked abundance data for both sources.
calculate_rank_abundance <- function(species_dist, P) {
  # 1. Calculate OBSERVED data using base R's table() to avoid n() error
  observed_counts <- table(species_dist$species)
  observed_data <- as.data.frame(observed_counts, stringsAsFactors = FALSE)
  names(observed_data) <- c("Species", "Abundance")

  # We can still use dplyr for safe operations like arranging and ranking
  observed_data <- observed_data %>%
    arrange(desc(Abundance)) %>%
    mutate(
      Rank = row_number(),
      Source = "Observed"
    ) %>%
    select(Rank, Abundance, Source)

  # 2. Calculate THEORETICAL data (this part was already safe)
  theoretical_abundances <- generate_fisher_log_series(
    P$N_SPECIES, P$N_INDIVIDUALS, P$DOMINANT_FRACTION, P$FISHER_ALPHA, P$FISHER_X
  )

  theoretical_data <- tibble(
    Abundance = sort(as.numeric(theoretical_abundances), decreasing = TRUE),
    Rank = seq_along(theoretical_abundances), # or 1:length(theoretical_abundances)
    Source = "Theoretical"
  )

  # 3. Combine and return
  bind_rows(observed_data, theoretical_data)
}


#' @title Plot Rank-Abundance Curve
#' @description Generates a ggplot object showing the rank-abundance curve(s).
#' @param rank_abundance_data A data frame produced by `calculate_rank_abundance`.
#' @return A `ggplot` object.
plot_rank_abundance <- function(rank_abundance_data) {
  ggplot(rank_abundance_data, aes(x = Rank, y = Abundance, color = Source)) +
    geom_line(aes(linetype = Source), linewidth = 1.1) +
    geom_point(aes(shape = Source), size = 3, fill = "white", stroke = 1.2) +
    # This is the original fix for the labels argument
    scale_y_log10(
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_color_manual(values = c("Observed" = "black", "Theoretical" = "#e41a1c")) +
    scale_linetype_manual(values = c("Observed" = "solid", "Theoretical" = "dashed")) +
    scale_shape_manual(values = c("Observed" = 21, "Theoretical" = 22)) +
    labs(
      title = "Species-Abundance Distribution (SAD)",
      subtitle = "Comparison of observed data vs. theoretical Fisher's log-series model",
      x = "Species Rank (Most to Least Abundant)",
      y = "Abundance (Log Scale)",
      color = "Distribution",
      linetype = "Distribution",
      shape = "Distribution"
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
}


#' @title Calculate Occupancy-Abundance Data
#' @description Calculates the total abundance and occupancy (number of sites)
#'   for each species from a site-by-species abundance matrix.
#' @param abund_matrix A data frame with sites in rows and species in columns,
#'   containing abundance counts. The first column should be named "site".
#' @return A data frame with columns for Species, TotalAbundance, and Occupancy.
calculate_occupancy_abundance <- function(abund_matrix) {
  # Exclude the 'site' column for calculations using base R indexing
  abund_numeric <- abund_matrix[, -which(names(abund_matrix) == "site")]

  # Calculate total abundance and occupancy for each species (column)
  oa_data <- data.frame(
    Species = names(abund_numeric),
    TotalAbundance = colSums(abund_numeric),
    Occupancy = colSums(abund_numeric > 0),
    row.names = NULL # Ensure row names are not carried over
  )

  return(oa_data)
}

#' @title Plot Occupancy-Abundance Relationship
#' @description Generates a ggplot object showing the relationship between
#'   species abundance and occupancy, typically on a log-log scale.
#' @param oa_data A data frame produced by `calculate_occupancy_abundance`.
#' @return A `ggplot` object.
plot_occupancy_abundance <- function(oa_data) {
  # Add a check for empty data
  if (nrow(oa_data) == 0 || all(oa_data$TotalAbundance == 0)) {
    return(
      ggplot() +
        labs(
          title = "Occupancy-Abundance Relationship",
          subtitle = "No data to plot."
        ) +
        theme_void()
    )
  }

  ggplot(oa_data, aes(x = TotalAbundance, y = Occupancy)) +
    geom_point(alpha = 0.6, size = 3, color = "#2c7fb8") +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed", linewidth = 0.8) +
    scale_x_log10(labels = scales::label_log()) +
    scale_y_log10(labels = scales::label_log()) +
    labs(
      title = "Occupancy-Abundance Relationship",
      subtitle = "Widespread species tend to be more abundant",
      x = "Total Abundance (Log Scale)",
      y = "Number of Sites Occupied (Log Scale)"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold")
    )
}


#' @title Calculate Species-Area Relationship Data
#' @description Computes a species accumulation curve using the 'random' method
#'   from `vegan::specaccum`.
#' @param abund_matrix A site-by-species data frame of species abundances.
#'   The first column should be named "site".
#' @return A data frame with columns for Sites, Richness (mean), and SD (standard deviation).
calculate_species_area <- function(abund_matrix) {
  # Exclude the 'site' column for calculations
  abund_numeric <- abund_matrix[, -which(names(abund_matrix) == "site")]

  # Calculate the species accumulation curve using 100 random permutations
  sar_curve <- vegan::specaccum(abund_numeric, method = "random", permutations = 100)

  # Extract the results into a clean data frame for plotting
  sar_data <- data.frame(
    Sites = sar_curve$sites,
    Richness = sar_curve$richness,
    SD = sar_curve$sd
  )

  return(sar_data)
}


#' @title Plot Species-Area Relationship
#' @description Generates a ggplot object showing the species accumulation curve.
#' @param sar_data A data frame produced by `calculate_species_area`.
#' @return A `ggplot` object.
plot_species_area <- function(sar_data) {
  ggplot(sar_data, aes(x = Sites, y = Richness)) +
    # Shaded region represents the standard deviation
    geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), fill = "#41b6c4", alpha = 0.4) +
    geom_line(color = "#08519c", linewidth = 1.2) +
    labs(
      title = "Species-Area Relationship (SAR)",
      subtitle = "Cumulative species richness as sampling area increases",
      x = "Number of Quadrats (Area)",
      y = "Cumulative Number of Species"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold")
    )
}


#' @title Calculate Distance-Decay Data
#' @description Computes pairwise geographic distances and community dissimilarities
#'   between all sites.
#' @param abund_matrix A site-by-species data frame of species abundances.
#' @param site_coords A data frame with site coordinates (must contain 'x' and 'y' columns).
#' @return A data frame with 'Distance' and 'Dissimilarity' columns for plotting.
calculate_distance_decay <- function(abund_matrix, site_coords) {
  # Get numeric coordinates and abundance matrices
  coords <- site_coords[, c("x", "y")]
  abund_numeric <- abund_matrix[, -which(names(abund_matrix) == "site")]

  # Calculate Euclidean geographic distance matrix
  geo_dist <- dist(coords, method = "euclidean")

  # Calculate Sørensen dissimilarity matrix
  # This is equivalent to Bray-Curtis on presence-absence data
  comm_dissim <- vegan::vegdist(abund_numeric, method = "bray", binary = TRUE)

  # Combine the two distance vectors into a single data frame
  decay_data <- data.frame(
    Distance = as.vector(geo_dist),
    Dissimilarity = as.vector(comm_dissim)
  )

  return(decay_data)
}


#' @title Plot Distance-Decay Relationship
#' @description Generates a scatter plot of community dissimilarity versus
#'   geographic distance.
#' @param decay_data A data frame produced by `calculate_distance_decay`.
#' @return A `ggplot` object.
plot_distance_decay <- function(decay_data) {
  ggplot(decay_data, aes(x = Distance, y = Dissimilarity)) +
    # Use semi-transparent points to handle overplotting
    geom_point(alpha = 0.3, shape = 16, color = "#253494") +
    geom_smooth(method = "loess", color = "#e31a1c", se = TRUE, linewidth = 1.1) +
    # Ensure y-axis is scaled from 0 to 1
    ylim(0, 1) +
    labs(
      title = "Distance-Decay of Community Similarity",
      subtitle = "Community similarity decreases as geographic distance increases",
      x = "Geographic Distance Between Quadrats",
      y = "Community Dissimilarity (Sørensen Index)"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold")
    )
}


#' @title Calculate Rarefaction Curves
#' @description Computes species rarefaction curves for each site (row) in an
#'   abundance matrix.
#' @param abund_matrix A site-by-species data frame of species abundances.
#' @return A tidy data frame with columns for SiteID, SampleSize, and
#'   RarefiedRichness, suitable for plotting with ggplot2.
calculate_rarefaction <- function(abund_matrix) {
  # Exclude the 'site' column for calculations
  abund_numeric <- abund_matrix[, -which(names(abund_matrix) == "site")]
  site_ids <- abund_matrix$site

  # vegan::rarecurve calculates the expected richness for each step in sample size
  # It returns a list, with one element per site
  rarefaction_list <- vegan::rarecurve(abund_numeric, step = 1)

  # Process the list into a single, tidy data frame
  output_list <- list()
  for (i in seq_along(rarefaction_list)) {
    # The vector of rarefied richness values for the current site
    richness_values <- rarefaction_list[[i]]
    # The corresponding sample sizes are stored in the 'Subsample' attribute
    sample_sizes <- attr(richness_values, "Subsample")

    output_list[[i]] <- data.frame(
      SiteID = as.factor(site_ids[i]),
      SampleSize = sample_sizes,
      RarefiedRichness = richness_values
    )
  }

  # Combine the list of data frames into one
  rarefaction_data <- do.call(rbind, output_list)

  return(rarefaction_data)
}


#' @title Plot Rarefaction Curves
#' @description Generates a ggplot object showing rarefaction curves for all sites.
#' @param rarefaction_data A data frame produced by `calculate_rarefaction`.
#' @return A `ggplot` object.
plot_rarefaction <- function(rarefaction_data) {
  # Calculate the max richness to set a consistent color scale
  max_richness <- max(rarefaction_data$RarefiedRichness, na.rm = TRUE)

  ggplot(rarefaction_data, aes(x = SampleSize, y = RarefiedRichness, group = SiteID, color = SiteID)) +
    geom_line(linewidth = 0.8) +
    scale_color_viridis_d(option = "plasma") +
    labs(
      title = "Rarefaction Curves for Each Site",
      subtitle = "Comparing species richness at equivalent sample sizes",
      x = "Number of Individuals Sampled (Sample Size)",
      y = "Expected Number of Species (Rarefied Richness)",
      color = "Quadrat ID"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
}


#' @title Generate Advanced Analysis Panel
#' @description Gathers all 2D advanced analysis plots and arranges them into a
#'   single panel using patchwork.
#' @param res A list containing all simulation artifacts (data frames, stats, etc.).
#' @return A combined `patchwork` ggplot object.
generate_advanced_panel <- function(res) {
  # --- 1. Generate all individual plot objects ---

  # To make the panel cleaner, we'll reduce font sizes and remove some legends
  theme_panel <- theme(text = element_text(size = 11), legend.title = element_text(size = 10), legend.text = element_text(size = 9))

  p_rank <- plot_rank_abundance(calculate_rank_abundance(res$species_dist, res$P)) +
    labs(subtitle = NULL) + theme_panel + theme(legend.position = "bottom")

  p_oa <- plot_occupancy_abundance(calculate_occupancy_abundance(res$abund_matrix)) +
    labs(subtitle = NULL) + theme_panel

  p_sar <- plot_species_area(calculate_species_area(res$abund_matrix)) +
    labs(subtitle = NULL) + theme_panel

  p_decay <- plot_distance_decay(calculate_distance_decay(res$abund_matrix, res$site_coords)) +
    labs(subtitle = NULL) + theme_panel

  p_rare <- plot_rarefaction(calculate_rarefaction(res$abund_matrix)) +
    labs(subtitle = NULL) + theme_panel + guides(color = "none") # Legend too busy for panel

  # --- 2. Combine plots using patchwork ---

  final_panel <- (p_rank | p_oa) /
    (p_sar | p_decay) /
    (p_rare | patchwork::plot_spacer()) +
    plot_annotation(
      title = "Advanced Ecological Analysis Panel",
      theme = theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
    )

  return(final_panel)
}


# -----------------------------------------------------------------------------
# 5. ANALYSIS & DATA EXTRACTION FUNCTIONS
# -----------------------------------------------------------------------------

#' @title Create Site-by-Species Abundance Matrix
#' @description Calculates the abundance of each species in each quadrat.
#' @param species_dist An `sf` object of species points.
#' @param quadrats An `sf` object of quadrat polygons.
#' @param all_species_names A character vector of all possible species names for columns.
#' @return A data frame with `site` (quadrat_id) and columns for each species' abundance.
create_abundance_matrix <- function(species_dist, quadrats, all_species_names) {
  # Perform a single spatial intersection for all quadrats
  intersections <- st_intersection(species_dist, quadrats)

  if (nrow(intersections) == 0) {
    # Return an empty matrix if no individuals were sampled
    abund_df <- data.frame(site = quadrats$quadrat_id)
    for (sp in all_species_names) abund_df[[sp]] <- 0
    return(abund_df)
  }

  # Count species per quadrat and pivot to wide format
  abund_df <- intersections %>%
    st_drop_geometry() %>%
    count(quadrat_id, species, name = "abundance") %>%
    pivot_wider(names_from = species, values_from = abundance, values_fill = 0)

  # Ensure all species are present as columns and join with all quadrat IDs
  missing_species <- setdiff(all_species_names, names(abund_df))
  for (sp in missing_species) {
    abund_df[[sp]] <- 0
  }

  # Ensure all sites are present, even empty ones
  result_df <- data.frame(quadrat_id = quadrats$quadrat_id) %>%
    left_join(abund_df, by = "quadrat_id") %>%
    mutate(across(-quadrat_id, ~ replace_na(., 0))) %>%
    select(site = quadrat_id, all_of(all_species_names)) %>%
    arrange(site)

  return(result_df)
}


#' @title Calculate Mean Environmental Conditions per Quadrat
#' @description Extracts the mean value of environmental variables for each quadrat.
#' @param env_grid A data frame of gridded environmental data.
#' @param quadrats An `sf` object of quadrat polygons.
#' @param domain_crs The coordinate reference system of the domain.
#' @return A data frame with `site` (quadrat_id) and mean env values.
calculate_quadrat_environment <- function(env_grid, quadrats, domain_crs) {
  env_sf <- st_as_sf(env_grid, coords = c("x", "y"), crs = domain_crs)

  # Join grid points to the quadrats they fall within
  joined_data <- st_join(quadrats, env_sf)

  # Calculate mean environmental values for each quadrat
  site_env <- joined_data %>%
    st_drop_geometry() %>%
    group_by(site = quadrat_id) %>%
    summarise(
      across(
        c(temperature_C, elevation_m, rainfall_mm),
        ~ mean(., na.rm = TRUE)
      ),
      .groups = "drop"
    )

  # Ensure all sites are present
  full_site_env <- data.frame(site = quadrats$quadrat_id) %>%
    left_join(site_env, by = "site")

  return(full_site_env)
}


# -----------------------------------------------------------------------------
# 6. CORE PLOTTING & REPORTING FUNCTIONS
# -----------------------------------------------------------------------------

#' @title Create Master Plot of Spatial Simulation
#' @description Generates a ggplot object showing the domain, species distributions,
#'   quadrats, and optionally an environmental gradient overlay.
#' @param domain `sf` polygon of the study area.
#' @param species `sf` points of individuals.
#' @param quadrats `sf` polygons of sampling quadrats.
#' @param P Configuration list.
#' @param show_gradient Logical, whether to show an environmental gradient.
#' @param env_gradients Data frame of environmental grid data.
#' @param gradient_type The name of the gradient to display.
#' @return A `ggplot` object.
plot_spatial_sampling <- function(domain, species, quadrats, P,
                                  show_gradient = FALSE, env_gradients = NULL,
                                  gradient_type = "temperature") {
  # Define a color palette robust to varying numbers of species
  all_species_names <- LETTERS[1:P$N_SPECIES]
  species_colors <- rev(colorspace::sequential_hcl(P$N_SPECIES, palette = "RdPu"))
  names(species_colors) <- all_species_names

  # Basic plot structure
  p <- ggplot() +
    geom_sf(data = domain, fill = "grey70", color = P$FOREGROUND_COLOUR, linewidth = 0.5, alpha = 0.4)

  # Add environmental gradient layer if requested
  if (show_gradient && !is.null(env_gradients) && gradient_type %in% names(env_gradients)) {
    p <- p +
      metR::geom_contour_fill(data = env_gradients, aes(x = x, y = y, z = .data[[gradient_type]]), alpha = 0.5) +
      scale_fill_viridis(option = "viridis", name = str_to_title(gsub("_", " ", gradient_type)), guide = "colorbar")
  }

  # Add species points and quadrat layers
  p <- p +
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
      legend.position = "none" # Remove individual legends for panel plot
    ) +
    labs(
      title = if (show_gradient) str_replace_all(str_to_title(gradient_type), "_", " ") else "Species Distribution",
      subtitle = paste(P$N_INDIVIDUALS, "Individuals |", P$N_SPECIES, "Species |", nrow(quadrats), "Quadrats")
    )

  return(p)
}


#' @title Generate a Full, Formatted Analysis Report
#' @description Creates a string containing a full analysis of the simulation
#'   results, formatted exactly as requested.
#' @param res A list containing all simulation artifacts (data frames, stats, etc.).
#' @return A single character string with the formatted report.
generate_full_report <- function(res) {
  # --- convenience converters for units --------------------------------------
  .opt_to_units <- function(opt, gname) {
    # map normalized [0,1] to the reporting scales used elsewhere
    switch(gname,
      "temperature" = opt * 30 - 2, # °C  (-2 to 28)
      "elevation" = opt * 2000, # m   (0 to 2000)
      "rainfall" = opt * 700 + 200, # mm  (200 to 900)
      NA_real_
    )
  }

  .tol_to_units <- function(tol, gname) {
    # tolerance is a width on the normalized axis; convert to width in units
    switch(gname,
      "temperature" = tol * 30, # °C width
      "elevation" = tol * 2000, # m width
      "rainfall" = tol * 700, # mm width
      NA_real_
    )
  }

  .fmt_grad_line <- function(rows, gname, units_label) {
    if (nrow(rows) == 0) {
      return("None")
    }
    paste(
      sprintf(
        "%s (opt=%.2f [%.1f %s], tol=%.2f [%.1f %s])",
        rows$species,
        rows$optimum,
        .opt_to_units(rows$optimum, gname),
        units_label,
        rows$tol,
        .tol_to_units(rows$tol, gname),
        units_label
      ),
      collapse = ", "
    )
  }

  # --- Environmental Gradients ------------------------------------------------
  env_report <- c("\nEnvironmental Gradients:")
  temp_range <- range(res$env_gradients$temperature_C, na.rm = TRUE)
  elev_range <- range(res$env_gradients$elevation_m, na.rm = TRUE)
  rain_range <- range(res$env_gradients$rainfall_mm, na.rm = TRUE)

  # Prefer the tidy P$GRADIENT table; otherwise fall back to legacy vectors
  if (!is.null(res$P$GRADIENT)) {
    G <- res$P$GRADIENT
    # per-gradient pretty strings with both normalized and physical units
    temp_species_str <- .fmt_grad_line(G[G$gradient == "temperature", , drop = FALSE], "temperature", "°C")
    elev_species_str <- .fmt_grad_line(G[G$gradient == "elevation", , drop = FALSE], "elevation", "m")
    rain_species_str <- .fmt_grad_line(G[G$gradient == "rainfall", , drop = FALSE], "rainfall", "mm")
  } else {
    # --- Legacy fallback (keeps your old behaviour) ---------------------------
    temp_resp <- which(res$P$GRADIENT_ASSIGNMENTS == "temperature")
    elev_resp <- which(res$P$GRADIENT_ASSIGNMENTS == "elevation")
    rain_resp <- which(res$P$GRADIENT_ASSIGNMENTS == "rainfall")

    # Interpret 3-vector GRADIENT_OPTIMA as by-gradient; otherwise index by species position
    legacy_opt <- res$P$GRADIENT_OPTIMA
    opt_by <- if (length(legacy_opt) == 3) "gradient" else "species"

    get_opt <- function(idx, gradient_name) {
      if (opt_by == "gradient") {
        # temperature,elevation,rainfall assumed in that semantic space
        g_map <- c(temperature = 1L, elevation = 2L, rainfall = 3L)
        legacy_opt[g_map[[gradient_name]]]
      } else {
        legacy_opt[idx]
      }
    }

    # Strings as before (normalized → physical unit shown only for opt)
    temp_species_str <- if (length(temp_resp) > 0) {
      paste(sprintf(
        "%s (optimum %.1f °C)",
        res$P$GRADIENT_SPECIES[temp_resp],
        (get_opt(temp_resp, "temperature") * 30 - 2)
      ), collapse = ", ")
    } else {
      "None"
    }

    elev_species_str <- if (length(elev_resp) > 0) {
      paste(sprintf(
        "%s (optimum %.0f m)",
        res$P$GRADIENT_SPECIES[elev_resp],
        (get_opt(elev_resp, "elevation") * 2000)
      ), collapse = ", ")
    } else {
      "None"
    }

    rain_species_str <- if (length(rain_resp) > 0) {
      paste(sprintf(
        "%s (optimum %.0f mm)",
        res$P$GRADIENT_SPECIES[rain_resp],
        (get_opt(rain_resp, "rainfall") * 700 + 200)
      ), collapse = ", ")
    } else {
      "None"
    }
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

  # --- Gradient Correlations --------------------------------------------------
  cor_mat <- cor(res$env_gradients[, c("temperature_C", "elevation_m", "rainfall_mm")], use = "complete.obs")
  interp <- if (max(abs(cor_mat[upper.tri(cor_mat)])) < 0.3) "Gradients are approximately orthogonal (low correlation)" else "Some gradient correlation detected"
  corr_report <- c(
    "\nGradient Correlations:",
    sprintf("  Temperature-Elevation: r=%.3f", cor_mat[1, 2]),
    sprintf("  Temperature-Rainfall:  r=%.3f", cor_mat[1, 3]),
    sprintf("  Elevation-Rainfall:    r=%.3f", cor_mat[2, 3]),
    paste0("  Interpretation: ", interp)
  )

  # --- Species Abundance Distribution ----------------------------------------
  sad <- as.data.frame(table(res$species_dist$species)) %>%
    `colnames<-`(c("Species", "Count")) %>%
    arrange(desc(Count)) %>%
    mutate(
      Percent = 100 * Count / sum(Count),
      Role = case_when(
        !is.null(res$P$GRADIENT) & Species %in% res$P$GRADIENT$species ~ {
          g <- res$P$GRADIENT$gradient[match(Species, res$P$GRADIENT$species)]
          sprintf("[%s-RESPONSIVE]", toupper(g))
        },
        is.null(res$P$GRADIENT) & Species %in% res$P$GRADIENT_SPECIES ~ {
          g <- res$P$GRADIENT_ASSIGNMENTS[match(Species, res$P$GRADIENT_SPECIES)]
          sprintf("[%s-RESPONSIVE]", toupper(g))
        },
        Species == "A" ~ "[DOMINANT - clustered]",
        TRUE ~ "[SUBORDINATE]"
      )
    )
  sad_report <- c("\nSpecies Abundance Distribution:")
  sad_report <- c(
    sad_report,
    sprintf("  Species %s: %3d individuals (%5.1f%%) %s", sad$Species, sad$Count, sad$Percent, sad$Role)
  )

  # --- Spatial Distribution of Alpha Diversity --------------------------------
  alpha_report <- c("\nSpatial Distribution of Alpha Diversity:")
  for (i in 1:nrow(res$quadrats)) {
    spp_in_q <- suppressWarnings(st_intersection(res$species_dist, res$quadrats[i, ]))
    if (nrow(spp_in_q) > 0) {
      abunds <- as.data.frame(table(spp_in_q$species))
      abunds_str <- paste(sprintf("%s(%s)", abunds$Var1, abunds$Freq), collapse = ", ")
      alpha_report <- c(
        alpha_report,
        sprintf("  Quadrat %2d: α = %2d species | N = %3d individuals", i, n_distinct(spp_in_q$species), nrow(spp_in_q)),
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

  # --- Diversity Partitioning -------------------------------------------------
  abund_data <- res$abund_matrix %>% select(-site)
  richness_data <- abund_data %>% mutate(richness = rowSums(. > 0), n_ind = rowSums(.))
  shannon_H <- diversity(abund_data, index = "shannon")
  simpson_D <- diversity(abund_data, index = "simpson") # vegan's "simpson" is 1-D
  mean_alpha <- mean(richness_data$richness)
  se_alpha <- sd(richness_data$richness) / sqrt(nrow(richness_data))
  gamma_div <- n_distinct(res$species_dist$species)
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
    sprintf("Shannon's H' (mean): %.3f ± %.3f SE", mean(shannon_H), sd(shannon_H) / sqrt(length(shannon_H))),
    sprintf("Simpson's (1-D, mean): %.3f ± %.3f SE", mean(simpson_D), sd(simpson_D) / sqrt(length(simpson_D))),
    sprintf("Gamma (regional species pool): %d species", gamma_div),
    sprintf("Beta (Whittaker): %.2f", beta_whittaker),
    sprintf("Beta (additive): %.2f", beta_additive),
    sprintf("Mean pairwise beta (Sørensen): %.3f", mean_sorensen),
    sprintf("Mean quadrat abundance: %.1f ± %.1f", mean(richness_data$n_ind), sd(richness_data$n_ind)),
    sprintf("Abundance variation (CV): %.3f", sd(richness_data$n_ind) / mean(richness_data$n_ind))
  )

  # --- Spatial Autocorrelation & Fisher Validation ---------------------------
  space_dist <- dist(res$site_coords[, c("x", "y")])
  richness_dist <- dist(richness_data$richness)
  mantel_interp <- "  No significant spatial autocorrelation."
  if (nrow(res$site_coords) >= 4 && var(richness_data$richness) > 0) {
    mantel_test <- suppressWarnings(cor.test(space_dist, richness_dist, method = "pearson"))
    mantel_r <- mantel_test$estimate
    mantel_p <- mantel_test$p.value
    if (!is.na(mantel_p) && mantel_p < 0.05) {
      mantel_interp <- if (mantel_r > 0) {
        "  Significant positive autocorrelation: environmental filtering or dispersal limitation."
      } else {
        "  Significant negative autocorrelation: competitive exclusion or strong heterogeneity."
      }
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

  # --- Fisher's Log-Series Validation ----------------------------------------
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
  r_squared <- if (var(obs_abund, na.rm = TRUE) > 0) cor(obs_abund, theory_abund, use = "complete.obs")^2 else 1
  eff_alpha <- fisher.alpha(table(res$species_dist$species))

  fisher_report <- c(
    "\nFisher's Log Series Model Validation:",
    sprintf("  RMSE: %.2f", rmse),
    sprintf("  R-squared: %.3f", r_squared),
    sprintf("  Max residual: %.1f", max(abs(residuals), na.rm = TRUE)),
    sprintf("  Specified alpha: %.2f", res$P$FISHER_ALPHA),
    sprintf("  Effective alpha from data: %.2f", eff_alpha)
  )

  # --- Combine & return ------------------------------------------------------
  full_report <- c(
    "========== ANALYSIS REPORT ==========",
    env_report, corr_report, sad_report, alpha_report, div_report, spat_report, fisher_report,
    "\nSIMULATION COMPLETED SUCCESSFULLY."
  )

  paste(full_report, collapse = "\n")
}


# -----------------------------------------------------------------------------
# 7. MAIN SIMULATION ORCHESTRATOR
# -----------------------------------------------------------------------------

#' @title Run the Full Spatial Sampling Simulation
#' @description The main function that orchestrates the entire simulation.
#' @param init_file Path to the `.txt` configuration file.
#' @param output_prefix A string for the base of all output filenames.
run_spatial_simulation <- function(init_file = "simul_init.txt",
                                   interactions_file = NULL,
                                   output_prefix = "simulation_output") {
  # --- 1. Setup (paths, report sink) ---
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

  # --- 2. Core simulation block ---
  results_list <- NULL
  tryCatch(
    {
      P <- load_config(init_file)

      # Load interactions from separate file
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
      site_env <- calculate_quadrat_environment(env_gradients, quadrats, st_crs(domain))
      site_coords <- suppressWarnings(st_coordinates(st_centroid(quadrats))) |>
        as.data.frame() |>
        dplyr::mutate(site = quadrats$quadrat_id) |>
        dplyr::select(site, x = X, y = Y)

      # --- 3. Save outputs ---
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
        P = P,
        domain = domain,
        species_dist = species_dist,
        quadrats = quadrats,
        env_gradients = env_gradients,
        abund_matrix = abund_matrix,
        site_coords = site_coords
      )
    },
    error = function(e) {
      cat("\nERROR in core simulation:\n", e$message, "\n", sep = "")
    }
  )

  # --- 4. Advanced analysis (optional) ---
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

  # --- 5. Final reporting ---
  if (!is.null(results_list)) {
    cat("\n", generate_full_report(results_list), "\n", sep = "")
    cat("Outputs saved to: ", normalizePath(output_dir), "\n", sep = "")
  }

  invisible(NULL)
}
