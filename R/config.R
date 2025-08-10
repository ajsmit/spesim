#' Parse a parameter file into raw key-value pairs (internal)
#'
#' @description
#' Lightweight reader for \code{KEY = value} lines. It preserves quoted
#' strings (e.g., hex colours), supports comma-separated vectors, simple
#' \code{c(...)} vectors, and \code{name:value} entries. No defaults or
#' validation are applied here - this is a low-level helper used by
#' \code{\link{load_config}()}.
#'
#' @param filename Character scalar; path to the configuration file.
#'
#' @return A named list of raw values (characters, numerics, logicals, or
#'   atomic vectors) suitable for merging into defaults. Names are kept
#'   exactly as given in the file (case preserved); the caller may
#'   normalise/uppercase if desired.
#'
#' @keywords internal
#' @noRd
.parse_init_file <- function(filename) {
  if (!file.exists(filename)) stop("Parameter file not found: ", filename)
  raw <- readLines(filename, warn = FALSE)

  # Keep hex colours like "#22223b" by only removing comments that start
  # *after* some whitespace and a '#'. We do not remove inline '#' inside quotes.
  strip_comments <- function(x) sub("\\s+#.*$", "", x)
  lines <- trimws(vapply(raw, strip_comments, character(1)))
  lines <- lines[nzchar(lines)]

  # Gather key=value lines, allowing multi-line c(...)
  kvs <- list()
  i <- 1L
  while (i <= length(lines)) {
    ln <- lines[i]
    if (!grepl("=", ln, fixed = TRUE)) {
      i <- i + 1L
      next
    }
    parts <- strsplit(ln, "=", fixed = TRUE)[[1]]
    key <- trimws(parts[1])
    val <- trimws(paste(parts[-1], collapse = "="))

    if (grepl("^c\\s*\\(", val)) {
      # accumulate until balanced parentheses
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
    # strip matching outer quotes
    v <- gsub('^"(.*)"$', "\\1", v)
    v <- gsub("^'(.*)'$", "\\1", v)

    # c(...) literal
    if (grepl("^c\\s*\\(", v)) {
      out <- tryCatch(eval(parse(text = v)), error = function(e) e)
      if (inherits(out, "error")) {
        stop("Could not parse vector literal: ", v, "\n", out$message)
      }
      return(out)
    }

    # comma-separated items
    if (grepl(",", v, fixed = TRUE)) {
      items <- trimws(strsplit(v, ",", fixed = TRUE)[[1]])
      nums <- suppressWarnings(as.numeric(items))
      return(if (!anyNA(nums)) nums else items)
    }

    # logical
    if (v %in% c("TRUE", "FALSE")) {
      return(as.logical(v))
    }

    # integer
    if (grepl("^[+-]?[0-9]+$", v)) {
      return(as.integer(v))
    }

    # numeric (incl. scientific)
    if (grepl("^[+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?$", v)) {
      return(as.numeric(v))
    }

    # name:value single entry
    if (grepl("^[^:]+:.+$", v)) {
      kv <- strsplit(v, ":", fixed = TRUE)[[1]]
      nm <- trimws(kv[1])
      vv <- trimws(kv[2])
      num <- suppressWarnings(as.numeric(vv))
      out <- if (!is.na(num)) c(num) else c(vv)
      names(out) <- nm
      return(out)
    }

    # fallback: raw string
    v
  }

  out <- list()
  for (kv in kvs) out[[kv$key]] <- parse_value(kv$val)
  out
}


#' Load Simulation Configuration (with defaults & validation)
#'
#' @description
#' Reads a text configuration file (KEY = value), merges values with package
#' defaults, validates options, resolves gradient specifications into a tidy
#' table, materialises quadrat sizes, and sets the random seed. This is the
#' **only** supported public entry point for reading configs; low-level parsing
#' is handled internally.
#'
#' @details
#' Supported value forms:
#' \itemize{
#'   \item Scalars: \code{300}, \code{TRUE}, \code{0.15}, \code{"#22223b"}.
#'   \item Comma vectors: \code{A,B,C} or \code{0.1, 0.2, 0.3}.
#'   \item R-style vectors: \code{c(A, B, C)} or \code{c(0.1, 0.2)} (multi-line OK).
#'   \item Named entries: \code{A:0.55}, \code{temperature:0.12}.
#' }
#' Gradient keys supported: \code{temperature}, \code{elevation}, \code{rainfall}.
#' Use \code{GRADIENT_SPECIES}, \code{GRADIENT_ASSIGNMENTS}, plus either
#' \code{GRADIENT_OPTIMA} and \code{GRADIENT_TOLERANCE} as scalar, per-species
#' named, or per-gradient named values.
#'
#' @param init_file Character scalar; path to the configuration file.
#'
#' @return A named list \code{P} with fully-resolved parameters used by the
#' simulator (including a tidy \code{P$GRADIENT} tibble and \code{P$QUADRAT_SIZE}).
#'
#' @seealso \code{\link{run_spatial_simulation}}
#' @examples
#' \dontrun{
#' P <- load_config("inst/examples/spesim_init_complete.txt")
#' }
#' @export
load_config <- function(init_file) {
  message("========== INITIALISING SPATIAL SAMPLING SIMULATION ==========")

  # Defaults
  P <- list(
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
    ADVANCED_ANALYSIS = FALSE,

    # Spatial point-process (advanced)
    SPATIAL_PROCESS_A = "poisson", # "poisson" | "thomas"
    A_PARENT_INTENSITY = NA_real_,
    A_MEAN_OFFSPRING = 10,
    A_CLUSTER_SCALE = 1,
    SPATIAL_PROCESS_OTHERS = "poisson", # "poisson" | "strauss" | "geyer"
    OTHERS_BETA = NA_real_,
    OTHERS_GAMMA = NA_real_,
    OTHERS_R = 1,
    OTHERS_S = 2
  )

  # Parse file (internal helper)
  raw_vals <- .parse_init_file(init_file)

  # Merge, allowing keys in any case in the file (normalise to UPPER)
  if (length(raw_vals)) {
    # map upper names to raw names
    raw_names_upper <- toupper(names(raw_vals))
    dup <- duplicated(raw_names_upper)
    if (any(dup)) {
      stop(
        "Duplicate keys (case-insensitive) in config: ",
        paste(unique(raw_names_upper[dup]), collapse = ", ")
      )
    }
    names(raw_vals) <- raw_names_upper
    for (nm in intersect(names(P), names(raw_vals))) {
      P[[nm]] <- raw_vals[[nm]]
    }
  }

  # Normalise some strings
  P$SPATIAL_PROCESS_A <- tolower(as.character(P$SPATIAL_PROCESS_A))
  P$SPATIAL_PROCESS_OTHERS <- tolower(as.character(P$SPATIAL_PROCESS_OTHERS))

  if (!P$SPATIAL_PROCESS_A %in% c("poisson", "thomas")) {
    stop("SPATIAL_PROCESS_A must be 'poisson' or 'thomas'.")
  }
  if (!P$SPATIAL_PROCESS_OTHERS %in% c("poisson", "strauss", "geyer")) {
    stop("SPATIAL_PROCESS_OTHERS must be 'poisson', 'strauss', or 'geyer'.")
  }

  # Safe numeric coercion for known numeric fields (no warnings)
  num_fields <- c(
    "SEED", "N_INDIVIDUALS", "N_SPECIES",
    "DOMINANT_FRACTION", "FISHER_ALPHA", "FISHER_X",
    "SAMPLING_RESOLUTION", "ENVIRONMENTAL_NOISE",
    "MAX_CLUSTERS_DOMINANT", "CLUSTER_SPREAD_DOMINANT",
    "INTERACTION_RADIUS", "N_QUADRATS",
    "N_TRANSECTS", "N_QUADRATS_PER_TRANSECT", "TRANSECT_ANGLE",
    "VORONOI_SEED_FACTOR",
    "POINT_SIZE", "POINT_ALPHA", "QUADRAT_ALPHA",
    "A_PARENT_INTENSITY", "A_MEAN_OFFSPRING", "A_CLUSTER_SCALE",
    "OTHERS_BETA", "OTHERS_GAMMA", "OTHERS_R", "OTHERS_S"
  )
  for (f in intersect(names(P), num_fields)) {
    # only coerce if not already numeric (avoid warnings)
    if (!is.numeric(P[[f]]) && !is.logical(P[[f]]) && !is.null(P[[f]])) {
      suppressWarnings(P[[f]] <- as.numeric(P[[f]]))
    }
  }

  # Logical flags
  P$ADVANCED_ANALYSIS <- as.logical(P$ADVANCED_ANALYSIS)

  # Colours (validate gracefully)
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
  P$BACKGROUND_COLOUR <- .validate_colour(P$BACKGROUND_COLOUR, "white")
  P$FOREGROUND_COLOUR <- .validate_colour(P$FOREGROUND_COLOUR, "#22223b")
  P$QUADRAT_COLOUR <- .validate_colour(P$QUADRAT_COLOUR, "black")

  # --- Gradients -------------------------------------------------------
  gs <- as.character(P$GRADIENT_SPECIES)
  ga <- as.character(P$GRADIENT_ASSIGNMENTS)

  if (length(gs) != length(ga)) {
    stop("GRADIENT_SPECIES and GRADIENT_ASSIGNMENTS must have the same length.")
  }
  if (length(gs) && !all(ga %in% c("temperature", "elevation", "rainfall"))) {
    stop("Unknown gradient(s): ", paste(setdiff(ga, c("temperature", "elevation", "rainfall")), collapse = ", "))
  }

  # Utilities expected elsewhere in the package:
  # - .parse_named_pairs_numeric(x)  (returns numeric named vector or scalar)
  # - `%||%`                         (lhs if not NULL/length>0 else rhs)
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
  qs <- QUADRAT_SIZES[[as.character(P$QUADRAT_SIZE_OPTION)]]
  if (is.null(qs)) stop("Invalid QUADRAT_SIZE_OPTION. Use small|medium|large.")
  P$QUADRAT_SIZE <- qs

  # Reproducibility
  set.seed(as.integer(P$SEED))
  P
}


#' Load interspecific interaction settings (radius + matrix)
#'
#' @description
#' Reads a separate interactions config file and returns a local-interaction
#' radius along with an \eqn{S \times S} interaction matrix for \eqn{S}
#' species (rows = focal, columns = neighbour). If the file is missing or
#' incomplete, a neutral matrix of ones and radius 0 are returned.
#'
#' @details
#' The interactions file is parsed as simple \code{KEY = value} pairs. The
#' following keys are recognized:
#' \itemize{
#'   \item \code{INTERACTION_RADIUS}: non-negative numeric; distance threshold
#'         within which neighbours influence assignment probabilities.
#'   \item \code{MATRIX_CSV}: path to a CSV containing a (sub)matrix with row
#'         names and column names corresponding to species labels (e.g. A..Z).
#'         Any missing rows/columns are filled with \code{1}.
#'   \item \code{EDGELIST_CSV}: path to a CSV with columns
#'         \code{focal, neighbor, value}; entries are inserted into the
#'         appropriate cells of the full matrix.
#'   \item \code{AUTO}: optional flag (\code{TRUE}/\code{FALSE}) enabling a
#'         simple built-in example structure if no CSV is provided.
#' }
#' If both \code{MATRIX_CSV} and \code{EDGELIST_CSV} are supplied,
#' \code{MATRIX_CSV} takes precedence.
#'
#' @param interactions_file Character scalar or \code{NULL}. Path to the
#'   interactions configuration file. If \code{NULL} or not found, the
#'   returned radius is 0 and all coefficients are \code{1}.
#' @param n_species Integer. Number of species; defines the size of the
#'   returned matrix and the expected label set \code{LETTERS[1:n_species]}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{radius}}{Numeric length-1; the interaction radius (0 disables interactions).}
#'   \item{\code{matrix}}{Numeric matrix of size \eqn{S \times S} with dimnames
#'         \code{LETTERS[1:S]}. Non-finite values are rejected.}
#' }
#'
#' @examples
#' \dontrun{
#' I <- load_interactions("config/interactions.txt", n_species = 10)
#' I$radius
#' I$matrix[1:3, 1:3]
#' }
#'
#' @seealso \code{\link{load_config}()}
#' @export
load_interactions <- function(interactions_file, n_species) {
  spp_names <- LETTERS[1:n_species]
  out_radius <- 0
  IM <- matrix(1.0, nrow = n_species, ncol = n_species, dimnames = list(spp_names, spp_names))

  if (is.null(interactions_file) || !file.exists(interactions_file)) {
    if (!is.null(interactions_file) && !file.exists(interactions_file)) {
      warning("interactions_file not found: ", interactions_file, " - interactions disabled.")
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
