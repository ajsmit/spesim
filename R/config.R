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
#'   \item R-style character vectors for inline edgelists: e.g.
#'      \code{INTERACTIONS_EDGELIST = c("A,B-D,0.8", "C,A,1.2", "E,*,0.95")}.
#'   \item Named entries: \code{A:0.55}, \code{temperature:0.12}.
#' }
#' You may specify INTERACTIONS_FILE = "path/to/interactions.txt" to point to
#' a separate interactions config, or embed rules directly with
#' \code{INTERACTIONS_EDGELIST = c("A,B-D,0.8", "C,A,1.2")}. If neither is provided,
#' interactions default to neutral (all 1.0) with \code{INTERACTION_RADIUS = 0}.
#'
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
    INTERACTIONS_FILE = NULL,
    INTERACTIONS_EDGELIST = NULL,
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

  # Keep interaction pointers / edgelists as character
  if (!is.null(P$INTERACTIONS_EDGELIST)) {
    P$INTERACTIONS_EDGELIST <- as.character(P$INTERACTIONS_EDGELIST)
  }
  if (!is.null(P$INTERACTIONS_FILE)) {
    P$INTERACTIONS_FILE <- as.character(P$INTERACTIONS_FILE)
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


# --- helpers ---------------------------------------------------------------
.parse_neighbor_spec <- function(neighbor_field, focal, spp_names) {
  # neighbor_field: e.g. "B-D,L" (cell is already parsed as one string by read.csv)
  # Supports commas/semicolons/pipes as separators inside the cell.
  if (is.na(neighbor_field) || !nzchar(neighbor_field)) {
    return(character(0))
  }
  s <- gsub("\\s+", "", neighbor_field)
  if (s == "*") {
    return(setdiff(spp_names, focal)) # everyone except focal
  }
  parts <- unlist(strsplit(s, "[,;|]+"))
  out <- character(0)
  for (p in parts) {
    if (grepl("^[A-Z]-[A-Z]$", p)) {
      a <- substr(p, 1, 1)
      b <- substr(p, 3, 3)
      rng <- LETTERS[which(LETTERS == a):which(LETTERS == b)]
      out <- c(out, rng)
    } else {
      out <- c(out, p)
    }
  }
  unique(out)
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
  radius <- 0
  IM <- matrix(1.0, n_species, n_species, dimnames = list(spp_names, spp_names))

  if (is.null(interactions_file) || !file.exists(interactions_file)) {
    warning("interactions_file not found; interactions disabled.")
    return(list(radius = radius, matrix = IM))
  }

  ext <- tolower(tools::file_ext(interactions_file))

  # -------------------- NEW CSV FORMAT (Option 1) --------------------
  if (ext %in% c("csv")) {
    df <- tryCatch(
      utils::read.csv(interactions_file, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) stop("Failed to read interactions CSV: ", e$message)
    )

    needed <- c("focal", "neighbor", "value")
    if (!all(needed %in% names(df))) {
      stop("Interactions CSV must have columns: focal, neighbor, value")
    }

    # Optional metadata row for radius:
    # _RADIUS_, "", 50
    meta_rows <- which(toupper(df$focal) == "_RADIUS_")
    if (length(meta_rows) > 0) {
      val <- suppressWarnings(as.numeric(df$value[meta_rows[1]]))
      if (is.finite(val) && val >= 0) radius <- val
      df <- df[-meta_rows, , drop = FALSE]
    }

    if (nrow(df) == 0) {
      return(list(radius = radius, matrix = IM))
    }

    for (i in seq_len(nrow(df))) {
      f <- as.character(df$focal[i])
      nb_raw <- as.character(df$neighbor[i])
      v <- suppressWarnings(as.numeric(df$value[i]))

      if (!nzchar(f) || !is.finite(v)) next
      if (!(f %in% spp_names)) {
        warning("Ignoring row ", i, ": focal '", f, "' not in species A..", LETTERS[n_species])
        next
      }

      targets <- .parse_neighbor_spec(nb_raw, focal = f, spp_names = spp_names)
      targets <- intersect(targets, spp_names)
      if (length(targets) == 0) next

      IM[f, targets] <- v
    }

    if (any(!is.finite(IM))) stop("INTERACTION_MATRIX contains non-finite values.")
    return(list(radius = radius, matrix = IM))
  }

  # -------------------- BACKWARD-COMPATIBLE .txt (key=value) -----------------
  raw <- readLines(interactions_file, warn = FALSE)
  raw <- trimws(sub("#.*$", "", raw))
  raw <- raw[nzchar(raw)]
  if (length(raw) == 0) {
    return(list(radius = radius, matrix = IM))
  }

  kv <- strsplit(raw, "=", fixed = TRUE)
  K <- toupper(trimws(vapply(kv, `[`, "", 1)))
  V <- trimws(vapply(kv, function(x) paste(x[-1], collapse = "="), "")) # allow '=' in RHS
  params <- stats::setNames(V, K)

  # radius
  if ("INTERACTION_RADIUS" %in% names(params)) {
    rnum <- suppressWarnings(as.numeric(params[["INTERACTION_RADIUS"]]))
    if (length(rnum) == 1 && is.finite(rnum) && rnum >= 0) radius <- rnum
  }

  # CSV pointers (legacy)
  if ("MATRIX_CSV" %in% names(params) && nzchar(params[["MATRIX_CSV"]])) {
    p <- params[["MATRIX_CSV"]]
    M <- as.matrix(utils::read.csv(p, row.names = 1, check.names = FALSE))
    fullM <- matrix(1.0, n_species, n_species, dimnames = list(spp_names, spp_names))
    rn <- intersect(rownames(M), spp_names)
    cn <- intersect(colnames(M), spp_names)
    fullM[rn, cn] <- as.numeric(M[rn, cn, drop = FALSE])
    IM <- fullM
  } else if ("EDGELIST_CSV" %in% names(params) && nzchar(params[["EDGELIST_CSV"]])) {
    p <- params[["EDGELIST_CSV"]]
    E <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
    need_cols <- c("focal", "neighbor", "value")
    if (!all(need_cols %in% names(E))) stop("EDGELIST_CSV must have columns: focal, neighbor, value")
    for (i in seq_len(nrow(E))) {
      f <- as.character(E$focal[i])
      n <- as.character(E$neighbor[i])
      v <- suppressWarnings(as.numeric(E$value[i]))
      if (f %in% spp_names && n %in% spp_names && is.finite(v)) IM[f, n] <- v
    }
  } else if (tolower(params[["AUTO"]] %||% "false") %in% c("true", "1", "yes")) {
    if (all(c("D", "E") %in% spp_names)) IM["E", "D"] <- 3.0
  }

  if (any(!is.finite(IM))) stop("INTERACTION_MATRIX contains non-finite values.")
  list(radius = radius, matrix = IM)
}


#' Parse inline interaction rules from the init file
#'
#' Converts `INTERACTIONS_EDGELIST` (CSV-like string; optional `INTERACTION_RADIUS`)
#' into the resolved interaction list (`radius`, `matrix`).
#'
#' @param rules Character scalar or character vector of CSV-like lines with columns
#'   `focal,neighbour,value` and optional wildcard/range syntax (see vignette).
#' @param radius Optional numeric radius override; if `NULL`, defaults to 0.
#' @param n_species Integer S; used to define the full SxS matrix.
#' @return A list with `radius` (numeric) and `matrix` (SxS numeric with dimnames).
#' @examples
#' \dontrun{
#' load_interactions_inline(
#'   rules = c("A,B-D,0.8", "C,A,1.2", "E,*,0.95"),
#'   radius = 2, n_species = 10
#' )
#' }
#' @export
load_interactions_inline <- function(rules, n_species, radius = 0) {
  spp <- LETTERS[1:n_species]
  M <- matrix(1, n_species, n_species, dimnames = list(spp, spp))
  r <- as.numeric(radius %||% 0)

  # permit an optional "radius=..." line up front
  head_like <- trimws(rules)
  if (length(head_like) && grepl("^radius\\s*=", head_like[1], ignore.case = TRUE)) {
    r <- as.numeric(sub("^radius\\s*=\\s*", "", head_like[1], ignore.case = TRUE))
    rules <- rules[-1]
  }

  # reuse the same expand/shorthand you already implemented for CSV rows
  expand_targets <- function(neighbor) {
    neighbor <- trimws(neighbor)
    if (neighbor == "*") {
      return(setdiff(spp, focal))
    } # filled per-row below
    # A,B-D,L  -> expand ranges and singles
    items <- unlist(strsplit(neighbor, "\\s*,\\s*"))
    out <- character(0)
    for (it in items) {
      if (grepl("-", it)) {
        ends <- unlist(strsplit(it, "\\s*-\\s*"))
        a <- match(ends[1], LETTERS)
        b <- match(ends[2], LETTERS)
        if (is.na(a) || is.na(b)) next
        out <- c(out, LETTERS[seq(min(a, b), max(a, b))])
      } else {
        out <- c(out, it)
      }
    }
    unique(out)
  }

  for (ln in rules) {
    if (!nzchar(trimws(ln))) next
    # expected: focal,neighbors,value   (commas separate the 3 fields only)
    parts <- unlist(strsplit(ln, "\\s*,\\s*"))
    if (length(parts) != 3) stop("Bad inline interaction row: '", ln, "'. Expect: focal,neighbors,value")
    focal <- parts[1]
    neighbors_raw <- parts[2]
    value <- suppressWarnings(as.numeric(parts[3]))
    if (!(focal %in% spp) || !is.finite(value)) stop("Invalid focal/value in row: '", ln, "'")
    if (neighbors_raw == "*") {
      targets <- setdiff(spp, focal)
    } else {
      targets <- expand_targets(neighbors_raw)
    }
    targets <- intersect(targets, spp)
    if (length(targets)) M[focal, targets] <- value
  }

  list(radius = r, matrix = M)
}


#' Validate an interaction specification
#'
#' Checks that the resolved interaction list (radius + matrix) is well formed
#' for the given species set. Warns or errors on shape/NA/finite/name issues.
#'
#' @param I A list with elements `radius` (numeric scalar) and `matrix` (S×S numeric).
#' @param spp_names Character vector of expected species names
#'   (e.g., \verb{LETTERS[1:S]}).
#' @param stop_on_error Logical; stop on validation failure (`TRUE`) or just warn (`FALSE`).
#' @return Invisibly returns `TRUE` on success.
#' @examples
#' \dontrun{
#' I <- list(radius = 2, matrix = diag(10))
#' validate_interactions(I, LETTERS[1:10])
#' }
#' @export
validate_interactions <- function(I, spp_names, stop_on_error = FALSE) {
  msgs <- character()
  ok <- TRUE

  # structure
  if (!is.list(I) || !all(c("radius", "matrix") %in% names(I))) {
    msgs <- c(msgs, "Object must be a list with elements 'radius' and 'matrix'.")
    ok <- FALSE
  } else {
    r <- I$radius
    M <- I$matrix

    # radius checks
    if (!is.numeric(r) || length(r) != 1 || !is.finite(r) || r < 0) {
      msgs <- c(msgs, "INTERACTION_RADIUS must be a single non-negative, finite number.")
      ok <- FALSE
    }

    # matrix shape
    S <- length(spp_names)
    if (!is.matrix(M) || any(dim(M) != S)) {
      msgs <- c(msgs, sprintf(
        "Matrix must be %dx%d (got %sx%s).", S, S,
        if (is.matrix(M)) paste(dim(M), collapse = "x") else "NA",
        if (is.matrix(M)) paste(dim(M), collapse = "x") else "NA"
      ))
      ok <- FALSE
    } else {
      # names
      rn <- rownames(M)
      cn <- colnames(M)
      if (is.null(rn) || is.null(cn) || !all(rn == spp_names) || !all(cn == spp_names)) {
        msgs <- c(msgs, "Matrix row/col names must exactly match 'spp_names' in order.")
        ok <- FALSE
      }

      # values
      if (any(!is.finite(M))) {
        msgs <- c(msgs, "Matrix contains non-finite values.")
        ok <- FALSE
      }

      # quick sanity notes (not fatal)
      if (any(diag(M) != 1)) {
        msgs <- c(msgs, "Note: diagonal has values != 1 (self-effects present).")
      }
      if (any(M <= 0)) {
        msgs <- c(msgs, "Warning: matrix contains non-positive values (could zero out probabilities).")
      }
      # species with no listed effects either way (all 1s in row & col)
      barren <- spp_names[(rowSums(M != 1) == 0) & (colSums(M != 1) == 0)]
      if (length(barren)) {
        msgs <- c(msgs, sprintf("Note: no non-1 interactions for species: %s.", paste(barren, collapse = ", ")))
      }
    }
  }

  if (!ok && stop_on_error) stop(paste(msgs, collapse = "\n"))
  if (length(msgs)) message(paste(msgs, collapse = "\n"))
  invisible(list(ok = ok, messages = msgs))
}


#' Pretty-print a compact interaction matrix summary
#'
#' Prints the radius and a sorted list of non-1.0 coefficients from the
#' interaction matrix.
#'
#' @param I A list with `radius` and `matrix`.
#' @param digits Integer; number of digits to print.
#' @param top_n Integer; show at most this many non-1.0 entries.
#' @return Invisibly returns `NULL`.
#' @examples
#' \dontrun{
#' print_interactions(I, digits = 3, top_n = 20)
#' }
#' @export
print_interactions <- function(I, digits = 3, top_n = NULL) {
  stopifnot(is.list(I), "radius" %in% names(I), "matrix" %in% names(I))
  r <- I$radius
  M <- I$matrix
  spp <- rownames(M)
  S <- nrow(M)

  # basics
  cat("---- Interactions Summary ----\n")
  cat(sprintf("Species: %d (%s..%s)\n", S, spp[1], spp[S]))
  cat(sprintf("Radius : %s\n", format(r, digits = digits)))

  # sparsity & symmetry
  nz <- which(abs(M - 1) > .Machine$double.eps, arr.ind = TRUE)
  K <- nrow(nz)
  cat(sprintf("Non-1 entries: %d (%.1f%% of %d)\n", K, 100 * K / (S * S), S * S))
  if (K > 0) {
    symm_err <- mean(abs(M - t(M)))
    cat(sprintf("Asymmetry (mean |M - t(M)|): %.*f\n", digits, symm_err))
  }

  # small matrix view
  if (S <= 15) {
    cat("\nMatrix ('.' = 1):\n")
    prettyM <- apply(M, 2, function(col) ifelse(abs(col - 1) < 1e-12, ".", format(round(col, digits), nsmall = 0)))
    dimnames(prettyM) <- dimnames(M)
    print(noquote(prettyM))
  }

  # edgelist of deviations from 1
  if (K > 0) {
    df <- data.frame(
      focal = spp[nz[, 1]],
      neighbor = spp[nz[, 2]],
      value = as.numeric(M[nz]),
      delta = abs(as.numeric(M[nz]) - 1),
      stringsAsFactors = FALSE
    )
    df <- df[order(-df$delta, df$focal, df$neighbor), ]
    if (!is.null(top_n)) df <- utils::head(df, top_n)

    cat("\nNon-1 entries (sorted by |value-1|):\n")
    print(
      utils::head(transform(df[, c("focal", "neighbor", "value")],
        value = round(value, digits)
      ), if (is.null(top_n)) nrow(df) else top_n),
      row.names = FALSE
    )
  }
  cat("------------------------------\n")
}
