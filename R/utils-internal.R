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
    stop("Row-delimited INTERACTION_MATRIX must be S rows of S numeric values (SxS).")
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

#' Detect availability of the fast Thomas engine
#'
#' @keywords internal
.has_cpp_thomas <- function() {
  f <- get0("rthomas_bbox_cpp", envir = asNamespace("spesim"), inherits = FALSE)
  isTRUE(is.function(f))
}

#' Detect availability of the fast Strauss engine
#'
#' @keywords internal
.has_cpp_strauss <- function() {
  f <- get0("rstrauss_bbox_cpp", envir = asNamespace("spesim"), inherits = FALSE)
  isTRUE(is.function(f))
}

#' Detect availability of the fast Geyer engine
#'
#' @keywords internal
.has_cpp_geyer <- function() {
  f <- get0("rgeyer_bbox_cpp", envir = asNamespace("spesim"), inherits = FALSE)
  isTRUE(is.function(f))
}
