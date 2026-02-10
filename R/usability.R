#' Validate a spesim configuration (without running the simulation)
#'
#' @description
#' Checks a configuration file or parameter list for common problems and returns
#' a structured report (errors + warnings) without stopping.
#'
#' This is aimed at teaching workflows where students iterate on init files and
#' need fast feedback before running a full simulation.
#'
#' @param config Either a path to an init file (character scalar) or an in-memory
#'   parameter list `P`.
#' @param interactions_file Optional interactions file path; if provided, it is
#'   validated as well.
#' @param strict Logical; if `TRUE`, treat warnings as errors.
#' @return A list with components `ok` (logical), `errors` (character),
#'   `warnings` (character), and (when possible) `P` (resolved config).
#' @export
validate_config <- function(config, interactions_file = NULL, strict = FALSE) {
  errors <- character()
  warns <- character()
  P <- NULL

  tryCatch(
    {
      if (is.character(config) && length(config) == 1L) {
        P <- load_config(config)
      } else if (is.list(config)) {
        P <- config
      } else {
        stop("`config` must be a path to an init file or an in-memory parameter list.")
      }

      # minimal required keys
      req <- c("N_INDIVIDUALS", "N_SPECIES", "SAMPLING_SCHEME")
      missing <- setdiff(req, names(P))
      if (length(missing)) {
        errors <- c(errors, paste0("Missing required parameters: ", paste(missing, collapse = ", ")))
      }

      # scheme sanity
      if (!is.null(P$SAMPLING_SCHEME)) {
        scheme <- tolower(as.character(P$SAMPLING_SCHEME))
        ok_schemes <- c("random", "tiled", "systematic", "transect", "voronoi")
        if (!scheme %in% ok_schemes) {
          errors <- c(errors, sprintf("SAMPLING_SCHEME must be one of: %s", paste(ok_schemes, collapse = ", ")))
        }
      }

      # interactions validation (best-effort)
      if (!is.null(interactions_file)) {
        I <- load_interactions(interactions_file, n_species = P$N_SPECIES)
        v <- validate_interactions(I, LETTERS[1:P$N_SPECIES], stop_on_error = FALSE)
        if (!isTRUE(v$ok)) {
          errors <- c(errors, paste0("Interactions spec invalid:\n", paste(v$messages, collapse = "\n")))
        }
      }
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    error = function(e) {
      errors <<- c(errors, conditionMessage(e))
      NULL
    }
  )

  if (isTRUE(strict) && length(warns)) {
    errors <- c(errors, paste0("(strict=TRUE) ", warns))
    warns <- character()
  }

  ok <- length(errors) == 0 && (length(warns) == 0 || !isTRUE(strict))
  list(ok = ok, errors = unique(errors), warnings = unique(warns), P = P)
}


#' List shipped example configuration files
#'
#' @description
#' Prints (and returns) the example init files shipped in `inst/examples/`.
#'
#' @return A tibble with columns `name` and `path`.
#' @export
list_examples <- function() {
  ex_dir <- system.file("examples", package = "spesim")
  if (!nzchar(ex_dir)) {
    stop("Could not locate installed examples directory for spesim.")
  }
  files <- list.files(ex_dir, pattern = "\\.(txt|csv)$", full.names = TRUE)
  out <- tibble::tibble(
    name = basename(files),
    path = files
  )
  print(out)
  invisible(out)
}


#' Copy a shipped example init file into your working directory
#'
#' @description
#' Convenience helper for workshops: copy an example init file from the package
#' into the current directory (or a specified location) so students can edit it.
#'
#' @param example Basename of an example file (e.g. `"spesim_init_basic.txt"`).
#' @param to Destination path. Default is `example` in the current working directory.
#' @param overwrite Logical; overwrite an existing file?
#' @return The destination path (invisibly).
#' @export
spesim_use_example <- function(example = "spesim_init_basic.txt", to = example, overwrite = FALSE) {
  src <- system.file("examples", example, package = "spesim")
  if (!nzchar(src)) {
    stop("Example not found in package: ", example, "\nRun `list_examples()` to see what is available.")
  }
  ok <- file.copy(src, to, overwrite = overwrite)
  if (!ok) stop("Failed to copy example to: ", to)
  message("Copied example init file to: ", normalizePath(to, winslash = "/", mustWork = FALSE))
  invisible(to)
}


#' Open spesim documentation pages in your browser
#'
#' @description
#' Opens pkgdown documentation in a web browser. This is a small usability
#' helper for teaching settings.
#'
#' @param page One of: `"home"`, `"workflow"`, `"quadrat-placement"`,
#'   `"public-api"`, or `"recipes"`.
#' @param url_base Optional base URL. Defaults to the package URL in DESCRIPTION.
#' @export
spesim_open_docs <- function(page = c("home", "workflow", "quadrat-placement", "public-api", "recipes"),
                             url_base = "https://ajsmit.github.io/spesim/") {
  page <- match.arg(page)
  path <- switch(page,
    home = "index.html",
    workflow = "articles/spesim-workflow.html",
    `quadrat-placement` = "articles/spesim-quadrat-placement.html",
    `public-api` = "articles/spesim-public-api.html",
    recipes = "articles/index.html"
  )
  utils::browseURL(paste0(url_base, path))
  invisible(TRUE)
}


#' A consistent ggplot2 theme for spesim teaching plots
#'
#' @param base_size Base font size.
#' @export
theme_spesim <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold")
    )
}
