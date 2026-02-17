#' `spesim_result` objects
#'
#' @description
#' Objects returned by [spesim_run()]. A `spesim_result` is a named list with a
#' stable set of components representing a complete simulated community +
#' sampling design.
#'
#' @details
#' Components:
#' \describe{
#'   \item{`P`}{Resolved parameter list (see [load_config()]).}
#'   \item{`domain`}{`sf` polygon study area.}
#'   \item{`species_dist`}{`sf` POINT layer of individuals with a `species` column.}
#'   \item{`quadrats`}{`sf` POLYGON layer of sampling quadrats (class `spesim_quadrats`).}
#'   \item{`env_gradients`}{Data frame of gridded environmental fields.}
#'   \item{`abund_matrix`}{Site × species abundance table (first column `site`).}
#'   \item{`site_env`}{Per-quadrat mean environment table.}
#'   \item{`site_coords`}{Data frame with columns `site`, `x`, `y` (quadrat centroids).}
#'   \item{`files`}{Optional list of written file paths when `write_outputs = TRUE`.}
#' }
#'
#' Methods:
#' - `print()` and `summary()` methods are provided.
#'
#' @name spesim_result
#' @aliases spesim_quadrats summary_spesim_result
NULL
