#' Plot a sampling domain with quadrats
#'
#' @description
#' Convenience plotting helper for quadrat-placement workflows.
#' Draws the sampling domain outline and a set of quadrat polygons (optionally
#' labelled by quadrat id) using **ggplot2**.
#'
#' @param domain An sf polygon/multipolygon object representing the sampling
#'   domain.
#' @param quadrats An sf polygon/multipolygon object representing quadrats.
#'   Typically returned by [place_quadrats()], [place_quadrats_tiled()],
#'   [place_quadrats_systematic()], [place_quadrats_transect()], or
#'   [place_quadrats_voronoi()].
#' @param title Optional character scalar used as the plot title.
#' @param label_ids Logical. If `TRUE` (default), label each quadrat using the
#'   column given by `id_col`.
#' @param id_col Character scalar naming the quadrat id column in `quadrats`.
#'   Defaults to `"quadrat_id"`.
#' @param domain_colour,quadrat_colour Character scalars giving line colours for
#'   the domain outline and quadrat outlines.
#' @param quadrat_fill Character scalar giving the quadrat fill colour.
#' @param quadrat_alpha Numeric in (0,1]. Alpha transparency for quadrat fill.
#' @param domain_linewidth,quadrat_linewidth Numeric line widths.
#' @param label_size Numeric text size for quadrat id labels.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' library(spesim)
#' set.seed(1)
#' dom <- create_sampling_domain()
#' qs  <- place_quadrats(dom, n_quadrats = 12, quadrat_size = c(1.5, 1.5))
#' plot_quadrats(dom, qs, title = "Random quadrats")
#' }
#'
#' @export
plot_quadrats <- function(
  domain,
  quadrats,
  title = NULL,
  label_ids = TRUE,
  id_col = "quadrat_id",
  domain_colour = "#22223b",
  quadrat_colour = "#22223b",
  quadrat_fill = "#4a4e69",
  quadrat_alpha = 0.18,
  domain_linewidth = 0.7,
  quadrat_linewidth = 0.6,
  label_size = 3
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_quadrats() requires the 'ggplot2' package. Please install it.")
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("plot_quadrats() requires the 'sf' package. Please install it.")
  }

  if (!inherits(domain, "sf")) stop("'domain' must be an sf object.")
  if (!inherits(quadrats, "sf")) stop("'quadrats' must be an sf object.")

  # Handle empty quadrats gracefully
  if (nrow(quadrats) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = domain, fill = NA, colour = domain_colour, linewidth = domain_linewidth) +
      ggplot2::coord_sf(datum = NA) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
    if (!is.null(title)) p <- p + ggplot2::labs(title = title)
    return(p)
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = domain, fill = NA, colour = domain_colour, linewidth = domain_linewidth) +
    ggplot2::geom_sf(
      data = quadrats,
      fill = quadrat_fill,
      alpha = quadrat_alpha,
      colour = quadrat_colour,
      linewidth = quadrat_linewidth
    ) +
    ggplot2::coord_sf(datum = NA) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())

  if (!is.null(title)) p <- p + ggplot2::labs(title = title)

  if (isTRUE(label_ids)) {
    if (!(id_col %in% names(quadrats))) {
      stop("'quadrats' does not contain the id column '", id_col, "'.")
    }
    cent <- suppressWarnings(sf::st_centroid(quadrats))
    cent[["label"]] <- quadrats[[id_col]]
    p <- p + ggplot2::geom_sf_text(
      data = cent,
      ggplot2::aes(label = label),
      size = label_size,
      colour = quadrat_colour
    )
  }

  p
}
