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
#' @param label_ids Logical or character. If `TRUE`, label each quadrat using the
#'   column given by `id_col`. If `FALSE`, do not label. If `"auto"` (default),
#'   labels are shown only when the number of quadrats is at most
#'   `label_threshold`.
#' @param label_threshold Integer. Only used when `label_ids = "auto"`.
#' @param id_col Character scalar naming the quadrat id column in `quadrats`.
#'   Defaults to `"quadrat_id"`.
#' @param show_voronoi Logical or character. If `TRUE`, and the `quadrats` object
#'   carries a `"voronoi_cells"` attribute (as returned by
#'   [place_quadrats_voronoi()] with `show_voronoi = TRUE`), plot those Voronoi cell
#'   boundaries underneath the quadrats. If `"auto"` (default), Voronoi overlays
#'   are shown only when the needed attributes are present.
#' @param show_seeds Logical or character. If `TRUE` and Voronoi overlays are
#'   enabled, also plot the Voronoi seed points if present as a `"voronoi_seeds"`
#'   attribute. If `"auto"` (default), seeds are shown when present.
#' @param domain_colour,quadrat_colour Character scalars giving line colours for
#'   the domain outline and quadrat outlines.
#' @param quadrat_fill Character scalar giving the quadrat fill colour.
#' @param quadrat_alpha Numeric in (0,1]. Alpha transparency for quadrat fill.
#' @param domain_linewidth,quadrat_linewidth Numeric line widths.
#' @param voronoi_linewidth Numeric line width for Voronoi cell boundaries.
#' @param label_size Numeric text size for quadrat id labels.
#' @param base_size Base font size passed to `ggplot2::theme_minimal()`.
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
  label_ids = c("auto", "yes", "no"),
  label_threshold = 40,
  id_col = "quadrat_id",
  show_voronoi = c("auto", "yes", "no"),
  show_seeds = c("auto", "yes", "no"),
  domain_colour = "#22223b",
  quadrat_colour = "#22223b",
  quadrat_fill = "#4a4e69",
  quadrat_alpha = 0.18,
  domain_linewidth = 0.7,
  quadrat_linewidth = 0.6,
  voronoi_linewidth = 0.5,
  label_size = 3,
  base_size = 12
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_quadrats() requires the 'ggplot2' package. Please install it.")
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("plot_quadrats() requires the 'sf' package. Please install it.")
  }

  if (!inherits(domain, "sf")) stop("'domain' must be an sf object.")
  if (!inherits(quadrats, "sf")) stop("'quadrats' must be an sf object.")

  label_ids <- match.arg(label_ids)
  show_labels <- switch(
    label_ids,
    auto = nrow(quadrats) <= label_threshold,
    yes = TRUE,
    no = FALSE
  )

  # Backwards-compatible handling for show_voronoi/show_seeds
  if (is.logical(show_voronoi)) show_voronoi <- if (isTRUE(show_voronoi)) "yes" else "no"
  if (is.logical(show_seeds)) show_seeds <- if (isTRUE(show_seeds)) "yes" else "no"
  show_voronoi <- match.arg(show_voronoi)
  show_seeds <- match.arg(show_seeds)

  # Handle empty quadrats gracefully
  if (nrow(quadrats) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = domain, fill = NA, colour = domain_colour, linewidth = domain_linewidth) +
      ggplot2::coord_sf(datum = NA, expand = FALSE) +
      theme_spesim(base_size = base_size)
    if (!is.null(title)) p <- p + ggplot2::labs(title = title)
    return(p)
  }

  # Optional Voronoi background
  vor_cells <- attr(quadrats, "voronoi_cells")
  seeds <- attr(quadrats, "voronoi_seeds")

  show_vor <- switch(
    show_voronoi,
    auto = !is.null(vor_cells),
    yes = TRUE,
    no = FALSE
  )
  show_seed_pts <- switch(
    show_seeds,
    auto = !is.null(seeds),
    yes = TRUE,
    no = FALSE
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = domain, fill = NA, colour = domain_colour, linewidth = domain_linewidth)

  if (isTRUE(show_vor) && !is.null(vor_cells)) {
    p <- p + ggplot2::geom_sf(data = vor_cells, fill = NA, colour = "grey70", linewidth = voronoi_linewidth)
  }
  if (isTRUE(show_vor) && isTRUE(show_seed_pts) && !is.null(seeds)) {
    p <- p + ggplot2::geom_sf(data = sf::st_as_sf(seeds), colour = "grey40", alpha = 0.7, size = 0.8)
  }

  p <- p +
    ggplot2::geom_sf(
      data = quadrats,
      fill = quadrat_fill,
      alpha = quadrat_alpha,
      colour = quadrat_colour,
      linewidth = quadrat_linewidth
    ) +
    ggplot2::coord_sf(datum = NA, expand = FALSE) +
    theme_spesim(base_size = base_size)

  if (!is.null(title)) p <- p + ggplot2::labs(title = title)

  if (isTRUE(show_labels)) {
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
