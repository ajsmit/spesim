#' Plot the spatial sampling simulation (domain, individuals, quadrats, optional gradient)
#'
#' @description
#' Creates a publication‑quality map of the simulated landscape. The domain outline is
#' drawn first; simulated individuals (POINT geometries) are overlaid and colored by
#' species; sampling quadrats (POLYGONs) are drawn with IDs. Optionally an environmental
#' raster/heat field is shown underneath for context.
#'
#' @param domain An \code{sf} polygon or multipolygon representing the study area.
#'   Geometry should be valid and in the same CRS as \code{species} and \code{quadrats}.
#'
#' @param species An \code{sf} POINT object of individuals with a character/factor column
#'   named \code{species}. All points must lie within \code{domain} for correct visual
#'   interpretation (the function does not clip).
#'
#' @param quadrats An \code{sf} POLYGON/MULTIPOLYGON object with a numeric/integer column
#'   \code{quadrat_id}. Quadrats are drawn as outlines and labeled at their centroids.
#'
#' @param P A named list of plotting/summary parameters. The following fields are read
#'   if present (each has a default if missing):
#'   \itemize{
#'     \item \code{POINT_SIZE} (\code{numeric}, default \code{0.2}) — point size for individuals,
#'     \item \code{POINT_ALPHA} (\code{numeric} in \eqn{(0,1)}, default \code{1}) — point transparency,
#'     \item \code{QUADRAT_COLOUR} (\code{character}, default \code{“black”}) — quadrat outline/label colour,
#'     \item \code{BACKGROUND_COLOUR} (\code{character}, default \code{“white”}) — plot background,
#'     \item \code{FOREGROUND_COLOUR} (\code{character}, default \code{”#22223b”}) — domain outline/title colour,
#'     \item \code{N_SPECIES} (\code{integer}) — used to size the palette; inferred from data if absent,
#'     \item \code{N_INDIVIDUALS} (\code{integer}) — used only for the subtitle count if available.
#'   }
#'
#' @param show_gradient Logical; if \code{TRUE}, an environmental surface is drawn beneath the
#'   geometries using \code{geom_raster()}. Default \code{FALSE}.
#'
#' @param env_gradients \emph{Required when} \code{show_gradient = TRUE}. A data frame with columns
#'   \code{x}, \code{y} (grid coordinates, in the same CRS units as \code{domain}) and one or more
#'   numeric columns representing environmental variables (e.g., \code{temperature_C},
#'   \code{elevation_m}, \code{rainfall_mm}). The column named by \code{gradient_type} must exist and be numeric.
#'
#' @param gradient_type Character scalar naming the column in \code{env_gradients} to plot when
#'   \code{show_gradient = TRUE}. Default \code{“temperature_C”}. The legend title is derived from this
#'   name (underscores replaced with spaces and title‑cased).
#'
#' @return A \code{ggplot} object that can be further modified (themes, scales, etc.).
#'
#' @details
#' Colors for species are drawn from \pkg{colorspace} \code{sequential_hcl} palette “RdPu” (reversed)
#' and mapped to the unique species present. All layers are plotted in the current display CRS of
#' the provided \code{sf} objects; ensure consistent CRS across inputs. When \code{show_gradient = TRUE},
#' the raster is drawn first with partial transparency so the domain outline and points remain visible.
#'
#' @seealso \code{\link{run_spatial_simulation}()}, \code{\link{create_sampling_domain}()},
#'   \code{\link{create_environmental_gradients}()}
#'
#' @examples
#' \dontrun{
#' p <- plot_spatial_sampling(domain, species_sf, quadrats_sf, P,
#'                            show_gradient = TRUE,
#'                            env_gradients = env_df,
#'                            gradient_type = “elevation_m”)
#' p + ggplot2::theme(legend.position = “right”)
#' }
#' @export
plot_spatial_sampling <- function(domain,
                                  species,
                                  quadrats,
                                  P,
                                  show_gradient = FALSE,
                                  env_gradients = NULL,
                                  gradient_type = "temperature_C") {
  # --- safe pulls from P with defaults
  POINT_SIZE <- if (!is.null(P$POINT_SIZE)) as.numeric(P$POINT_SIZE) else 0.2
  POINT_ALPHA <- if (!is.null(P$POINT_ALPHA)) as.numeric(P$POINT_ALPHA) else 1.0
  QUADRAT_COLOUR <- if (!is.null(P$QUADRAT_COLOUR)) as.character(P$QUADRAT_COLOUR) else "black"
  BACKGROUND_COLOUR <- if (!is.null(P$BACKGROUND_COLOUR)) as.character(P$BACKGROUND_COLOUR) else "white"
  FOREGROUND_COLOUR <- if (!is.null(P$FOREGROUND_COLOUR)) as.character(P$FOREGROUND_COLOUR) else "#22223b"
  N_SPECIES <- if (!is.null(P$N_SPECIES)) as.integer(P$N_SPECIES) else length(unique(species$species))

  # palette and labels
  spp_levels <- sort(unique(as.character(species$species)))
  if (length(spp_levels) < N_SPECIES) N_SPECIES <- length(spp_levels)
  pal <- rev(colorspace::sequential_hcl(N_SPECIES, palette = "RdPu"))
  names(pal) <- spp_levels

  # base ggplot with optional raster gradient underneath
  p <- ggplot()

  if (isTRUE(show_gradient)) {
    if (is.null(env_gradients) || !(gradient_type %in% names(env_gradients))) {
      stop(
        "When show_gradient = TRUE, provide env_gradients with a '",
        gradient_type, "' column."
      )
    }
    # build a small data.frame to avoid tidy-eval inside aes()
    gdf <- data.frame(
      x = env_gradients$x,
      y = env_gradients$y,
      value = env_gradients[[gradient_type]]
    )
    # use raster-like fill to avoid extra dependencies
    p <- p +
      geom_raster(
        data = gdf,
        mapping = aes(x = x, y = y, fill = value),
        alpha = 0.6,
        interpolate = TRUE
      ) +
      viridis::scale_fill_viridis(
        name = gsub("_", " ", tools::toTitleCase(gradient_type)),
        option = "viridis"
      )
  }

  # domain outline on top of gradient
  p <- p +
    geom_sf(data = domain, fill = "grey90", color = FOREGROUND_COLOUR, linewidth = 0.5, alpha = 0.4)

  # species points
  p <- p +
    geom_sf(
      data = species,
      mapping = aes(color = species),
      size = POINT_SIZE,
      alpha = POINT_ALPHA
    ) +
    scale_color_manual(values = pal, name = "Species", drop = FALSE)

  # quadrats (outline + id labels)
  p <- p +
    geom_sf(data = quadrats, fill = NA, color = QUADRAT_COLOUR, linewidth = 0.4) +
    geom_sf_text(
      data = quadrats,
      mapping = aes(label = quadrat_id),
      color = QUADRAT_COLOUR, size = 2.5, fontface = "bold"
    )

  # titles and theme
  title_text <- if (isTRUE(show_gradient)) {
    gsub("_", " ", tools::toTitleCase(gradient_type))
  } else {
    "Species Distribution"
  }
  subtitle_text <- paste0(
    if (!is.null(P$N_INDIVIDUALS)) P$N_INDIVIDUALS else nrow(species), " Individuals | ",
    length(spp_levels), " Species | ",
    nrow(quadrats), " Quadrats"
  )

  p +
    coord_sf(expand = FALSE) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = BACKGROUND_COLOUR, color = NA),
      plot.title = element_text(color = FOREGROUND_COLOUR, size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(color = FOREGROUND_COLOUR, size = 10, hjust = 0.5, margin = margin(b = 10))
    ) +
    labs(title = title_text, subtitle = subtitle_text)
}
