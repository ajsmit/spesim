#' Build a site-by-species abundance matrix from point–polygon overlaps
#'
#' @description
#' Aggregates simulated individuals (\code{species_dist}, POINTS) into sampling
#' units (\code{quadrats}, POLYGONS) and returns a wide site × species table
#' (one row per quadrat; one column per species). Counts are computed via
#' spatial intersection; species absent from a quadrat are filled with zeros.
#'
#' @param species_dist An \code{sf} POINT layer of individuals that contains a
#'   character (or factor) column named \code{species}. Geometry must be in the
#'   same CRS as \code{quadrats}.
#' @param quadrats An \code{sf} POLYGON (or MULTIPOLYGON) layer of sampling
#'   units that contains an integer (or character) column named \code{quadrat_id}
#'   identifying each site. Geometry must be in the same CRS as \code{species_dist}.
#' @param all_species_names Character vector giving the complete set of species
#'   to include as columns (e.g., \code{LETTERS[1:S]}). Any species not present
#'   in the intersections are still created as columns and zero-filled.
#'
#' @details
#' The function uses \code{\link[sf]{st_intersection}()} to associate each
#' individual with the quadrat polygon that contains it, then tallies counts by
#' \code{(quadrat_id, species)} and pivots to wide format. If no points fall in
#' any quadrat, the result is a data frame with one row per quadrat and zeros in
#' all species columns. Missing species columns (relative to \code{all_species_names})
#' are added and zero-filled to ensure a consistent schema.
#'
#' @return
#' A base \code{data.frame} with one row per quadrat and columns:
#' \describe{
#'   \item{\code{site}}{Copy of \code{quadrat_id}.}
#'   \item{\code{<species>}}{Non‑negative integer counts for each species in
#'         \code{all_species_names}.}
#' }
#'
#' @section CRS:
#' Inputs must share the same coordinate reference system. If they differ,
#' reproject beforehand (e.g., \code{species_dist <- sf::st_transform(species_dist, sf::st_crs(quadrats))}).
#'
#' @seealso
#' \code{\link[sf]{st_intersection}}, \code{\link[tidyr]{pivot_wider}},
#' \code{\link{calculate_quadrat_environment}}
#'
#' @examples
#' \dontrun{
#' spp <- simulate_points() # must contain geometry + 'species'
#' quads <- make_quadrats() # must contain geometry + 'quadrat_id'
#' A <- create_abundance_matrix(spp, quads, all_species_names = LETTERS[1:10])
#' head(A)
#' }
#' @export
create_abundance_matrix <- function(species_dist, quadrats, all_species_names) {
  intersections <- sf::st_intersection(species_dist, quadrats)
  if (nrow(intersections) == 0) {
    abund_df <- data.frame(site = quadrats$quadrat_id)
    for (sp in all_species_names) abund_df[[sp]] <- 0
    return(abund_df)
  }
  abund_df <- intersections |>
    sf::st_drop_geometry() |>
    dplyr::count(quadrat_id, species, name = "abundance") |>
    tidyr::pivot_wider(names_from = species, values_from = abundance, values_fill = 0)

  missing_species <- setdiff(all_species_names, names(abund_df))
  for (sp in missing_species) abund_df[[sp]] <- 0

  data.frame(quadrat_id = quadrats$quadrat_id) |>
    dplyr::left_join(abund_df, by = "quadrat_id") |>
    dplyr::mutate(dplyr::across(-quadrat_id, ~ tidyr::replace_na(., 0))) |>
    dplyr::select(site = quadrat_id, dplyr::all_of(all_species_names)) |>
    dplyr::arrange(site)
}

#' Summarise mean environmental conditions per quadrat
#'
#' @description
#' Converts a gridded environmental table (\code{env_grid}) into an \code{sf}
#' point layer, joins it to quadrat polygons, and computes per‑quadrat mean
#' values (temperature, elevation, rainfall). Quadrats with no overlapping
#' grid points are retained and returned with \code{NA} means.
#'
#' @param env_grid A regular (or irregular) data frame with numeric columns
#'   \code{x}, \code{y} giving point coordinates, and environmental columns
#'   \code{temperature_C}, \code{elevation_m}, \code{rainfall_mm}. Coordinates
#'   are assumed to be in the CRS specified by \code{domain_crs}.
#' @param quadrats An \code{sf} POLYGON (or MULTIPOLYGON) layer with a
#'   \code{quadrat_id} column. Must be in the same CRS as \code{domain_crs}.
#' @param domain_crs A coordinate reference system for \code{env_grid} points.
#'   Can be an integer EPSG code, a PROJ4string/WKT, or an \code{sf} \code{crs}
#'   object. This CRS should match that of \code{quadrats}.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item converts \code{env_grid} to \code{sf} points via
#'         \code{\link[sf]{st_as_sf}} using \code{coords = c("x","y")},
#'   \item performs a spatial join \code{\link[sf]{st_join}} of points into
#'         quadrats (default predicate: \code{st_intersects}),
#'   \item groups by \code{quadrat_id} and returns the mean of each
#'         environmental variable (with \code{na.rm = TRUE}),
#'   \item left‑joins back to the full set of quadrats to keep empty sites.
#' }
#' If you prefer a different join logic (e.g., nearest neighbour), adapt the
#' join call to specify a different predicate or use \code{st_nearest_feature}.
#'
#' @return
#' A \code{data.frame} with one row per quadrat and columns:
#' \code{site}, \code{temperature_C}, \code{elevation_m}, \code{rainfall_mm}.
#' Means are numeric; sites with no overlapping points will have \code{NA} in
#' the environmental columns.
#'
#' @section CRS:
#' \code{env_grid} is interpreted in \code{domain_crs}; \code{quadrats} must
#' already be in the same CRS. Reproject beforehand as needed.
#'
#' @seealso
#' \code{\link[sf]{st_as_sf}}, \code{\link[sf]{st_join}}, \code{\link{create_environmental_gradients}}
#'
#' @examples
#' \dontrun{
#' env <- create_environmental_gradients(domain, resolution = 50, noise_level = 0.05)
#' E <- calculate_quadrat_environment(env, quadrats, domain_crs = sf::st_crs(domain))
#' head(E)
#' }
#' @export
calculate_quadrat_environment <- function(env_grid, quadrats, domain_crs) {
  env_sf <- sf::st_as_sf(env_grid, coords = c("x", "y"), crs = domain_crs)
  joined_data <- sf::st_join(quadrats, env_sf)
  site_env <- joined_data |>
    sf::st_drop_geometry() |>
    dplyr::group_by(site = quadrat_id) |>
    dplyr::summarise(
      dplyr::across(c(temperature_C, elevation_m, rainfall_mm), ~ mean(., na.rm = TRUE)),
      .groups = "drop"
    )
  dplyr::left_join(data.frame(site = quadrats$quadrat_id), site_env, by = "site")
}
