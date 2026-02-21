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
  # Performance-optimized version of the original function.
  #
  # The original implementation was correct but slow, primarily because
  # sf::st_intersection() is computationally expensive. It has to calculate
  # a full geometric intersection, allocate memory for the results, and then
  # dplyr::count() + pivot_wider() performs costly data reshaping.
  #
  # This version is significantly faster by using a sparse intersection matrix,
  # which avoids creating new geometries. It directly tallies abundances into a
  # pre-allocated matrix, which is a far more direct and memory-efficient way
  # to build a site-by-species table.
  #
  # Steps:
  # 1. Use sf::st_intersects() to get a sparse matrix indicating which points
  #    are in which quadrats. The result is a list of integer vectors, where
  #    each element `i` of the list contains the indices of the points inside
  #    quadrat `i`.
  # 2. Pre-allocate an abundance matrix with sites (quadrats) as rows and
  #    species as columns, initialized to all zeros.
  # 3. Iterate through the intersection list. For each quadrat, get the species
  #    of the points that fell inside it, and use table() to efficiently count
  #    the occurrences of each species.
  # 4. Update the corresponding row in the pre-allocated matrix with these counts.
  # 5. Convert the final matrix to a data.frame and ensure column names and
  #    types match the original function's output.

  # st_intersects is much faster as it doesn't create new geometries
  intersects <- sf::st_intersects(quadrats, species_dist)

  # Pre-allocate the abundance matrix for efficiency
  abund_matrix <- matrix(0,
    nrow = nrow(quadrats),
    ncol = length(all_species_names),
    dimnames = list(quadrats$quadrat_id, all_species_names)
  )

  # Get the species factor levels to ensure consistency
  species_factor <- factor(species_dist$species, levels = all_species_names)

  # Efficiently populate the matrix
  for (i in seq_along(intersects)) {
    point_indices <- intersects[[i]]
    if (length(point_indices) > 0) {
      species_in_quadrat <- species_factor[point_indices]
      # table() is a very fast way to count factor levels
      counts <- table(species_in_quadrat)
      abund_matrix[i, names(counts)] <- counts
    }
  }

  # Convert to data.frame, add site column, and ensure correct column order
  abund_df <- as.data.frame(abund_matrix)
  abund_df$site <- rownames(abund_df)
  
  # Reorder columns to have 'site' first, then all species names
  abund_df <- abund_df[, c("site", all_species_names)]
  
  # Ensure the site column is of the same type as quadrat_id and arrange
  abund_df$site <- as.character(abund_df$site)
  quadrat_ids_df <- data.frame(site = as.character(quadrats$quadrat_id))
  
  # Use a merge to ensure all sites are present and in the correct order
  res <- merge(quadrat_ids_df, abund_df, by = "site", all.x = TRUE)
  
  # Replace NA with 0 for sites that had no intersecting points
  res[is.na(res)] <- 0
  
  # Final arrange to match original output
  res <- res[order(res$site), ]

  # Ensure row names are reset
  rownames(res) <- NULL

  res
}

#' Summarise mean environmental conditions per quadrat
#'
#' @description
#' Converts a gridded environmental table (\code{env_grid}) into an \code{sf}
#' point layer, joins it to quadrat polygons, and computes per‑quadrat mean
#' values for all numeric environmental columns. Quadrats with no overlapping
#' grid points are retained and returned with \code{NA} means.
#'
#' @param env_grid A regular (or irregular) data frame with numeric columns
#'   \code{x}, \code{y} giving point coordinates, plus one or more numeric
#'   environmental columns. Coordinates are assumed to be in the CRS specified
#'   by \code{domain_crs}.
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
#'   \item groups by \code{quadrat_id} and returns the mean of each numeric
#'         environmental variable (with \code{na.rm = TRUE}),
#'   \item left‑joins back to the full set of quadrats to keep empty sites.
#' }
#' If you prefer a different join logic (e.g., nearest neighbour), adapt the
#' join call to specify a different predicate or use \code{st_nearest_feature}.
#'
#' @return
#' A \code{data.frame} with one row per quadrat and columns:
#' \code{site} plus one column per summarised numeric environmental variable.
#' Means are numeric; sites with no overlapping points will have \code{NA}.
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
  # Performance-optimized version of the original function.
  #
  # The original implementation used sf::st_join(), which is convenient but
  # can be slow for large datasets because it involves a full spatial join
  # operation. This creates a large intermediate data structure before the
  # summary statistics are computed.
  #
  # This version is faster because it avoids the full join.
  #
  # Steps:
  # 1. Convert the input grid to an `sf` object, same as before.
  # 2. Use sf::st_intersects() to find which grid points fall inside which
  #    quadrats. This returns a sparse list of indices, which is much more
  #    efficient than a full join.
  # 3. Identify the numeric columns from the environmental grid that need to be
  #    summarised.
  # 4. Pre-allocate a result matrix to store the mean environmental values for
  #    each quadrat. This is more memory-efficient than building up a result
  #    incrementally.
  # 5. Iterate over each quadrat. For each one, subset the environmental data
  #    to include only the points that fall inside it, and then efficiently
  #    calculate the column-wise means. This is faster than dplyr's
  #    group_by() + summarise() for this specific task.
  # 6. Convert the results matrix to a data.frame and ensure it includes all
  #    quadrats, filling with NA for any that had no overlapping points.

  env_sf <- sf::st_as_sf(env_grid, coords = c("x", "y"), crs = domain_crs)
  intersects <- sf::st_intersects(quadrats, env_sf)

  num_cols <- setdiff(names(env_grid), c("x", "y"))
  if (length(num_cols) == 0) {
    return(data.frame(site = quadrats$quadrat_id))
  }

  # Pre-allocate a matrix for the results for efficiency
  site_env_matrix <- matrix(NA_real_,
    nrow = nrow(quadrats),
    ncol = length(num_cols),
    dimnames = list(quadrats$quadrat_id, num_cols)
  )

  # Extract numeric data for faster subsetting
  env_data <- env_grid[, num_cols, drop = FALSE]

  # Loop through intersections and compute means
  for (i in seq_along(intersects)) {
    point_indices <- intersects[[i]]
    if (length(point_indices) > 0) {
      # colMeans is highly optimized for this task
      site_env_matrix[i, ] <- colMeans(env_data[point_indices, , drop = FALSE], na.rm = TRUE)
    }
  }

  # Convert to data.frame and add site column
  site_env_df <- as.data.frame(site_env_matrix)
  site_env_df$site <- rownames(site_env_df)

  # Ensure all sites from quadrats are present, even if they had no points.
  # quadrat_id may be integer; coerce to character to match create_abundance_matrix.
  res <- merge(data.frame(site = as.character(quadrats$quadrat_id)), site_env_df, by = "site", all.x = TRUE)
  
  # Reorder to match original output if necessary
  res <- res[, c("site", num_cols)]

  res
}
