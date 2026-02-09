# Build a site-by-species abundance matrix from point–polygon overlaps

Aggregates simulated individuals (`species_dist`, POINTS) into sampling
units (`quadrats`, POLYGONS) and returns a wide site × species table
(one row per quadrat; one column per species). Counts are computed via
spatial intersection; species absent from a quadrat are filled with
zeros.

## Usage

``` r
create_abundance_matrix(species_dist, quadrats, all_species_names)
```

## Arguments

- species_dist:

  An `sf` POINT layer of individuals that contains a character (or
  factor) column named `species`. Geometry must be in the same CRS as
  `quadrats`.

- quadrats:

  An `sf` POLYGON (or MULTIPOLYGON) layer of sampling units that
  contains an integer (or character) column named `quadrat_id`
  identifying each site. Geometry must be in the same CRS as
  `species_dist`.

- all_species_names:

  Character vector giving the complete set of species to include as
  columns (e.g., `LETTERS[1:S]`). Any species not present in the
  intersections are still created as columns and zero-filled.

## Value

A base `data.frame` with one row per quadrat and columns:

- `site`:

  Copy of `quadrat_id`.

- `<species>`:

  Non‑negative integer counts for each species in `all_species_names`.

## Details

The function uses
[`st_intersection()`](https://r-spatial.github.io/sf/reference/geos_binary_ops.html)
to associate each individual with the quadrat polygon that contains it,
then tallies counts by `(quadrat_id, species)` and pivots to wide
format. If no points fall in any quadrat, the result is a data frame
with one row per quadrat and zeros in all species columns. Missing
species columns (relative to `all_species_names`) are added and
zero-filled to ensure a consistent schema.

## CRS

Inputs must share the same coordinate reference system. If they differ,
reproject beforehand (e.g.,
`species_dist <- sf::st_transform(species_dist, sf::st_crs(quadrats))`).

## See also

[`st_intersection`](https://r-spatial.github.io/sf/reference/geos_binary_ops.html),
[`pivot_wider`](https://tidyr.tidyverse.org/reference/pivot_wider.html),
[`calculate_quadrat_environment`](https://ajsmit.github.io/spesim/reference/calculate_quadrat_environment.md)

## Examples

``` r
if (FALSE) { # \dontrun{
spp <- simulate_points() # must contain geometry + 'species'
quads <- make_quadrats() # must contain geometry + 'quadrat_id'
A <- create_abundance_matrix(spp, quads, all_species_names = LETTERS[1:10])
head(A)
} # }
```
