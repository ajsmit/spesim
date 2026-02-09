# Summarise mean environmental conditions per quadrat

Converts a gridded environmental table (`env_grid`) into an `sf` point
layer, joins it to quadrat polygons, and computes per‑quadrat mean
values (temperature, elevation, rainfall). Quadrats with no overlapping
grid points are retained and returned with `NA` means.

## Usage

``` r
calculate_quadrat_environment(env_grid, quadrats, domain_crs)
```

## Arguments

- env_grid:

  A regular (or irregular) data frame with numeric columns `x`, `y`
  giving point coordinates, and environmental columns `temperature_C`,
  `elevation_m`, `rainfall_mm`. Coordinates are assumed to be in the CRS
  specified by `domain_crs`.

- quadrats:

  An `sf` POLYGON (or MULTIPOLYGON) layer with a `quadrat_id` column.
  Must be in the same CRS as `domain_crs`.

- domain_crs:

  A coordinate reference system for `env_grid` points. Can be an integer
  EPSG code, a PROJ4string/WKT, or an `sf` `crs` object. This CRS should
  match that of `quadrats`.

## Value

A `data.frame` with one row per quadrat and columns: `site`,
`temperature_C`, `elevation_m`, `rainfall_mm`. Means are numeric; sites
with no overlapping points will have `NA` in the environmental columns.

## Details

The function:

1.  converts `env_grid` to `sf` points via
    [`st_as_sf`](https://r-spatial.github.io/sf/reference/st_as_sf.html)
    using `coords = c("x","y")`,

2.  performs a spatial join
    [`st_join`](https://r-spatial.github.io/sf/reference/st_join.html)
    of points into quadrats (default predicate: `st_intersects`),

3.  groups by `quadrat_id` and returns the mean of each environmental
    variable (with `na.rm = TRUE`),

4.  left‑joins back to the full set of quadrats to keep empty sites.

If you prefer a different join logic (e.g., nearest neighbour), adapt
the join call to specify a different predicate or use
`st_nearest_feature`.

## CRS

`env_grid` is interpreted in `domain_crs`; `quadrats` must already be in
the same CRS. Reproject beforehand as needed.

## See also

[`st_as_sf`](https://r-spatial.github.io/sf/reference/st_as_sf.html),
[`st_join`](https://r-spatial.github.io/sf/reference/st_join.html),
[`create_environmental_gradients`](https://ajsmit.github.io/spesim/reference/create_environmental_gradients.md)

## Examples

``` r
if (FALSE) { # \dontrun{
env <- create_environmental_gradients(domain, resolution = 50, noise_level = 0.05)
E <- calculate_quadrat_environment(env, quadrats, domain_crs = sf::st_crs(domain))
head(E)
} # }
```
