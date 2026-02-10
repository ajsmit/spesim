# Systematic Grid Quadrat Placement

Places axis‑aligned rectangular quadrats at the centers of a regular
grid spanning the domain’s bounding box, and retains only those quadrats
whose full geometry lies inside the domain polygon. The grid density is
chosen to approximate n_quadrats total placements by respecting the
domain’s aspect ratio (so cells are roughly square in map space while
matching the requested count).

## Usage

``` r
place_quadrats_systematic(domain, n_quadrats, quadrat_size)
```

## Arguments

- domain:

  An sf polygon/multipolygon representing the sampling region. Must
  carry a valid projected CRS suitable for linear measurements.

- n_quadrats:

  Integer (≥1). Approximate number of quadrats to place. The actual
  number returned is the count of grid‑centered quadrats that fully fit
  inside domain.

- quadrat_size:

  Numeric vector of length 2, c(width, height), giving the quadrat
  dimensions in the same linear units as domain (e.g., meters). Both
  values must be positive and finite.

## Value

An sf object (polygons) with:

- quadrat_id:

  Sequential integer identifier for each placed quadrat.

- geometry:

  Axis‑aligned rectangular polygon for the quadrat.

## Details

- The grid is defined over the domain’s bounding box using nx × ny grid
  centers, where nx ≈ sqrt(n_quadrats / aspect_ratio) and ny ≈
  aspect_ratio × nx, with aspect_ratio = height / width of the bounding
  box. This produces roughly nx \* ny ≈ n_quadrats centers.

- For each grid center, an axis‑aligned rectangle of size quadrat_size =
  c(width, height) is created. Quadrats are kept only if they are fully
  within the domain (using sf::st_within()).

- If the domain is very irregular or narrow, the number of retained
  quadrats can be less than n_quadrats. The function returns all valid
  quadrats and emits a warning if none fit.

- The domain should use a projected (planar) CRS so that distances and
  sizes are in linear units (e.g., meters). If using lon/lat, reproject
  first (e.g., to UTM).

## See also

[`place_quadrats`](https://ajsmit.github.io/spesim/reference/place_quadrats.md)
(random, non‑overlapping),
[`place_quadrats_tiled`](https://ajsmit.github.io/spesim/reference/place_quadrats_tiled.md)
(systematic tiling by cell size),
[`place_quadrats_transect`](https://ajsmit.github.io/spesim/reference/place_quadrats_transect.md)
(parallel transects),
[`place_quadrats_voronoi`](https://ajsmit.github.io/spesim/reference/place_quadrats_voronoi.md)
(Voronoi‑based centers)

## Examples

``` r
if (FALSE) { # \dontrun{
library(sf)
dom <- create_sampling_domain()
qs  <- place_quadrats_systematic(
  domain = dom,
  n_quadrats = 24,
  quadrat_size = c(1.5, 1.5)
)
plot(st_geometry(dom), col = "grey95", border = "grey60")
plot(st_geometry(qs), add = TRUE, border = "black", lwd = 1)
} # }
```
