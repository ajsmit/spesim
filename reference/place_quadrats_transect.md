# Parallel‑Transect Quadrat Placement

Places axis‑aligned rectangular quadrats at evenly spaced positions
along a set of parallel transect lines crossing the domain at a
specified compass angle. The function first shrinks the domain by half
the quadrat diagonal to form a safe sampling area, ensuring that every
returned quadrat lies fully inside the original domain.

## Usage

``` r
place_quadrats_transect(
  domain,
  n_transects,
  n_quadrats_per_transect,
  quadrat_size,
  angle
)
```

## Arguments

- domain:

  An sf polygon/multipolygon defining the sampling region. Must carry a
  valid projected CRS suitable for linear units (e.g., meters).

- n_transects:

  Integer (≥1). Number of parallel transect lines to generate.

- n_quadrats_per_transect:

  Integer (≥1). Number of quadrats to place along each transect after
  clipping to the safe area.

- quadrat_size:

  Numeric vector of length 2, c(width, height), giving quadrat
  dimensions in the same units as the domain CRS. Both values must be
  positive and finite.

- angle:

  Numeric. Compass bearing in degrees for transect orientation (0 =
  North/up, 90 = East/right). Values outside \[0, 360) are allowed and
  are interpreted modulo 360.

## Value

An sf object (polygons) with:

- quadrat_id:

  Sequential integer identifier for each placed quadrat.

- geometry:

  Axis‑aligned rectangular polygon for the quadrat.

## Details

- The transect orientation is given by angle in compass degrees (0 =
  North, 90 = East, 180 = South, 270 = West).

- n_transects parallel lines are positioned across the domain’s extent,
  then each line is cropped to the safe area. Along each retained
  segment, n_quadrats_per_transect centers are placed by linear
  interpolation from segment start to end (a single quadrat uses the
  midpoint).

- Each center generates an axis‑aligned rectangular quadrat of size
  quadrat_size = c(width, height) in the domain’s units.

- Use a projected (planar) CRS for domain (e.g., meters) so that
  distances and sizes are meaningful. If the domain is too small
  relative to quadrat_size, the safe area may be empty and no quadrats
  are returned.

## See also

[`place_quadrats`](https://ajsmit.github.io/spesim/reference/place_quadrats.md)
(random, non‑overlapping),
[`place_quadrats_systematic`](https://ajsmit.github.io/spesim/reference/place_quadrats_systematic.md)
(regular grid over bounding box),
[`place_quadrats_tiled`](https://ajsmit.github.io/spesim/reference/place_quadrats_tiled.md)
(systematic tiling by cell size),
[`place_quadrats_voronoi`](https://ajsmit.github.io/spesim/reference/place_quadrats_voronoi.md)
(Voronoi‑based centers)

## Examples

``` r
if (FALSE) { # \dontrun{
library(sf)
dom <- create_sampling_domain()
qs <- place_quadrats_transect(
  domain = dom,
  n_transects = 3,
  n_quadrats_per_transect = 6,
  quadrat_size = c(1.5, 1.5),
  angle = 45
)
plot(st_geometry(dom), col = "grey95", border = "grey60")
plot(st_geometry(qs), add = TRUE, border = "black", lwd = 1)
} # }
```
