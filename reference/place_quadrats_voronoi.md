# Voronoi‑based Quadrat Placement

Places rectangular quadrats by first generating a dense set of random
seed points inside the domain, computing the Voronoi tessellation, and
then identifying those Voronoi cells whose largest inscribed circle is
big enough to contain the requested quadrat (by diagonal). For each
suitable cell, the quadrat is placed at the cell’s “pole of
inaccessibility” (center of the inscribed circle), which tends to yield
well‑spaced sampling locations even in irregular polygons.

## Usage

``` r
place_quadrats_voronoi(
  domain,
  n_quadrats,
  quadrat_size,
  voronoi_seed_factor,
  show_voronoi = FALSE
)
```

## Arguments

- domain:

  An sf polygon/multipolygon defining the sampling region. Must carry a
  valid projected CRS suitable for linear measurements.

- n_quadrats:

  Integer (≥1). Target number of quadrats to place. The actual number
  returned may be smaller if the domain cannot accommodate enough
  suitable Voronoi cells.

- quadrat_size:

  Numeric vector of length 2, c(width, height), giving quadrat
  dimensions in the same linear units as domain (e.g., meters). Both
  values must be positive and finite.

- voronoi_seed_factor:

  Numeric (≥1 recommended). Multiplier controlling how many random seeds
  are used to build the Voronoi diagram: n_seeds = n_quadrats \*
  voronoi_seed_factor. Use larger values (e.g., 8–20) in very irregular
  or narrow domains to improve coverage.

- show_voronoi:

  Logical. If TRUE, the returned sf object (quadrats) will carry
  additional attributes that can be used to visualise the underlying
  Voronoi tessellation used to select candidate cells.

## Value

An sf object (polygons) with:

- quadrat_id:

  Sequential integer identifier of the placed quadrats.

- geometry:

  Axis‑aligned rectangular polygon for each quadrat.

If `show_voronoi = TRUE`, the returned object also has attributes:

- voronoi_cells:

  An sf object with the Voronoi polygons clipped to the domain.

- voronoi_seeds:

  An sfc POINT collection of the random seed points.

- inscribed_circles:

  An sfc POLYGON collection of inscribed-circle polygons.

## Details

- The function assumes a projected (planar) CRS for domain so that
  distances/areas are in linear units. If your data are in
  longitude/latitude, reproject (e.g., to UTM) before calling.

- The number of initial seeds is n_quadrats \* voronoi_seed_factor.
  Larger values explore the domain more densely and typically find more
  suitable cells, at the cost of extra computation.

- Suitability is decided by comparing each cell’s inscribed‑circle
  radius to half the requested quadrat diagonal: radius \>=
  sqrt(width^2 + height^2) / 2.

- If fewer than n_quadrats suitable cells exist, all suitable cells are
  used and a warning is issued.

- Quadrats are axis‑aligned rectangles centered at the selected cell
  centers. Because the cell centers can be near each other, axis‑aligned
  rectangles may still overlap slightly in tight spaces—even when
  circles do not. If strict non‑overlap is required, post‑filtering is
  recommended.

- Internally uses sf::st_voronoi() for tessellation and
  st_inscribed_circle() (provided by lwgeom) to compute maximal
  in‑circle centers/radii.

## See also

[`place_quadrats`](https://ajsmit.github.io/spesim/reference/place_quadrats.md)
(random, non‑overlapping),
[`place_quadrats_tiled`](https://ajsmit.github.io/spesim/reference/place_quadrats_tiled.md)
(systematic tiling),
[`place_quadrats_systematic`](https://ajsmit.github.io/spesim/reference/place_quadrats_systematic.md)
(grid centers within domain),
[`place_quadrats_transect`](https://ajsmit.github.io/spesim/reference/place_quadrats_transect.md)
(parallel transects)

## Examples

``` r
if (FALSE) { # \dontrun{
library(sf)
set.seed(2)
dom <- create_sampling_domain()
qs  <- place_quadrats_voronoi(
  domain = dom,
  n_quadrats = 20,
  quadrat_size = c(1.5, 1.5),
  voronoi_seed_factor = 12,
  show_voronoi = TRUE
)
plot(st_geometry(dom), col = "grey95", border = "grey60")
plot(st_geometry(attr(qs, "voronoi_cells")), add = TRUE, border = "grey80")
plot(st_geometry(qs), add = TRUE, border = "black", lwd = 1)
} # }
```
