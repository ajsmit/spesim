# Tiled (Systematic Cell) Quadrat Placement

Places non-overlapping rectangular quadrats using a regular grid that
tiles the extent of the input domain. Candidate cells are generated with
sf::st_make_grid() using the requested quadrat width and height. Cells
that fall entirely inside the domain (via sf::st_within(..., sparse =
FALSE)) are kept as valid locations, from which up to n_quadrats are
sampled at random (set a seed for reproducibility). Because candidates
come from a grid, selected quadrats never overlap.

## Usage

``` r
place_quadrats_tiled(domain, n_quadrats, quadrat_size)
```

## Arguments

- domain:

  An sf polygon or multipolygon object defining the sampling region.
  Must have a valid CRS suitable for linear measurements (projected, not
  geographic).

- n_quadrats:

  Integer (≥1). The target number of quadrats to place. If fewer valid
  grid cells exist inside domain, the function will place the maximum
  possible and adjust this value downward with a warning.

- quadrat_size:

  Numeric vector of length 2, c(width, height), giving the quadrat
  dimensions in the same units as domain (e.g., meters). Both values
  must be positive and finite.

## Value

An sf object (polygons) with two columns:

- quadrat_id:

  Sequential integer identifier of the placed quadrats.

- geometry:

  Polygon geometry of each quadrat.

## Details

- The function assumes a planar (projected) CRS where distances and
  areas are in linear units; if your domain is in longitude/latitude,
  reproject it first (e.g., to UTM).

- quadrat_size is passed to sf::st_make_grid(cellsize = ...), so the two
  numbers are interpreted as width and height in the domain’s units.

- If fewer than n_quadrats valid cells fit entirely inside the domain,
  the function places as many as possible and emits a warning.

- Returned quadrats are labeled sequentially in the order they are
  sampled (quadrat_id = 1, 2, …, n_quadrats).

## See also

[`place_quadrats`](https://ajsmit.github.io/spesim/reference/place_quadrats.md)
(random non-overlapping),
[`place_quadrats_systematic`](https://ajsmit.github.io/spesim/reference/place_quadrats_systematic.md)
(grid centers within domain),
[`place_quadrats_transect`](https://ajsmit.github.io/spesim/reference/place_quadrats_transect.md)
(along parallel transects),
[`place_quadrats_voronoi`](https://ajsmit.github.io/spesim/reference/place_quadrats_voronoi.md)
(via Voronoi cells)

## Examples

``` r
if (FALSE) { # \dontrun{
library(sf)
set.seed(1)
dom <- create_sampling_domain()
qs  <- place_quadrats_tiled(dom, n_quadrats = 12, quadrat_size = c(1.5, 1.5))
plot(st_geometry(dom), col = "grey95", border = "grey60")
plot(st_geometry(qs), add = TRUE, border = "black", lwd = 1)
} # }
```
