# Random Non‑Overlapping Quadrat Placement

Places axis‑aligned rectangular quadrats at random locations inside
domain, rejecting candidates that either (i) are not fully contained or
(ii) overlap any previously accepted quadrat. A simple
rejection‑sampling loop is used with a finite attempt budget.

## Usage

``` r
place_quadrats(domain, n_quadrats, quadrat_size)
```

## Arguments

- domain:

  An sf polygon/multipolygon defining the sampling region. Must carry a
  valid projected CRS suitable for linear units.

- n_quadrats:

  Integer (≥1). Target number of non‑overlapping quadrats to place.

- quadrat_size:

  Numeric vector of length 2, c(width, height), giving quadrat
  dimensions in the same units as the domain CRS. Both values must be
  positive and finite.

## Value

An sf object (polygons) with:

- quadrat_id:

  Sequential integer identifier for each placed quadrat.

- geometry:

  Axis‑aligned rectangular polygon for the quadrat.

If no valid placement is possible, returns an empty sf with the same CRS
and emits a warning.

## Details

- Each candidate center is drawn uniformly from the domain’s bounding
  box, shrunk so the full rectangle of size quadrat_size = c(width,
  height) fits inside. The resulting polygon is accepted only if it is
  entirely within domain and non‑overlapping with all accepted quadrats.

- The loop stops when n_quadrats are placed or after n_quadrats \* 100
  attempts (a warning is issued if the target is not met).

- CRS: Use a projected CRS (e.g., meters). quadrat_size is interpreted
  in the CRS linear units.

- Reproducibility: Set a seed (e.g., set.seed(...)) before calling to
  obtain repeatable placements.

## See also

[`place_quadrats_systematic`](https://ajsmit.github.io/spesim/reference/place_quadrats_systematic.md)
(regular grid),
[`place_quadrats_tiled`](https://ajsmit.github.io/spesim/reference/place_quadrats_tiled.md)
(systematic tiling by cell size),
[`place_quadrats_voronoi`](https://ajsmit.github.io/spesim/reference/place_quadrats_voronoi.md)
(Voronoi‑based centers),
[`place_quadrats_transect`](https://ajsmit.github.io/spesim/reference/place_quadrats_transect.md)
(parallel transects)

## Examples
