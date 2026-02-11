# Spatstat diagnostics for one or more species

Computes edge-corrected summary diagnostics (L(r)-r and pair-correlation
extrema) inside the true irregular domain window.

This is primarily intended for method testing and teaching: it helps
confirm whether the realised pattern is clustered/inhibited/CSR-like,
beyond a simple nearest-neighbour heuristic.

## Usage

``` r
spesim_spatstat_diagnostics(species_dist, domain, species = "all", rmax = NULL)
```

## Arguments

- species_dist:

  sf points with a `species` column.

- domain:

  sf polygon window.

- species:

  vector of species labels.

- rmax:

  maximum distance for summaries (optional).

## Value

data.frame with one row per species.
