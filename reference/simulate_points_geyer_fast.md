# Fast Geyer saturation simulator (bbox Metropolis–Hastings)

Fast Geyer saturation simulator (bbox Metropolis–Hastings)

## Usage

``` r
simulate_points_geyer_fast(
  domain,
  n_target,
  r,
  gamma,
  sat,
  sweeps = 2000,
  burnin = 200,
  thin = 1
)
```

## Arguments

- domain:

  sf polygon/multipolygon defining the window (only bbox is used).

- n_target:

  integer, number of points to return.

- r:

  interaction radius.

- gamma:

  interaction parameter (\<1 inhibition, \>1 clustering).

- sat:

  saturation count (positive integer).

- sweeps:

  MH sweeps (per point).

- burnin:

  burn-in sweeps.

- thin:

  unused (reserved).

## Value

`sf` POINT layer with `n_target` points.
