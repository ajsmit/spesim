# Fast Strauss simulator (Rcpp backend) returning an sf point layer

Fixed-\\n\\ Metropolis–Hastings sampler for the Strauss process with
interaction radius `r` and parameter `gamma` (\\0\<\gamma\le 1\\). Uses
a cell grid, cached neighbor counts, squared-distance checks, and light
auto-tuning for the proposal radius. This is a fast, self-contained
alternative to spatstat for baseline point generation.

## Usage

``` r
simulate_points_strauss_fast(
  domain,
  n_target,
  r,
  gamma,
  sweeps = 2000,
  burnin = 200,
  thin = 1,
  step0 = NA_real_,
  target_acc = 0.35,
  tune_every = 200
)
```

## Arguments

- domain:

  An sf polygon/multipolygon defining the window.

- n_target:

  Integer; number of points to return.

- r:

  Interaction radius (map units).

- gamma:

  Interaction parameter (\\\le 1\\ for inhibition).

- sweeps:

  Number of MH sweeps (proposals per point). Default `2000`.

- burnin:

  Number of initial sweeps reserved for adaptation. Default `200`.

- thin:

  Currently unused (kept for API symmetry). Default `1`.

- step0:

  Optional initial proposal radius (defaults to `0.5 * r`).

- target_acc:

  Target acceptance for auto-tuning. Default `0.35`.

- tune_every:

  Sweeps between tuning updates. Default `200`.

## Value

An sf POINT layer with `n_target` rows.
