# Simulate locations via a Geyer-like saturation process (spatstat-free)

Implements a simple sequential acceptance sampler: propose candidate
points uniformly inside `domain`. Let `k` be the number of
already-accepted points within radius `r`. Accept with probability
\\\gamma^{\min(k, s)}\\, where \\\gamma\\ is the interaction parameter
and \\s\\ is the saturation count. For \\\gamma \> 1\\, the rule induces
clustering up to \\s\\ neighbours (values \\\< 1\\ yield inhibition).

## Usage

``` r
simulate_points_geyer(
  domain,
  n_target,
  beta = NULL,
  gamma = 1.5,
  r = 1,
  sat = 2
)
```

## Arguments

- domain:

  sf polygon/multipolygon.

- n_target:

  Integer target number of points.

- beta:

  Unused placeholder (kept for API compatibility).

- gamma:

  Interaction parameter (numeric). \\\gamma \> 1\\ encourages
  clustering; \\\gamma \< 1\\ encourages inhibition.

- r:

  Interaction radius (map units).

- sat:

  Saturation count (non-negative integer).

## Value

sf POINT layer (no attributes).
