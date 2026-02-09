# Simulate locations via a Strauss-like (soft inhibition) process

Implements a simple sequential acceptance sampler: propose candidate
points uniformly inside `domain`. Let `k` be the number of
already-accepted points within radius `r`. Accept with probability
\\\gamma^{k}\\ (with \\0 \< \gamma \le 1\\); otherwise reject. Repeat
until `n_target` points or a max iteration budget.

## Usage

``` r
simulate_points_strauss(domain, n_target, beta = NULL, gamma = 0.2, r = 1)
```

## Arguments

- domain:

  sf polygon/multipolygon.

- n_target:

  Integer target number of points.

- beta:

  Unused placeholder (kept for API compatibility).

- gamma:

  Inhibition strength (0–1). Smaller values increase inhibition.

- r:

  Interaction radius (map units).

## Value

sf POINT layer (no attributes).
