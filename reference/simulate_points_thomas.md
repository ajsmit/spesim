# Simulate point locations via a Thomas (Gaussian Neyman–Scott) process

Parent–offspring clustering without external dependencies: parents are
sampled uniformly inside `domain`; each parent produces
\\\text{Pois}(\mu)\\ offspring displaced by isotropic Gaussian noise
with sd `sigma`. Offspring are clipped to the domain. The result is
thinned/replicated to match `n_target`.

## Usage

``` r
simulate_points_thomas(domain, n_target, kappa = NULL, mu = 10, sigma = 1)
```

## Arguments

- domain:

  sf polygon/multipolygon.

- n_target:

  Integer target number of points.

- kappa:

  Optional numeric parent intensity (parents per unit area). If `NULL`,
  a value is derived to hit `n_target` given `mu`.

- mu:

  Mean offspring per parent.

- sigma:

  Cluster scale (Gaussian sd, in map units).

## Value

sf POINT layer (no attributes).
