# Fast Thomas Process (Rcpp-backed) in an arbitrary polygon

Simulate \\n\\target\\ points from a Thomas (Gaussian Neyman–Scott)
process, using a fast C++ generator in the domain's bounding box and
then filtering to the polygon with sf. This is a drop-in replacement for
the spatstat-based Thomas simulator and is substantially faster for
large simulations.

## Usage

``` r
rthomas_fast(
  domain,
  n_target,
  kappa = NULL,
  mu = 10,
  sigma = 1,
  oversample = 1.3
)
```

## Arguments

- domain:

  An sf polygon/multipolygon defining the study area (projected CRS
  recommended).

- n_target:

  Integer, desired number of retained points inside `domain`.

- kappa:

  Optional parent intensity (parents per unit area). If `NULL`, it is
  derived from `n_target`, `mu`, and domain area.

- mu:

  Mean offspring per parent (Poisson).

- sigma:

  Gaussian cluster scale (standard deviation of offspring displacement;
  map units).

- oversample:

  Scalar \\\> 1\\. Multiplier used to request enough raw offspring in
  the bbox so that, after polygon filtering, you still retain about
  `n_target`. Default `1.3`.

## Value

An sf POINT layer with exactly `n_target` rows (if feasible), or fewer
if `kappa`/`mu` are too small to generate enough points.

## Details

Let \\A_D\\ be the polygon area and \\A_B\\ its bbox area. The expected
fraction of bbox points that survive the polygon filter is \\f \approx
A_D / A_B\\. We request roughly `n_target / f` points from C++ (times
`oversample`), then filter, and finally thin or (if short) resample with
replacement to return exactly `n_target`. If you prefer strict “no
replacement”, set `oversample` higher or supply a larger `kappa`.

## Examples

``` r
if (FALSE) { # \dontrun{
dom <- create_sampling_domain()
pts <- rthomas_fast(dom, n_target = 2000, mu = 10, sigma = 1)
plot(sf::st_geometry(pts), pch = 16, cex = 0.4)
} # }
```
