# Estimate the effective sampling frame for quadrat centres

For an axis-aligned rectangle (quadrat) to lie fully within an irregular
polygon domain, the quadrat centre must lie inside a shrunk version of
the domain.

Exactly, this shrink is an erosion in an L-infinity metric (half-width
and half-height constraints). For simplicity and robustness in
teaching/method testing, this helper uses a conservative Euclidean
buffer based on the half-diagonal of the quadrat.

This provides an interpretable, single-number summary of boundary
constraints:

- `safe_area_fraction` = area(safe_domain) / area(domain)

Values close to 1 mean boundary exclusion is mild; values near 0 mean
that only a small interior region can host fully-contained quadrats.

## Usage

``` r
estimate_sampling_frame(domain, quadrat_size)
```

## Arguments

- domain:

  An `sf` polygon/multipolygon.

- quadrat_size:

  Numeric length-2, c(width, height), in the same units as the domain
  CRS.

## Value

A list with `safe_domain` (sf), `safe_area_fraction` (numeric), and
`buffer_dist` (numeric half-diagonal).
