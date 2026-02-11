# Audit environmental filtering (optima/tolerances vs realised abundances)

Audit environmental filtering (optima/tolerances vs realised abundances)

## Usage

``` r
audit_environmental_filtering(res)
```

## Arguments

- res:

  Simulation results list.

## Value

A data.frame; one row per gradient-responsive species, with columns:

- species:

  Species label.

- gradient:

  Gradient name (e.g. temperature/elevation/rainfall).

- n_sites:

  Number of sampled quadrats with finite environment + abundance.

- occupancy_sites:

  Number of quadrats with abundance \> 0 (a simple sparsity check).

- rho:

  Spearman correlation between abundance and closeness-to-optimum
  (higher is more consistent with filtering).

- slope:

  Slope from a simple linear model of abundance vs distance-from-optimum
  (sign is diagnostic only).

- qualitative:

  Teaching-friendly label (e.g. abundance_peaks_near_optimum,
  weak_or_no_signal, too_sparse_to_assess).
