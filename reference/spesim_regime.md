# Classify whether a spesim run matches the intended qualitative regime

This helper turns the numeric outputs of
[`spesim_audit()`](https://ajsmit.github.io/spesim/reference/spesim_audit.md)
into a compact, teaching-friendly judgement (green/amber/red) about
whether the realised simulation matches the qualitative regime the user
believes they are simulating.

The defaults are deliberately conservative: when sample sizes are too
small or signals are weak, the result tends toward **amber** rather than
over- confident green/red.

## Usage

``` r
spesim_regime(res, audit = NULL, thresholds = list())
```

## Arguments

- res:

  A `spesim_result` list (as returned by
  [`spesim_run()`](https://ajsmit.github.io/spesim/reference/spesim_run.md)).

- audit:

  Optional list returned by
  [`spesim_audit()`](https://ajsmit.github.io/spesim/reference/spesim_audit.md).
  If `NULL`, the function will compute an audit internally.

- thresholds:

  Optional named list overriding default thresholds. See Details.

## Value

A named list with sub-lists `spatial`, `filtering`, `sampling`. Each
contains:

- `grade`: one of `"green"`, `"amber"`, `"red"`, `"unknown"`

- `reason`: short explanation

- `metrics`: compact numeric context

## Details

Default thresholds (you can override via `thresholds=`):

- Spatial NN ratio vs CSR (per species; uses dominant A if present):

  - `clustered`: nn_ratio \< 0.85

  - `CSR_like`: 0.85 ≤ nn_ratio ≤ 1.15

  - `inhibited`: nn_ratio \> 1.15

  - minimum points for classification: `min_n = 50`

- Filtering (per gradient-responsive species):

  - minimum occupied sites: `min_occ = 5`

  - `good`: Spearman rho(-\|z\|, abundance) \> 0.2 AND slope(abund~z) \<
    0

  - `opposite`: rho \< -0.2 OR slope \> 0

- Sampling constraint:

  - if boundary exclusion estimate \> 50% → amber/red (scheme-dependent)

  - if random placement attempts are close to the cap (≥90% of max
    attempts) → amber
