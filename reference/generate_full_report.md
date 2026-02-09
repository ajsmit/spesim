# (Advanced) Build a human-readable multi-section report

**Most users do not need to call this directly.** When you run
`run_spatial_simulation(write_outputs = TRUE)`, a report is written to
disk automatically under your timestamped `OUTPUT_PREFIX`. This function
returns the same report as a single character string for advanced
workflows (e.g., embedding the text in another system, tests, or custom
pipelines).

The report summarises environmental ranges, gradient correlations,
species abundance distribution, per-quadrat alpha diversity, diversity
partitioning, simple spatial autocorrelation, Fisher log-series
validation, and computation notes (e.g., whether fast Rcpp engines were
used for the point processes).

## Usage

``` r
generate_full_report(res)
```

## Arguments

- res:

  A named list produced by the simulator with at least the elements:

  `P`

  :   Parameter list returned by
      [`load_config`](https://ajsmit.github.io/spesim/reference/load_config.md)().
      Must include keys such as `N_INDIVIDUALS`, `N_SPECIES`,
      `DOMINANT_FRACTION`, `FISHER_ALPHA`, `FISHER_X`, and, if used, a
      tibble `GRADIENT` with columns `species`, `gradient` (one of
      "temperature","elevation","rainfall"), `optimum`, and `tol`.

  `env_gradients`

  :   A data frame with columns `temperature_C`, `elevation_m`, and
      `rainfall_mm` containing the gridded environmental fields used for
      context.

  `species_dist`

  :   An `sf` POINT layer of individuals with a `species` column. Used
      for tallies and alpha snapshots.

  `quadrats`

  :   An `sf` POLYGON layer with `quadrat_id`; used to compute per-site
      alpha summaries.

  `abund_matrix`

  :   A site \\\times\\ species abundance table (first column `site`;
      remaining columns are species counts), typically returned by
      [`create_abundance_matrix`](https://ajsmit.github.io/spesim/reference/create_abundance_matrix.md)().

  `site_coords`

  :   A data frame with columns `site`, `x`, `y` giving quadrat
      centroids in the same CRS used for analysis.

## Value

A single character scalar containing the full report text. No files are
written by this function; callers typically append the string to a log
or include it in the simulation report sink.

## Details

Sections produced:

- **Environmental Gradients:** min/max/range for temperature (deg C),
  elevation (m), rainfall (mm), with a short pattern note and a list of
  gradient-responsive species including their optima/tolerances rendered
  in natural units.

- **Gradient Correlations:** pairwise Pearson correlations among the
  three gradients and a brief interpretation (orthogonal vs.
  correlated).

- **Species Abundance Distribution:** ranked counts with percentage and
  labels for dominant and gradient-responsive species.

- **Spatial Alpha Diversity:** per-quadrat species richness and a
  compact species list (with counts) for each quadrat.

- **Diversity Partitioning:** mean alpha richness (+/-SE), Shannon's H',
  Simpson's 1-D, gamma richness, Whittaker and additive beta, mean
  pairwise Sorensen dissimilarity (computed as binary Bray-Curtis on
  presence-absence), and simple abundance dispersion metrics.

- **Spatial Autocorrelation:** Pearson correlation between inter-quadrat
  distances and differences in richness (a Mantel-style proxy) with
  p-value and interpretation.

- **Fisher Log-series Validation:** RMSE, R^2, max residual between
  observed ranks and theoretical log-series abundances; compares
  configured `FISHER_ALPHA` with
  [`fisher.alpha`](https://vegandevs.github.io/vegan/reference/diversity.html)
  estimated from the data.

- **Computation Notes:** which baseline point-process models were
  requested for the dominant and other species (`SPATIAL_PROCESS_A`,
  `SPATIAL_PROCESS_OTHERS`), and whether the fast Rcpp engines were used
  when applicable:

  - Thomas: `rthomas_bbox_cpp` (fast) or spatstat-based fallback

  - Strauss: `rstrauss_bbox_cpp` (fast) or spatstat-based fallback

Internally, this routine relies on base summaries, sf for spatial
intersections, and vegan for diversity indices. It avoids
tidy-evaluation in favor of explicit column access to keep dependencies
minimal within a non-interactive reporting context.

## Typical usage

Most users should read the saved report from disk rather than call this
function. See
[`read_latest_report()`](https://ajsmit.github.io/spesim/reference/read_latest_report.md)
for a convenient helper.

## Requirements

The `res` object must be consistent (all components correspond to the
same simulation and coordinate reference system). Missing or empty
inputs will lead to sections being populated with `NA` summaries or
informative defaults where possible; irrecoverable inconsistencies raise
errors.

## See also

[`run_spatial_simulation`](https://ajsmit.github.io/spesim/reference/run_spatial_simulation.md),
[`read_latest_report`](https://ajsmit.github.io/spesim/reference/read_latest_report.md),
[`create_abundance_matrix`](https://ajsmit.github.io/spesim/reference/create_abundance_matrix.md),
[`calculate_quadrat_environment`](https://ajsmit.github.io/spesim/reference/calculate_quadrat_environment.md)
