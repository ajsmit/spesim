# Run a complete spesim simulation (recommended)

High-level, teaching-friendly wrapper around
[`run_spatial_simulation()`](https://ajsmit.github.io/spesim/reference/run_spatial_simulation.md).
It supports both file-driven and programmatic workflows, adds optional
progress messaging, and returns an S3 result object with stable
components.

Use this as the main entry point in workshops and teaching material.

## Usage

``` r
spesim_run(
  config,
  domain = NULL,
  interactions_file = NULL,
  output_prefix = NULL,
  write_outputs = FALSE,
  seed = NULL,
  quiet = FALSE,
  ...
)
```

## Arguments

- config:

  Either a path to an init file (character scalar), or an in-memory
  parameter list `P` (typically from
  [`load_config()`](https://ajsmit.github.io/spesim/reference/load_config.md)).

- domain:

  Optional `sf` polygon study area (see
  [`create_sampling_domain()`](https://ajsmit.github.io/spesim/reference/create_sampling_domain.md)).
  If `NULL`, a default domain is created.

- interactions_file:

  Optional path to an interactions config file.

- output_prefix:

  Base output prefix (timestamp is appended) when
  `write_outputs = TRUE`.

- write_outputs:

  Logical; write CSVs/figures/report to disk? Default `FALSE` for a
  smoother interactive experience.

- seed:

  Optional integer. If supplied, overrides `P$SEED` (and is applied via
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before simulation).

- quiet:

  Logical; suppress most messages.

- ...:

  Passed through to
  [`run_spatial_simulation()`](https://ajsmit.github.io/spesim/reference/run_spatial_simulation.md).

## Value

An object of class `spesim_result` (a named list) with components:

- `P`: resolved parameters

- `domain`: study area (`sf`)

- `species_dist`: simulated individuals (`sf` points)

- `quadrats`: sampling quadrats (`sf` polygons; class `spesim_quadrats`)

- `env_gradients`: environmental grid

- `abund_matrix`: site × species abundances

- `site_coords`: quadrat centroids

## See also

[`spesim_demo()`](https://ajsmit.github.io/spesim/reference/spesim_demo.md),
[`load_config()`](https://ajsmit.github.io/spesim/reference/load_config.md),
[`run_spatial_simulation()`](https://ajsmit.github.io/spesim/reference/run_spatial_simulation.md)
