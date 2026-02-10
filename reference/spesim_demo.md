# Run a compact spesim demonstration workflow

Convenience helper for teaching and quick exploration. It runs a small
simulation (in-memory by default) and optionally returns a set of
ready-made plots.

The goal is to provide a single, discoverable entry point for the common
“basic → intermediate → advanced” workflow described in the workflow
vignette.

## Usage

``` r
spesim_demo(
  level = c("basic", "intermediate", "advanced"),
  init_file = NULL,
  P = NULL,
  seed = NULL,
  write_outputs = FALSE,
  output_prefix = NULL,
  make_plots = FALSE,
  ...
)
```

## Arguments

- level:

  Character scalar. One of `"basic"`, `"intermediate"`, or `"advanced"`.

- init_file:

  Optional path to an init file. If `NULL`, uses the packaged example
  `examples/spesim_init_basic.txt`.

- P:

  Optional parameter list as returned by
  [`load_config()`](https://ajsmit.github.io/spesim/reference/load_config.md).
  If supplied, it takes precedence over `init_file`.

- seed:

  Optional integer seed. If provided, overrides any seed in `P`.

- write_outputs:

  Logical. Passed to
  [`run_spatial_simulation()`](https://ajsmit.github.io/spesim/reference/run_spatial_simulation.md).
  Default `FALSE` (fast; no files written).

- output_prefix:

  Optional output prefix passed to
  [`run_spatial_simulation()`](https://ajsmit.github.io/spesim/reference/run_spatial_simulation.md).
  Only used when `write_outputs = TRUE`.

- make_plots:

  Logical. If `TRUE`, returns a list of ggplot/patchwork objects in
  `out$plots`.

- ...:

  Additional arguments passed to
  [`run_spatial_simulation()`](https://ajsmit.github.io/spesim/reference/run_spatial_simulation.md).

## Value

A list with elements:

- `res`: the results list returned by
  [`run_spatial_simulation()`](https://ajsmit.github.io/spesim/reference/run_spatial_simulation.md).

- `plots` (optional): named list of plots (when `make_plots = TRUE`).

- `tables` (optional): named list of intermediate analysis tables (when
  `level != "basic"`).

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic (fast): quadrats + individuals
x <- spesim_demo("basic", make_plots = TRUE)
x$plots$quadrats

# Intermediate: adds classic teaching plots
y <- spesim_demo("intermediate", make_plots = TRUE)
y$plots$rank_abundance

# Advanced: adds the full multi-panel diagnostic
z <- spesim_demo("advanced", make_plots = TRUE)
z$plots$advanced_panel
} # }
```
