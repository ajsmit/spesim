# Run a method-testing bundle (simulate + audit + standard outputs)

Convenience helper for teaching and **method testing**.

`spesim_method_test()` runs a simulation, computes a lightweight audit
of the realised regime (via
[`spesim_audit()`](https://ajsmit.github.io/spesim/reference/spesim_audit.md)),
and optionally returns a small set of standard analysis tables and
plots.

Compared to
[`spesim_demo()`](https://ajsmit.github.io/spesim/reference/spesim_demo.md),
this helper is oriented toward *auditing and comparing settings* (e.g.,
sampling schemes) rather than a curated "basic → advanced" teaching
progression.

## Usage

``` r
spesim_method_test(
  init_file = NULL,
  P = NULL,
  seed = NULL,
  write_outputs = FALSE,
  output_prefix = NULL,
  diagnostics = "nn",
  make_tables = TRUE,
  make_plots = TRUE,
  ...
)
```

## Arguments

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

- diagnostics:

  Diagnostics to compute in
  [`spesim_audit()`](https://ajsmit.github.io/spesim/reference/spesim_audit.md).
  Default `"nn"` for deterministic, dependency-minimal behaviour.

- make_tables:

  Logical; if `TRUE`, include standard analysis tables.

- make_plots:

  Logical; if `TRUE`, include standard plots.

- ...:

  Additional arguments passed to
  [`run_spatial_simulation()`](https://ajsmit.github.io/spesim/reference/run_spatial_simulation.md).

## Value

A list with elements:

- `res`: results list returned by
  [`run_spatial_simulation()`](https://ajsmit.github.io/spesim/reference/run_spatial_simulation.md).

- `audit`: audit list returned by
  [`spesim_audit()`](https://ajsmit.github.io/spesim/reference/spesim_audit.md).

- `tables` (optional): named list of analysis tables (when
  `make_tables = TRUE`).

- `plots` (optional): named list of plots (when `make_plots = TRUE`).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- spesim_method_test(make_plots = TRUE)
x$audit
x$plots$spatial
} # }
```
