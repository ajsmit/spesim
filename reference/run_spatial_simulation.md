# Run the full spatial sampling simulation

End-to-end runner that can be called in two ways:

**A. File-driven (legacy / CLI-friendly)** Provide `init_file` (and
optionally `interactions_file`). The function reads configuration via
[`load_config()`](https://ajsmit.github.io/spesim/reference/load_config.md),
simulates the community, places quadrats, builds tables and figures,
writes them to disk under `output_prefix` (timestamped), and returns a
results list. **Interactions** can be supplied (i) via
`interactions_file`, (ii) by putting `INTERACTIONS_FILE = path/to/file`
in the init file, or (iii) by embedding inline rules with
`INTERACTIONS_EDGELIST = c("A,B-D,0.8", "C,A,1.2", "E,*,0.95")` directly
in the init file.

**B. Programmatic (as used in the README)** Provide an in-memory
parameter list `P` (typically from
[`load_config()`](https://ajsmit.github.io/spesim/reference/load_config.md))
and an `sf` polygon `domain` (or let the function create one).
Interactions are resolved from arguments or fields inside `P`: you may
pass an external `interactions_file`, set `P$INTERACTIONS_FILE`
(pointer), or set `P$INTERACTIONS_EDGELIST` (inline rules). By default
`write_outputs = TRUE`, so the same files are written; set
`write_outputs = FALSE` to skip writing and just get the results list
back.

Internally it:

1.  resolves parameters (`P`) and optional interspecific interactions,

2.  simulates a heterogeneous community (point process + environment),

3.  places sampling quadrats using the requested scheme,

4.  compiles site \\\times\\ species tables and per-site environment,

5.  (optionally) writes a 2\\\times\\2 map panel and an advanced
    multi-plot panel,

6.  (optionally) appends a human-readable text report, and

7.  returns a named results list.

## Usage

``` r
run_spatial_simulation(
  init_file = NULL,
  interactions_file = NULL,
  output_prefix = NULL,
  domain = NULL,
  P = NULL,
  write_outputs = TRUE,
  interactions_validate = TRUE,
  interactions_strict = TRUE,
  interactions_print = TRUE,
  interactions_top_n = 20
)
```

## Arguments

- init_file:

  Character path to a text configuration (KEY = value). If supplied,
  this takes precedence over `P`/`domain` inputs (file-driven mode). See
  [`load_config()`](https://ajsmit.github.io/spesim/reference/load_config.md)
  for supported keys.

- interactions_file:

  Optional character path for a separate interactions config. If `NULL`,
  the function will look inside the init file (or `P`) for either
  `INTERACTIONS_FILE` (pointer) or `INTERACTIONS_EDGELIST` (inline
  rules). Falls back to neutral interactions (all 1s; radius 0) if
  nothing is provided.

- output_prefix:

  Base path for outputs (timestamp appended). If `NULL`, derives from
  `P$OUTPUT_PREFIX` or defaults to `"simulation_output"`. Ignored when
  `write_outputs = FALSE`.

- domain:

  Optional `sf` polygon (study area). If `NULL`, a default domain is
  created by
  [`create_sampling_domain()`](https://ajsmit.github.io/spesim/reference/create_sampling_domain.md).
  Ignored in file-driven mode when `init_file` is provided (a domain is
  still created unless you supply your own).

- P:

  Optional fully materialized parameter list (often from
  [`load_config()`](https://ajsmit.github.io/spesim/reference/load_config.md)).
  If `init_file` is **not** given, `P` is used as-is.

- write_outputs:

  Logical; write CSVs/figures/report to disk? Default `TRUE`.

- interactions_validate:

  Logical; validate the resolved interaction matrix with
  [`validate_interactions()`](https://ajsmit.github.io/spesim/reference/validate_interactions.md)
  (default `TRUE`).

- interactions_strict:

  Logical; if validation finds problems, stop with an error (`TRUE`) or
  just warn (`FALSE`). Default `TRUE`.

- interactions_print:

  Logical; pretty-print a compact summary of the interaction matrix with
  [`print_interactions()`](https://ajsmit.github.io/spesim/reference/print_interactions.md)
  (default `TRUE`).

- interactions_top_n:

  Integer; cap the number of non-1.0 entries shown by the
  pretty-printer. Default `20`.

## Value

A named list with components:

- `P` – resolved parameters

- `domain` – `sf` polygon(s) of the study area

- `species_dist` – `sf` points with `species` and attached env columns

- `quadrats` – `sf` polygons of sampled quadrats

- `env_gradients` – data.frame grid of environmental fields

- `abund_matrix` – site \\\times\\ species table (first column `site`)

- `site_coords` – data.frame with `site, x, y` (quadrat centroids)

When `write_outputs = TRUE`, files are written to disk under
`output_prefix`.

## Files written (when `write_outputs = TRUE`)

- `*_abundances.csv` – site \\\times\\ species matrix

- `*_environments.csv` – mean environment per quadrat

- `*_quadrat_centroids.csv` – centroids table

- `*_fig_panel.png` – 2\\\times\\2 spatial panel

- `*_fig_advanced_panel.png` – advanced diagnostics (if enabled)

- `*_report.txt` – textual analysis report

## See also

[`load_config()`](https://ajsmit.github.io/spesim/reference/load_config.md),
[`load_interactions()`](https://ajsmit.github.io/spesim/reference/load_interactions.md),
[`load_interactions_inline()`](https://ajsmit.github.io/spesim/reference/load_interactions_inline.md),
[`validate_interactions()`](https://ajsmit.github.io/spesim/reference/validate_interactions.md),
[`print_interactions()`](https://ajsmit.github.io/spesim/reference/print_interactions.md),
[`generate_heterogeneous_distribution()`](https://ajsmit.github.io/spesim/reference/generate_heterogeneous_distribution.md),
[`create_abundance_matrix()`](https://ajsmit.github.io/spesim/reference/create_abundance_matrix.md),
[`calculate_quadrat_environment()`](https://ajsmit.github.io/spesim/reference/calculate_quadrat_environment.md),
[`plot_spatial_sampling()`](https://ajsmit.github.io/spesim/reference/plot_spatial_sampling.md),
[`generate_advanced_panel()`](https://ajsmit.github.io/spesim/reference/generate_advanced_panel.md),
[`generate_full_report()`](https://ajsmit.github.io/spesim/reference/generate_full_report.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## A) File-driven
run_spatial_simulation(
  init_file = "inst/examples/spesim_init_complete.txt",
  interactions_file = NULL,
  output_prefix = "out/run"
)

## B) Programmatic (as in README)
dom <- create_sampling_domain()
init <- system.file("examples", "spesim_init_complete.txt", package = "spesim")
P <- load_config(init)
res <- run_spatial_simulation(
  domain = dom, P = P,
  interactions_file = NULL,
  write_outputs = FALSE
)
str(res$abund_matrix)
} # }
```
