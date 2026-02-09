# Load Simulation Configuration (with defaults & validation)

Reads a text configuration file (KEY = value), merges values with
package defaults, validates options, resolves gradient specifications
into a tidy table, materialises quadrat sizes, and sets the random seed.
This is the **only** supported public entry point for reading configs;
low-level parsing is handled internally.

## Usage

``` r
load_config(init_file)
```

## Arguments

- init_file:

  Character scalar; path to the configuration file.

## Value

A named list `P` with fully-resolved parameters used by the simulator
(including a tidy `P$GRADIENT` tibble and `P$QUADRAT_SIZE`).

## Details

Supported value forms:

- Scalars: `300`, `TRUE`, `0.15`, `"#22223b"`.

- Comma vectors: `A,B,C` or `0.1, 0.2, 0.3`.

- R-style character vectors for inline edgelists: e.g.
  `INTERACTIONS_EDGELIST = c("A,B-D,0.8", "C,A,1.2", "E,*,0.95")`.

- Named entries: `A:0.55`, `temperature:0.12`.

You may specify INTERACTIONS_FILE = "path/to/interactions.txt" to point
to a separate interactions config, or embed rules directly with
`INTERACTIONS_EDGELIST = c("A,B-D,0.8", "C,A,1.2")`. If neither is
provided, interactions default to neutral (all 1.0) with
`INTERACTION_RADIUS = 0`.

Gradient keys supported: `temperature`, `elevation`, `rainfall`. Use
`GRADIENT_SPECIES`, `GRADIENT_ASSIGNMENTS`, plus either
`GRADIENT_OPTIMA` and `GRADIENT_TOLERANCE` as scalar, per-species named,
or per-gradient named values.

## See also

[`run_spatial_simulation`](https://ajsmit.github.io/spesim/reference/run_spatial_simulation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
P <- load_config("inst/examples/spesim_init_complete.txt")
} # }
```
