# Parse inline interaction rules from the init file

Converts `INTERACTIONS_EDGELIST` (CSV-like string; optional
`INTERACTION_RADIUS`) into the resolved interaction list (`radius`,
`matrix`).

## Usage

``` r
load_interactions_inline(rules, n_species, radius = 0)
```

## Arguments

- rules:

  Character scalar or character vector of CSV-like lines with columns
  `focal,neighbour,value` and optional wildcard/range syntax (see
  vignette).

- n_species:

  Integer S; used to define the full SxS matrix.

- radius:

  Optional numeric radius override; if `NULL`, defaults to 0.

## Value

A list with `radius` (numeric) and `matrix` (SxS numeric with dimnames).

## Examples

``` r
if (FALSE) { # \dontrun{
load_interactions_inline(
  rules = c("A,B-D,0.8", "C,A,1.2", "E,*,0.95"),
  radius = 2, n_species = 10
)
} # }
```
