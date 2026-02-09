# Rank-Abundance Data

Construct a tidy table for rank-abundance plotting by combining the
**observed** species counts in `species_dist` with the **theoretical**
Fisher log-series abundances implied by `P`.

## Usage

``` r
calculate_rank_abundance(species_dist, P)
```

## Arguments

- species_dist:

  An `sf` point object with a character column `species` giving
  per-individual species identities.

- P:

  A named list as returned by
  [`load_config`](https://ajsmit.github.io/spesim/reference/load_config.md)
  containing at least the fields used by
  [`generate_fisher_log_series`](https://ajsmit.github.io/spesim/reference/generate_fisher_log_series.md):
  `N_SPECIES`, `N_INDIVIDUALS`, `DOMINANT_FRACTION`, `FISHER_ALPHA`, and
  `FISHER_X`.

## Value

A data frame with columns:

- `Rank`:

  Integer species rank (1 = most abundant).

- `Abundance`:

  Non-negative integer abundance.

- `Source`:

  `"Observed"` or `"Theoretical"`.

## Details

Observed abundances are computed with base R
[`table()`](https://rdrr.io/r/base/table.html), ranked in descending
order, and annotated as `Source = "Observed"`. Theoretical abundances
are generated via
[`generate_fisher_log_series`](https://ajsmit.github.io/spesim/reference/generate_fisher_log_series.md)
using parameters in `P` (e.g., `N_SPECIES`, `N_INDIVIDUALS`,
`FISHER_ALPHA`, `FISHER_X`, `DOMINANT_FRACTION`) and annotated as
`Source = "Theoretical"`. The two data frames are row-bound.

## See also

[`plot_rank_abundance`](https://ajsmit.github.io/spesim/reference/plot_rank_abundance.md),
[`generate_fisher_log_series`](https://ajsmit.github.io/spesim/reference/generate_fisher_log_series.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ra <- calculate_rank_abundance(species_dist, P)
head(ra)
} # }
```
