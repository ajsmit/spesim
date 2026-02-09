# Species-Area (Accumulation) Data

Compute a species accumulation curve (SAR) using
`vegan::specaccum(method = "random")`.

## Usage

``` r
calculate_species_area(abund_matrix)
```

## Arguments

- abund_matrix:

  A site x species abundance data frame where the first column is `site`
  and remaining columns are species abundances.

## Value

A data frame with columns:

- `Sites`:

  Number of sites sampled.

- `Richness`:

  Mean cumulative species richness.

- `SD`:

  Standard deviation across permutations.

## Details

Uses 100 random permutations by default to estimate mean richness and
its standard deviation as a function of the number of sites.

## See also

[`plot_species_area`](https://ajsmit.github.io/spesim/reference/plot_species_area.md),
[`vegan::specaccum`](https://vegandevs.github.io/vegan/reference/specaccum.html)

## Examples

``` r
if (FALSE) { # \dontrun{
sar <- calculate_species_area(abund_matrix)
head(sar)
} # }
```
