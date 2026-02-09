# Plot Rarefaction Curves

Draws per-site rarefaction curves (expected richness as a function of
the number of individuals sampled) using the long-format output from
[`calculate_rarefaction`](https://ajsmit.github.io/spesim/reference/calculate_rarefaction.md).

## Usage

``` r
plot_rarefaction(rarefaction_data)
```

## Arguments

- rarefaction_data:

  A data frame with columns `SiteID` (factor/character site label),
  `SampleSize` (non-negative integer), and `RarefiedRichness` (expected
  species count).

## Value

A `ggplot` object.

## Details

Each site/quadrat is drawn as a separate line. Colors are mapped to
`SiteID` using a discrete viridis palette for readability.

## See also

[`calculate_rarefaction`](https://ajsmit.github.io/spesim/reference/calculate_rarefaction.md)

## Examples

``` r
if (FALSE) { # \dontrun{
rr <- calculate_rarefaction(abund_matrix)
plot_rarefaction(rr)
} # }
```
