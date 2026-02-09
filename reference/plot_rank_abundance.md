# Plot Rank-Abundance Curve

Draws observed and theoretical rank-abundance (SAD) curves on a
log-scaled abundance axis. The input should be the combined table
returned by
[`calculate_rank_abundance()`](https://ajsmit.github.io/spesim/reference/calculate_rank_abundance.md),
which contains one row per rank with a `Source` column distinguishing
*Observed* vs *Theoretical*.

## Usage

``` r
plot_rank_abundance(rank_abundance_data)
```

## Arguments

- rank_abundance_data:

  A data frame with (at minimum) the columns `Rank` (integer rank
  starting at 1), `Abundance` (non-negative numeric counts), and
  `Source` (factor/character with levels such as `"Observed"` and
  `"Theoretical"`).

## Value

A `ggplot` object.

## Details

The y-axis is displayed on a base-10 logarithmic scale. Lines and points
are styled by `Source` (solid vs dashed; different shapes and colors).

## See also

[`calculate_rank_abundance`](https://ajsmit.github.io/spesim/reference/calculate_rank_abundance.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ra <- calculate_rank_abundance(species_dist, P)
plot_rank_abundance(ra)
} # }
```
