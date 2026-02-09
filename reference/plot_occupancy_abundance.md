# Plot Occupancy-Abundance Relationship

Scatterplot of total abundance (per species) versus site occupancy
(number of quadrats in which the species occurs) with both axes on a log
scale. Adds an optional least-squares trend line for visual guidance.

## Usage

``` r
plot_occupancy_abundance(oa_data)
```

## Arguments

- oa_data:

  A data frame as returned by
  [`calculate_occupancy_abundance`](https://ajsmit.github.io/spesim/reference/calculate_occupancy_abundance.md)
  with columns `Species`, `TotalAbundance` (non-negative numeric), and
  `Occupancy` (non-negative integer).

## Value

A `ggplot` object.

## Details

If the input contains no observations (or all totals are zero), a
minimal placeholder plot is returned indicating that there is nothing to
draw. Otherwise both axes are shown on base-10 logarithmic scales.

## See also

[`calculate_occupancy_abundance`](https://ajsmit.github.io/spesim/reference/calculate_occupancy_abundance.md)

## Examples

``` r
if (FALSE) { # \dontrun{
oa <- calculate_occupancy_abundance(abund_matrix)
plot_occupancy_abundance(oa)
} # }
```
