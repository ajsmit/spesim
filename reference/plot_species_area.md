# Plot Species-Area Relationship (SAR)

Plots the mean species accumulation (species-area) curve with a ribbon
denoting +/-1 standard deviation across permutations, based on the
summary returned by
[`calculate_species_area`](https://ajsmit.github.io/spesim/reference/calculate_species_area.md).

## Usage

``` r
plot_species_area(sar_data)
```

## Arguments

- sar_data:

  A data frame with columns `Sites` (number of quadrats sampled),
  `Richness` (mean cumulative richness), and `SD` (standard deviation
  across permutations).

## Value

A `ggplot` object.

## Details

The curve represents expected richness as sampling effort increases;
uncertainty is displayed as a semi-transparent ribbon.

## See also

[`calculate_species_area`](https://ajsmit.github.io/spesim/reference/calculate_species_area.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sar <- calculate_species_area(abund_matrix)
plot_species_area(sar)
} # }
```
