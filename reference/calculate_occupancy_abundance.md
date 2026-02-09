# Occupancy-Abundance Table

Summarise total abundance and site occupancy (number of sites with
non-zero counts) for each species from a site x species matrix.

## Usage

``` r
calculate_occupancy_abundance(abund_matrix)
```

## Arguments

- abund_matrix:

  A data frame where the first column is `site` (site / quadrat
  identifier) and the remaining columns are species abundances
  (non-negative integers).

## Value

A data frame with columns:

- `Species`:

  Species (column) name.

- `TotalAbundance`:

  Column sum across sites.

- `Occupancy`:

  Number of sites with abundance \\\> 0\\.

## See also

[`plot_occupancy_abundance`](https://ajsmit.github.io/spesim/reference/plot_occupancy_abundance.md),
[`create_abundance_matrix`](https://ajsmit.github.io/spesim/reference/create_abundance_matrix.md)

## Examples

``` r
if (FALSE) { # \dontrun{
oa <- calculate_occupancy_abundance(abund_matrix)
head(oa)
} # }
```
