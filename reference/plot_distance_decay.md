# Plot Distance-Decay Relationship

Displays community dissimilarity versus geographic distance with a
smooth loess trend. The y-axis is constrained to \[0, 1\] to reflect the
range of Sorensen (binary Bray-Curtis) dissimilarity.

## Usage

``` r
plot_distance_decay(decay_data)
```

## Arguments

- decay_data:

  A data frame as returned by
  [`calculate_distance_decay`](https://ajsmit.github.io/spesim/reference/calculate_distance_decay.md)
  with columns `Distance` (pairwise Euclidean distances among quadrat
  centroids) and `Dissimilarity` (pairwise Sorensen dissimilarities).

## Value

A `ggplot` object.

## Details

Points show pairwise site comparisons; the loess smoother summarizes the
overall distance-decay pattern. Input should already exclude
self-comparisons.

## See also

[`calculate_distance_decay`](https://ajsmit.github.io/spesim/reference/calculate_distance_decay.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dd <- calculate_distance_decay(abund_matrix, site_coords)
plot_distance_decay(dd)
} # }
```
