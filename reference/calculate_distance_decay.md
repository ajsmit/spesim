# Distance-Decay Data

Pair geographic distances between sites with community dissimilarities
(Sorensen index computed as binary Bray-Curtis) for distance-decay
plots.

## Usage

``` r
calculate_distance_decay(abund_matrix, site_coords)
```

## Arguments

- abund_matrix:

  A site x species abundance data frame where the first column is `site`
  and the remaining columns are abundances.

- site_coords:

  A data frame with numeric columns `x` and `y` giving site coordinates
  (in a projected CRS) in the same row order as `abund_matrix`.

## Value

A data frame with two numeric columns:

- `Distance`:

  Euclidean distance between site pairs.

- `Dissimilarity`:

  Sorensen dissimilarity (0-1).

## See also

[`plot_distance_decay`](https://ajsmit.github.io/spesim/reference/plot_distance_decay.md),
[`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)

## Examples

``` r
if (FALSE) { # \dontrun{
dd <- calculate_distance_decay(abund_matrix, site_coords)
head(dd)
} # }
```
