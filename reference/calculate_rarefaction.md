# Rarefaction Curves Data

Generate per-site rarefaction curves (expected richness vs. sample size)
using
[`vegan::rarecurve`](https://vegandevs.github.io/vegan/reference/rarefy.html).

## Usage

``` r
calculate_rarefaction(abund_matrix)
```

## Arguments

- abund_matrix:

  A site x species abundance data frame where the first column is `site`
  and remaining columns are species counts.

## Value

A data frame with columns:

- `SiteID`:

  Factor identifying the site (from `abund_matrix$site`).

- `SampleSize`:

  Number of individuals subsampled.

- `RarefiedRichness`:

  Expected species richness at that sample size.

## Details

The result is returned in long (tidy) format with one row per site x
sample size point. The `SampleSize` values come from the `Subsample`
attribute provided by
[`vegan::rarecurve`](https://vegandevs.github.io/vegan/reference/rarefy.html).

## See also

[`plot_rarefaction`](https://ajsmit.github.io/spesim/reference/plot_rarefaction.md),
[`vegan::rarecurve`](https://vegandevs.github.io/vegan/reference/rarefy.html)

## Examples

``` r
if (FALSE) { # \dontrun{
rr <- calculate_rarefaction(abund_matrix)
head(rr)
} # }
```
