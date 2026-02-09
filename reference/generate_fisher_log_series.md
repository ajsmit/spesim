# Generate Fisher's log-series abundances with a dominant species

Creates a rank-abundance vector for *S* species following Fisher's
log-series, while allocating a fixed fraction of individuals to a single
dominant species ("A"). The remaining individuals are distributed across
species B... in decreasing expectation \\\propto \alpha x^{r} / r\\,
where \\r = 2, \ldots, S\\.

## Usage

``` r
generate_fisher_log_series(
  n_species,
  n_individuals,
  dominant_fraction,
  alpha,
  x
)
```

## Arguments

- n_species:

  Integer (\>= 1). Total number of species.

- n_individuals:

  Integer (\>= 1). Total number of individuals to allocate.

- dominant_fraction:

  Numeric in \\\[0,1\]\\. Fraction of individuals assigned to species
  `"A"`; the remainder are split among species `"B"` ... according to
  the log-series.

- alpha:

  Numeric (\> 0). Fisher's \\\alpha\\ parameter controlling tail weight.

- x:

  Numeric, typically close to 1. The log-series scaling parameter.

## Value

A named numeric vector of integer abundances whose names are `"A"`,
`"B"`, `"C"`, ... up to `n_species`. The vector sums to `n_individuals`.

## Details

The algorithm:

1.  Reserves `round(n_individuals * dominant_fraction)` individuals for
    species `"A"`.

2.  Computes relative expectations for ranks `2:n_species` using
    \\\alpha x^{r} / r\\ and scales them to the remaining individuals.

3.  Rounds to integers and adjusts the second species (if present) so
    that the final sum equals `n_individuals`.

4.  Names the vector with `LETTERS[1:n_species]` and drops any species
    whose rounded abundance is `0`.

## Note

Rounding may produce zeros for some tail species; these are omitted from
the return value. If you need a fixed-length vector including zeros,
reindex against `LETTERS[1:n_species]` after the call.

## See also

[`generate_heterogeneous_distribution`](https://ajsmit.github.io/spesim/reference/generate_heterogeneous_distribution.md)

## Examples

``` r
abund <- generate_fisher_log_series(
  n_species = 10, n_individuals = 1000,
  dominant_fraction = 0.3, alpha = 3, x = 0.95
)
sum(abund) # 1000
#> [1] 1000
```
