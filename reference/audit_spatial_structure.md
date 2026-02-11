# Audit spatial structure via nearest-neighbour distances

Audit spatial structure via nearest-neighbour distances

## Usage

``` r
audit_spatial_structure(species_dist, domain, species = "all", nn_k = 1)
```

## Arguments

- species_dist:

  `sf` points with a `species` column.

- domain:

  `sf` polygon used only for area estimation.

- species:

  Species labels to include; `"all"` uses all present.

- nn_k:

  neighbour order (1 = nearest neighbour).

## Value

data.frame with NN summary per species.
