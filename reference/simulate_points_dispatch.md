# Dispatch simulation by process kind

A small internal router that calls one of the simulators:

- `"poisson"` — homogeneous Poisson (uniform) via
  [`sf::st_sample()`](https://r-spatial.github.io/sf/reference/st_sample.html)

- `"thomas"` — parent–offspring clustering

- `"strauss"` — soft inhibition (sequential acceptance)

- `"geyer"` — saturation model (sequential acceptance)

## Usage

``` r
simulate_points_dispatch(kind, domain, n_target, args = list())
```

## Arguments

- kind:

  Character scalar: one of `"poisson"`, `"thomas"`, `"strauss"`,
  `"geyer"` (case-insensitive).

- domain:

  sf polygon/multipolygon sampling window.

- n_target:

  Integer target number of points.

- args:

  Named list of extra parameters, as described above.

## Value

sf POINT layer (no attributes).

## Details

Arguments for a given process can be supplied via `args` using either
the generic names (`beta`, `mu`, `sigma`, `gamma`, `r`, `sat`) or the
package’s config-style names (`A_PARENT_INTENSITY`, `A_MEAN_OFFSPRING`,
`A_CLUSTER_SCALE`, `OTHERS_BETA`, `OTHERS_GAMMA`, `OTHERS_R`,
`OTHERS_S`).
