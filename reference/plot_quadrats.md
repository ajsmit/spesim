# Plot a sampling domain with quadrats

Convenience plotting helper for quadrat-placement workflows. Draws the
sampling domain outline and a set of quadrat polygons (optionally
labelled by quadrat id) using **ggplot2**.

## Usage

``` r
plot_quadrats(
  domain,
  quadrats,
  title = NULL,
  label_ids = TRUE,
  id_col = "quadrat_id",
  domain_colour = "#22223b",
  quadrat_colour = "#22223b",
  quadrat_fill = "#4a4e69",
  quadrat_alpha = 0.18,
  domain_linewidth = 0.7,
  quadrat_linewidth = 0.6,
  label_size = 3
)
```

## Arguments

- domain:

  An sf polygon/multipolygon object representing the sampling domain.

- quadrats:

  An sf polygon/multipolygon object representing quadrats. Typically
  returned by
  [`place_quadrats()`](https://ajsmit.github.io/spesim/reference/place_quadrats.md),
  [`place_quadrats_tiled()`](https://ajsmit.github.io/spesim/reference/place_quadrats_tiled.md),
  [`place_quadrats_systematic()`](https://ajsmit.github.io/spesim/reference/place_quadrats_systematic.md),
  [`place_quadrats_transect()`](https://ajsmit.github.io/spesim/reference/place_quadrats_transect.md),
  or
  [`place_quadrats_voronoi()`](https://ajsmit.github.io/spesim/reference/place_quadrats_voronoi.md).

- title:

  Optional character scalar used as the plot title.

- label_ids:

  Logical. If `TRUE` (default), label each quadrat using the column
  given by `id_col`.

- id_col:

  Character scalar naming the quadrat id column in `quadrats`. Defaults
  to `"quadrat_id"`.

- domain_colour, quadrat_colour:

  Character scalars giving line colours for the domain outline and
  quadrat outlines.

- quadrat_fill:

  Character scalar giving the quadrat fill colour.

- quadrat_alpha:

  Numeric in (0,1\]. Alpha transparency for quadrat fill.

- domain_linewidth, quadrat_linewidth:

  Numeric line widths.

- label_size:

  Numeric text size for quadrat id labels.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
library(spesim)
set.seed(1)
dom <- create_sampling_domain()
qs  <- place_quadrats(dom, n_quadrats = 12, quadrat_size = c(1.5, 1.5))
plot_quadrats(dom, qs, title = "Random quadrats")
} # }
```
