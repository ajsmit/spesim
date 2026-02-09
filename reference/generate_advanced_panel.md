# Advanced Analysis Panel (rank/occupancy/SAR/decay/rarefaction + optional PPC)

Build a multi‑panel summary of key ecological diagnostics from a
simulation run: rank–abundance (SAD), occupancy–abundance, species–area
(accumulation), distance–decay, and per‑quadrat rarefaction curves. If
the point‑process packages spatstat.geom and spatstat.explore are
available, an additional row of **point–process diagnostics** is
appended, showing Ripley’s K (border correction), \\L(r)-r\\, and the
pair‑correlation function \\g(r)\\.

## Usage

``` r
generate_advanced_panel(res)
```

## Arguments

- res:

  A list produced by the main simulator containing at least:

  `P`

  :   Parameter list (used for labels/themes if needed).

  `species_dist`

  :   An `sf` POINT layer of individuals with a `species` column.

  `abund_matrix`

  :   Site \\\times\\ species abundance data frame (first column
      `site`).

  `site_coords`

  :   Data frame with `site, x, y` for quadrat centroids.

  `domain`

  :   An `sf` polygon/multipolygon of the study area (used for PPC
      window).

## Value

A `patchwork` object (a `ggplot` layout). You can print it or save it
with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Details

The optional point–process diagnostics are computed only when both
spatstat.geom (for spatial windows and point patterns) and
spatstat.explore (for summary functions such as `Kest` and `pcf`) are
installed. The domain is converted to a window (`owin`) using
[`spatstat.geom::as.owin()`](https://rdrr.io/pkg/spatstat.geom/man/as.owin.html)
on the `sf` polygon; the individuals are converted to a `ppp` with that
window. We then compute:

- `Kest(ppp, correction = "border")` and plot the border‑corrected K.

- \\L(r) = \sqrt{K(r)/\pi}\\ and \\L(r)-r\\ (\\\>0\\ suggests
  clustering).

- `pcf(ppp)` as an estimate of the pair‑correlation \\g(r)\\ (\\\>1\\
  suggests clustering; \\\<1\\ inhibition).

If conversion to `owin`/`ppp` fails, the PPC row is omitted gracefully.

## See also

[`calculate_rank_abundance`](https://ajsmit.github.io/spesim/reference/calculate_rank_abundance.md),
[`plot_rank_abundance`](https://ajsmit.github.io/spesim/reference/plot_rank_abundance.md),
[`calculate_occupancy_abundance`](https://ajsmit.github.io/spesim/reference/calculate_occupancy_abundance.md),
[`plot_occupancy_abundance`](https://ajsmit.github.io/spesim/reference/plot_occupancy_abundance.md),
[`calculate_species_area`](https://ajsmit.github.io/spesim/reference/calculate_species_area.md),
[`plot_species_area`](https://ajsmit.github.io/spesim/reference/plot_species_area.md),
[`calculate_distance_decay`](https://ajsmit.github.io/spesim/reference/calculate_distance_decay.md),
[`plot_distance_decay`](https://ajsmit.github.io/spesim/reference/plot_distance_decay.md),
[`calculate_rarefaction`](https://ajsmit.github.io/spesim/reference/calculate_rarefaction.md),
[`plot_rarefaction`](https://ajsmit.github.io/spesim/reference/plot_rarefaction.md)

## Examples

``` r
if (FALSE) { # \dontrun{
panel <- generate_advanced_panel(results_list)
ggplot2::ggsave("advanced_panel.png", panel, width = 12, height = 14, dpi = 300)
} # }
```
