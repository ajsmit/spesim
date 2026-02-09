# Plot the spatial sampling simulation (domain, individuals, quadrats, optional gradient)

Creates a publication‑quality map of the simulated landscape. The domain
outline is drawn first; simulated individuals (POINT geometries) are
overlaid and colored by species; sampling quadrats (POLYGONs) are drawn
with IDs. Optionally an environmental raster/heat field is shown
underneath for context.

## Usage

``` r
plot_spatial_sampling(
  domain,
  species,
  quadrats,
  P,
  show_gradient = FALSE,
  env_gradients = NULL,
  gradient_type = "temperature_C"
)
```

## Arguments

- domain:

  An `sf` polygon or multipolygon representing the study area. Geometry
  should be valid and in the same CRS as `species` and `quadrats`.

- species:

  An `sf` POINT object of individuals with a character/factor column
  named `species`. All points must lie within `domain` for correct
  visual interpretation (the function does not clip).

- quadrats:

  An `sf` POLYGON/MULTIPOLYGON object with a numeric/integer column
  `quadrat_id`. Quadrats are drawn as outlines and labeled at their
  centroids.

- P:

  A named list of plotting/summary parameters. The following fields are
  read if present (each has a default if missing):

  - `POINT_SIZE` (`numeric`, default `0.2`) — point size for
    individuals,

  - `POINT_ALPHA` (`numeric` in \\(0,1)\\, default `1`) — point
    transparency,

  - `QUADRAT_COLOUR` (`character`, default `“black”`) — quadrat
    outline/label colour,

  - `BACKGROUND_COLOUR` (`character`, default `“white”`) — plot
    background,

  - `FOREGROUND_COLOUR` (`character`, default `”#22223b”`) — domain
    outline/title colour,

  - `N_SPECIES` (`integer`) — used to size the palette; inferred from
    data if absent,

  - `N_INDIVIDUALS` (`integer`) — used only for the subtitle count if
    available.

- show_gradient:

  Logical; if `TRUE`, an environmental surface is drawn beneath the
  geometries using `geom_raster()`. Default `FALSE`.

- env_gradients:

  *Required when* `show_gradient = TRUE`. A data frame with columns `x`,
  `y` (grid coordinates, in the same CRS units as `domain`) and one or
  more numeric columns representing environmental variables (e.g.,
  `temperature_C`, `elevation_m`, `rainfall_mm`). The column named by
  `gradient_type` must exist and be numeric.

- gradient_type:

  Character scalar naming the column in `env_gradients` to plot when
  `show_gradient = TRUE`. Default `“temperature_C”`. The legend title is
  derived from this name (underscores replaced with spaces and
  title‑cased).

## Value

A `ggplot` object that can be further modified (themes, scales, etc.).

## Details

Colors for species are drawn from colorspace `sequential_hcl` palette
“RdPu” (reversed) and mapped to the unique species present. All layers
are plotted in the current display CRS of the provided `sf` objects;
ensure consistent CRS across inputs. When `show_gradient = TRUE`, the
raster is drawn first with partial transparency so the domain outline
and points remain visible.

## See also

[`run_spatial_simulation()`](https://ajsmit.github.io/spesim/reference/run_spatial_simulation.md),
[`create_sampling_domain()`](https://ajsmit.github.io/spesim/reference/create_sampling_domain.md),
[`create_environmental_gradients()`](https://ajsmit.github.io/spesim/reference/create_environmental_gradients.md)

## Examples
