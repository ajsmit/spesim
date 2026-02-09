# Create synthetic environmental gradients over a spatial domain

Generates three continuous, spatially autocorrelated environmental
gradients — temperature, elevation, and rainfall — over a regular grid
covering the specified polygonal study area. The gradients are expressed
in both normalised \\\[0, 1\]\\ form and scaled to physical units.

## Usage

``` r
create_environmental_gradients(domain, resolution, noise_level)
```

## Arguments

- domain:

  An `sf` polygon object defining the extent over which gradients will
  be computed. Only its bounding box is used in this step.

- resolution:

  Integer ≥ 2; number of grid steps per axis. The output will contain
  `resolution^2` rows.

- noise_level:

  Non-negative numeric; standard deviation of Gaussian noise added
  independently to each gradient before clipping.

## Value

A `data.frame` with columns:

- `x, y`:

  Grid point coordinates (same CRS as `domain`).

- `temperature, elevation, rainfall`:

  Normalised gradients in \\\[0, 1\]\\.

- `temperature_C`:

  Temperature in degrees Celsius (approx. range: -2 to 28).

- `elevation_m`:

  Elevation in metres (0–2000).

- `rainfall_mm`:

  Annual rainfall in millimetres (200–900).

## Details

The function:

1.  Derives the bounding box of `domain` and constructs a regular
    `resolution` × `resolution` grid of points.

2.  Normalises \\x\\ and \\y\\ coordinates to \\\[0, 1\]\\.

3.  Computes:

    - **Temperature**: weighted linear combination of normalised \\x\\
      (0.7) and \\y\\ (0.3) coordinates.

    - **Elevation**: decreases radially from domain centre.

    - **Rainfall**: weighted difference of \\x\\ and \\y\\ components,
      rescaled to \\\[0, 1\]\\.

4.  Adds Gaussian noise with standard deviation `noise_level` to each
    raw gradient before clipping to \\\[0, 1\]\\.

5.  Produces scaled versions in common physical units: `temperature_C`
    (°C), `elevation_m` (metres), and `rainfall_mm` (millimetres/year).

Note that gradients are computed for all grid points within the
*bounding box* of `domain`; no masking to the polygon interior is
applied here. Masking can be done later with
[`st_intersects`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html)
or similar.

## See also

[`st_bbox`](https://r-spatial.github.io/sf/reference/st_bbox.html),
[`st_intersects`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html)

## Examples

``` r
if (FALSE) { # \dontrun{
domain <- create_sampling_domain()
env <- create_environmental_gradients(domain, resolution = 50, noise_level = 0.05)
head(env)
} # }
```
