# Create an irregular sampling domain polygon

Generates a single, closed, irregular polygon intended to represent a
realistic, non-rectangular study area for spatial simulation. The
polygon outline is derived from a parameterized radius-angle curve with
added sinusoidal perturbations and random noise to produce an "organic"
shape.

**Units:** The default domain uses an arbitrary planar coordinate system
(i.e., it is *unitless*). If you need real-world units (m, km) or a
known CRS, supply your own `sf` polygon domain to the
simulation/plotting functions.

The output is returned as an `sf` object containing one polygon in its
geometry column. This can be used directly in plotting or as the spatial
domain within which individuals, quadrats, or transects are placed.

## Usage

``` r
create_sampling_domain()
```

## Value

An `sf` object with one row and a `POLYGON` geometry.

## Details

The shape is generated in polar coordinates \\(r, \theta)\\ using:

- a base radius of 10 units,

- sinusoidal radial modulations with periods of 3 and 5 cycles per
  revolution,

- independent uniform noise on both radius and Cartesian coordinates to
  break symmetry,

- vertical compression by a factor of 0.8 to induce aspect variation.

The 20 vertices are connected in sequence and closed by repeating the
first point. The result is wrapped in
[`sf::st_polygon()`](https://r-spatial.github.io/sf/reference/st.html)
and returned as a single-row `sf` object.

## See also

[`st_polygon`](https://r-spatial.github.io/sf/reference/st.html),
[`st_sf`](https://r-spatial.github.io/sf/reference/sf.html)

## Examples

``` r
if (FALSE) { # \dontrun{
domain <- create_sampling_domain()
plot(domain$geometry)
} # }
```
