# Create an Axis‑Aligned Quadrat Polygon from a Center

Internal helper that builds an axis‑aligned rectangular polygon centered
on center_point with dimensions size = c(width, height).

## Usage

``` r
create_quadrat_from_center(center_point, size)
```

## Arguments

- center_point:

  An sf POINT (sfg or sfc of length 1).

- size:

  numeric length-2, c(width, height).

## Value

An sf polygon geometry (sfg) representing the rectangle.

sfg POLYGON (single rectangle).
