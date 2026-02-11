# Plot environmental filtering responses (sanity check)

A teaching/method-testing helper that visualises whether a
gradient-responsive species shows higher abundance in quadrats whose
environments are near the configured optimum (with tolerance band).

This plot is designed to pair with
[`audit_environmental_filtering()`](https://ajsmit.github.io/spesim/reference/audit_environmental_filtering.md).
It is intentionally simple: it uses per-quadrat mean environments (as
computed by
[`calculate_quadrat_environment()`](https://ajsmit.github.io/spesim/reference/calculate_quadrat_environment.md))
and the site x species abundance matrix.

## Usage

``` r
plot_filtering_response(res, species = "all", add_smooth = TRUE)
```

## Arguments

- res:

  A `spesim_result` list returned by
  [`spesim_run()`](https://ajsmit.github.io/spesim/reference/spesim_run.md).

- species:

  Species labels to plot. Default `"all"` uses all gradient-responsive
  species in `res$P$GRADIENT`.

- add_smooth:

  Logical; add a loess smooth (teaching-friendly). Default `TRUE`.

## Value

A ggplot object.
