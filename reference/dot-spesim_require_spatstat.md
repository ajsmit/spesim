# Optional spatstat-based spatial diagnostics (edge-corrected)

These helpers are used by
[`spesim_audit()`](https://ajsmit.github.io/spesim/reference/spesim_audit.md)
when `diagnostics` includes `"spatstat"`. They are kept optional to
avoid a heavy dependency chain.

The aim is not to replace full point-process analysis, but to provide
*edge-corrected* summaries (in the correct irregular window) that help
users verify that the realised spatial regime matches the named process
they think they asked for.

## Usage

``` r
.spesim_require_spatstat()
```
