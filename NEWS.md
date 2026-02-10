# spesim 0.3.1

## Usability

- New recommended high-level runner: `spesim_run()` (S3 result object; better teaching workflow defaults).
- New helpers for teaching/workshops: `validate_config()`, `list_examples()`, `spesim_use_example()`, `spesim_open_docs()`.
- New plotting theme: `theme_spesim()`.
- `plot_quadrats()` now supports `show_voronoi = "auto"`/`"yes"`/`"no"` and uses tighter plotting defaults.
- Config parsing errors for `c(...)` vectors now include file + line context.

# spesim 0.3.0

* Add `spesim_demo()` convenience wrapper for teaching/quick exploration.
* Improve pkgdown build cleanliness (vignette titles, NEWS headings).

# spesim 0.2.0-dev

* Implemented Thomas (Rcpp), Strauss (Rcpp fast), and Geyer (Rcpp)
  point processes for clustered species.
