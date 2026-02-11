# Changelog

## spesim 0.3.3

### Method-testing audits

- Add
  [`spesim_audit()`](https://ajsmit.github.io/spesim/reference/spesim_audit.md)
  plus supporting helpers to quantify whether a realised run matches the
  intended qualitative regime:
  - nearest-neighbour spatial diagnostics per species
    (clustered/CSR-like/inhibited)
  - environment–abundance checks for gradient-responsive species
    (optima/tolerances vs realised quadrat environments)
  - sampling scheme boundary exclusion / rejection summaries (when
    available)
- Quadrat placement functions now attach a `placement_audit` attribute
  with rejection/exclusion metrics (scheme-dependent).
- The text report now includes a **Conceptual audit** section
  summarising these checks.

## spesim 0.3.2

### Documentation

- Add **Start here** article to improve pkgdown discoverability.
- Add a **model card** vignette clarifying what spesim does/does not
  model and how to interpret results.
- Add a **validation & sanity checks** vignette with suggested
  diagnostics (including optional spatstat checks).
- Clarify that the default sampling domain is synthetic/unitless (users
  can supply a real-world `sf` domain).

### Reproducibility

- Reset `set.seed(P$SEED)` before generating advanced panels and text
  reports when `write_outputs = TRUE`, improving determinism for
  vegan-based components.

## spesim 0.3.1

### Usability

- New recommended high-level runner:
  [`spesim_run()`](https://ajsmit.github.io/spesim/reference/spesim_run.md)
  (S3 result object; better teaching workflow defaults).
- New helpers for teaching/workshops:
  [`validate_config()`](https://ajsmit.github.io/spesim/reference/validate_config.md),
  [`list_examples()`](https://ajsmit.github.io/spesim/reference/list_examples.md),
  [`spesim_use_example()`](https://ajsmit.github.io/spesim/reference/spesim_use_example.md),
  [`spesim_open_docs()`](https://ajsmit.github.io/spesim/reference/spesim_open_docs.md).
- New plotting theme:
  [`theme_spesim()`](https://ajsmit.github.io/spesim/reference/theme_spesim.md).
- [`plot_quadrats()`](https://ajsmit.github.io/spesim/reference/plot_quadrats.md)
  now supports `show_voronoi = "auto"`/`"yes"`/`"no"` and uses tighter
  plotting defaults.
- Config parsing errors for `c(...)` vectors now include file + line
  context.

## spesim 0.3.0

- Add
  [`spesim_demo()`](https://ajsmit.github.io/spesim/reference/spesim_demo.md)
  convenience wrapper for teaching/quick exploration.
- Improve pkgdown build cleanliness (vignette titles, NEWS headings).

## spesim 0.2.0-dev

- Implemented Thomas (Rcpp), Strauss (Rcpp fast), and Geyer (Rcpp) point
  processes for clustered species.
