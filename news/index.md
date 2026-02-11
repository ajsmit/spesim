# Changelog

## spesim 0.3.14

### Build hygiene and vignette snippet checks

- Removed pkgdown/vignette build warnings by aligning vignette titles
  and setting `highlight: null` for html vignettes.
- Added a test to ensure all vignette code chunks are syntactically
  valid ([`knitr::purl()`](https://rdrr.io/pkg/knitr/man/knit.html) +
  [`parse()`](https://rdrr.io/r/base/parse.html)).

## spesim 0.3.13

### Performance notes

- New vignette: “Performance notes” (`spesim-performance`) with
  conservative interactive ranges and practical tips.

## spesim 0.3.12

### Method-testing vignette

- New vignette: “Method testing with spesim” (`spesim-method-testing`).

## spesim 0.3.11

### Method-testing bundle helper

- New
  [`spesim_method_test()`](https://ajsmit.github.io/spesim/reference/spesim_method_test.md)
  convenience function returns `res`, `audit`, and optional standard
  tables/plots for teaching and method-testing.

## spesim 0.3.10

### Report: configurable conceptual audit

- [`generate_full_report()`](https://ajsmit.github.io/spesim/reference/generate_full_report.md)
  gains `include_audit` and `audit_top_n` to control the report’s
  conceptual audit section.
- Report audit now runs `spesim_audit(..., diagnostics = "nn")` for
  deterministic, lightweight behaviour.

## spesim 0.3.9

### Filtering audits: occupancy-aware diagnostics

- [`audit_environmental_filtering()`](https://ajsmit.github.io/spesim/reference/audit_environmental_filtering.md)
  now reports `occupancy_sites` (number of quadrats with abundance \>
  0).
- Low occupancy is flagged as `too_sparse_to_assess`, helping prevent
  over-interpretation of noisy correlations.

## spesim 0.3.8

### Filtering diagnostics: plotting helper

- Add
  [`plot_filtering_response()`](https://ajsmit.github.io/spesim/reference/plot_filtering_response.md)
  to visualise environment–abundance relationships for
  gradient-responsive species.
- The plot shows the configured optimum (dashed line) and tolerance band
  (shaded) in natural units.
- Environmental gradients vignette updated with a filtering sanity-check
  section.

## spesim 0.3.7

### Sampling audits: effective sampling frame (area-based)

- Add
  [`estimate_sampling_frame()`](https://ajsmit.github.io/spesim/reference/estimate_sampling_frame.md)
  and include a conservative boundary constraint summary in
  `placement_audit`:
  - `safe_area_fraction` = area(safe_domain)/area(domain), where
    safe_domain is the domain buffered inward by half the quadrat
    diagonal.
  - `buffer_dist` (half-diagonal) and `quadrat_size` recorded for
    context.
- Quadrat placement audits and
  [`audit_sampling_scheme()`](https://ajsmit.github.io/spesim/reference/audit_sampling_scheme.md)
  now report the effective sampling frame fraction alongside
  exclusion/rejection rates.

## spesim 0.3.6

### Sampling audits: standardise `placement_audit`

- Quadrat placement functions now attach their `placement_audit`
  attribute via shared helpers and normalise a minimal common schema.
- [`audit_sampling_scheme()`](https://ajsmit.github.io/spesim/reference/audit_sampling_scheme.md)
  now validates/normalises `placement_audit` before reporting.
- Quadrat placement vignette updated with guidance on inspecting
  `placement_audit` for boundary exclusion.

## spesim 0.3.5

### Method-testing tutor: edge-corrected spatial diagnostics (optional)

- Extend
  [`spesim_audit()`](https://ajsmit.github.io/spesim/reference/spesim_audit.md)
  with `diagnostics = c("nn", "spatstat")`.
- Add optional **spatstat-based, edge-corrected** spatial summaries in
  the *true* irregular window:
  - Ripley-style departure summary via L(r)-r (border correction)
  - pair correlation function extrema (Ripley correction)
- New helper:
  [`spesim_spatstat_diagnostics()`](https://ajsmit.github.io/spesim/reference/spesim_spatstat_diagnostics.md)
  (used internally; requires `spatstat.geom` + `spatstat.explore`).

## spesim 0.3.4

### Method-testing tutor: regime classification

- Add
  [`spesim_regime()`](https://ajsmit.github.io/spesim/reference/spesim_regime.md)
  to convert audit metrics into teaching-friendly **green/amber/red**
  regime classifications.
- [`spesim_audit()`](https://ajsmit.github.io/spesim/reference/spesim_audit.md)
  now includes a `regime` element and `print.spesim_audit()` prints the
  classification summary.

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
