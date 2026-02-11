# spesim 0.3.8

## Filtering diagnostics: plotting helper

- Add `plot_filtering_response()` to visualise environment–abundance relationships for gradient-responsive species.
- The plot shows the configured optimum (dashed line) and tolerance band (shaded) in natural units.
- Environmental gradients vignette updated with a filtering sanity-check section.

# spesim 0.3.7

## Sampling audits: effective sampling frame (area-based)

- Add `estimate_sampling_frame()` and include a conservative boundary constraint summary in `placement_audit`:
  - `safe_area_fraction` = area(safe_domain)/area(domain), where safe_domain is the domain buffered inward by half the quadrat diagonal.
  - `buffer_dist` (half-diagonal) and `quadrat_size` recorded for context.
- Quadrat placement audits and `audit_sampling_scheme()` now report the effective sampling frame fraction alongside exclusion/rejection rates.

# spesim 0.3.6

## Sampling audits: standardise `placement_audit`

- Quadrat placement functions now attach their `placement_audit` attribute via shared helpers and normalise a minimal common schema.
- `audit_sampling_scheme()` now validates/normalises `placement_audit` before reporting.
- Quadrat placement vignette updated with guidance on inspecting `placement_audit` for boundary exclusion.

# spesim 0.3.5

## Method-testing tutor: edge-corrected spatial diagnostics (optional)

- Extend `spesim_audit()` with `diagnostics = c("nn", "spatstat")`.
- Add optional **spatstat-based, edge-corrected** spatial summaries in the *true* irregular window:
  - Ripley-style departure summary via L(r)-r (border correction)
  - pair correlation function extrema (Ripley correction)
- New helper: `spesim_spatstat_diagnostics()` (used internally; requires `spatstat.geom` + `spatstat.explore`).

# spesim 0.3.4

## Method-testing tutor: regime classification

- Add `spesim_regime()` to convert audit metrics into teaching-friendly **green/amber/red** regime classifications.
- `spesim_audit()` now includes a `regime` element and `print.spesim_audit()` prints the classification summary.

# spesim 0.3.3

## Method-testing audits

- Add `spesim_audit()` plus supporting helpers to quantify whether a realised run matches the intended qualitative regime:
  - nearest-neighbour spatial diagnostics per species (clustered/CSR-like/inhibited)
  - environment–abundance checks for gradient-responsive species (optima/tolerances vs realised quadrat environments)
  - sampling scheme boundary exclusion / rejection summaries (when available)
- Quadrat placement functions now attach a `placement_audit` attribute with rejection/exclusion metrics (scheme-dependent).
- The text report now includes a **Conceptual audit** section summarising these checks.

# spesim 0.3.2

## Documentation

- Add **Start here** article to improve pkgdown discoverability.
- Add a **model card** vignette clarifying what spesim does/does not model and how to interpret results.
- Add a **validation & sanity checks** vignette with suggested diagnostics (including optional spatstat checks).
- Clarify that the default sampling domain is synthetic/unitless (users can supply a real-world `sf` domain).

## Reproducibility

- Reset `set.seed(P$SEED)` before generating advanced panels and text reports when `write_outputs = TRUE`, improving determinism for vegan-based components.

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
