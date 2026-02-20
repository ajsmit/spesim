# spesim 0.5.1

## Performance

- **C++ point-in-polygon engine.** All polygon domain checks previously used
  `sf::st_as_sf()` + `sf::st_within()` on every call, incurring full CRS
  comparison and GEOS setup overhead. They now delegate to a compiled C++
  ray-casting implementation (`pip_cpp`), with ring coordinates extracted once
  and cached for the duration of any tight loop. Measured speedups at
  representative call sizes: **6–13× faster** for batches of 1–40 points,
  **112× faster** for the sequential-inhibition (Strauss) single-point loop.

- **Vectorised initial community placement.** `spesim_simulate_neutral_recruitment()`
  previously sampled one position per individual inside a `for` loop
  (J individual calls to `sf::st_within`). Placement is now done in a single
  vectorised call to `.sample_uniform_in_domain()`, eliminating J − 1 polygon
  coordinate extractions and sf round-trips.

- **Vectorised sf point construction.** `.as_sfc_points()` replaced a
  `lapply(sf::st_point(...))` loop with a single `sf::st_as_sf()` call,
  which is the recommended vectorised path.

- **Shared pip tester across all helpers.** `.sample_uniform_in_domain()`,
  `.simulate_strauss_points()`, and `.linear_pos_to_xy()` all accept an
  optional `pip_fn` argument so callers that already hold a pre-built tester
  avoid redundant ring-coordinate extraction.

# spesim 0.4.6

## Maintenance

- Added `rnaturalearthdata` to `Suggests` to satisfy CRAN dependency declaration (used in vignette via `requireNamespace()`).

# spesim 0.4.5

## Bug fixes + promised features now active

- Fixed `generate_fisher_log_series()` for `n_species = 1` (now sums to `n_individuals`).
- Fixed fast Strauss MH sampler state updates (accept/reject now consistent).
- Fixed Thomas and Strauss/Geyer fast helpers to be polygon-aware (bbox simulate + clip + thin/top-up).
- Implemented environmental filtering + local interspecific interactions during species assignment as documented.
- Distance–decay now aligns `site_coords` to `abund_matrix` by site id and uses a robust Sorensen dissimilarity for presence/absence.
- `ZSM_M` is now active: when supplied as a finite value in `(0, 1)`, `SAD_MODEL = "zsm"` switches from the theta-only Ewens (Chinese restaurant) sampler to a Moran-style death–birth with immigration heuristic; lower `m` produces a more uneven SAD. The NEWS 0.4.2 entry below stating `ZSM_M` is "ignored" is superseded by this release.

# spesim 0.4.2

## Vignettes / neutral SAD clarification

- Updated point-process recipes vignette with an explicit explanation of the current neutral-theory SAD helper (what it can and cannot do), plus worked examples and advanced panels.
- `SAD_MODEL = "zsm"` is now documented and implemented as a **theta-only Ewens sampler**; `ZSM_M` is currently ignored (accepted for forward compatibility).

# spesim 0.4.1

## Maintenance

- Build/check hygiene improvements (CRAN incoming note reductions).

# spesim 0.4.0

## Runner refactor + stable results

- `spesim_run()` is now the recommended, consolidated entry point and returns a stable S3 results object (`spesim_result`).
- `run_spatial_simulation()` is retained as a **deprecated** wrapper around `spesim_run()` for backwards compatibility.
- Vignettes and README updated to reflect the new workflow.

# spesim 0.3.19

## Startup message

- `library(spesim)` now prints the package version in the startup message.

# spesim 0.3.18

## Concave domain complexity control

- `create_sampling_domain(shape = "concave")` gains `concave_lobes` (2–5) to control large-scale outline complexity; lower values yield simpler concave polygons.

# spesim 0.3.17

## Default domain shape diversity

- `create_sampling_domain()` now supports multiple shape families (concave/irregular, convex hulls, ellipses/circles, and boxy rectangles/squares), and the default uses a mixed generator to better mimic real-world landscape patch shapes.

# spesim 0.3.16

## Default domain variability in spesim_run()

- `spesim_run()` now generates a new random default domain for each run when `domain = NULL` and no `seed` is supplied.
- When `seed` is supplied, the default domain is reproducible across runs.

## Plotting side-effects

- `calculate_rarefaction()` no longer triggers any base-graphics plotting side-effects (fixes duplicate rarefaction output when used inside `generate_advanced_panel()`).

# spesim 0.3.15

## Ordination teaching vignette (db-RDA)

- New vignette: constrained ordination on simulated data using distance-based RDA (db-RDA) via `vegan::capscale()`.

# spesim 0.3.14

## Build hygiene and vignette snippet checks

- Removed pkgdown/vignette build warnings by aligning vignette titles and setting `highlight: null` for html vignettes.
- Added a test to ensure all vignette code chunks are syntactically valid (`knitr::purl()` + `parse()`).

# spesim 0.3.13

## Performance notes

- New vignette: "Performance notes" (`spesim-performance`) with conservative interactive ranges and practical tips.

# spesim 0.3.12

## Method-testing vignette

- New vignette: "Method testing with spesim" (`spesim-method-testing`).

# spesim 0.3.11

## Method-testing bundle helper

- New `spesim_method_test()` convenience function returns `res`, `audit`, and optional standard tables/plots for teaching and method-testing.

# spesim 0.3.10

## Report: configurable conceptual audit

- `generate_full_report()` gains `include_audit` and `audit_top_n` to control the report's conceptual audit section.
- Report audit now runs `spesim_audit(..., diagnostics = "nn")` for deterministic, lightweight behaviour.

# spesim 0.3.9

## Filtering audits: occupancy-aware diagnostics

- `audit_environmental_filtering()` now reports `occupancy_sites` (number of quadrats with abundance > 0).
- Low occupancy is flagged as `too_sparse_to_assess`, helping prevent over-interpretation of noisy correlations.

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
