# spesim public API (stability guide)

This document explains which functions in **spesim** are intended for
**user-facing use** (the *public API*), and which helpers may change
more freely.

The short version:

- If it’s **exported** and documented under **Reference** on the pkgdown
  site, it is considered **public**.
- Internal helpers (usually starting with a dot,
  e.g. `.parse_init_file`) are not public and may change without notice.

## Public API groups

### A) Main workflow

These are the primary entry points most users should start with:

- \[load_config()\] — read init files into a parameter list
- \[run_spatial_simulation()\] — end-to-end simulation runner
- \[spesim_demo()\] — one-liner demo/teaching wrapper

### B) Sampling design (quadrats)

These are user-visible so they can be taught and tested directly:

- \[place_quadrats()\] (random, non-overlapping)
- \[place_quadrats_tiled()\]
- \[place_quadrats_systematic()\]
- \[place_quadrats_transect()\]
- \[place_quadrats_voronoi()\]

### C) Plotting helpers

- \[plot_spatial_sampling()\] — map: domain + individuals + quadrats (+
  optional gradient)
- \[plot_quadrats()\] — domain + quadrats (optionally labels; can
  overlay Voronoi cells)

### D) Intermediate analyses (tables) + plots

Tables:

- \[calculate_rank_abundance()\]
- \[calculate_occupancy_abundance()\]
- \[calculate_species_area()\]
- \[calculate_distance_decay()\]
- \[calculate_rarefaction()\]

Plots:

- \[plot_rank_abundance()\]
- \[plot_occupancy_abundance()\]
- \[plot_species_area()\]
- \[plot_distance_decay()\]
- \[plot_rarefaction()\]

### E) Advanced outputs

- \[generate_advanced_panel()\] — multi-plot diagnostic panel
- \[generate_full_report()\] — human-readable text report

## Stability expectations

- **Patch releases** (0.3.x): bug fixes, documentation fixes, small
  additive improvements; avoid breaking the public API.
- **Minor releases** (0.x+1.0): additive features, new arguments, new
  exported functions; behaviour may evolve but should remain broadly
  compatible.
- **Major releases** (1.0+): API cleanup and any intentional breaking
  changes.

If you rely on a function in a teaching module, prefer these patterns:

- Use \[spesim_demo()\] to get a stable “bundle” of outputs.
- When using lower-level functions directly, pin a GitHub commit/tag for
  reproducibility.

## Programmatic discovery

In R you can see exported functions with:

``` r
getNamespaceExports("spesim")
```
