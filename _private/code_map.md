# `spesim` Project Code Map

*Generated: 2025-08-09*

---

## Top-Level Layout

| Path | Purpose |
|------|---------|
| `R/` | Package source code (all exported and internal functions) |
| `man/` | Auto-generated Rd help files via roxygen2 |
| `inst/examples/` | Example initialisation files: `simul_init.txt`, `interactions_init.txt` |
| `tests/` | *testthat* unit tests scaffold |
| `_private/` | Developer notes (`TODO.md`, this file, etc.) |
| other | `DESCRIPTION`, `NAMESPACE`, `LICENSE`, `README.md`, build config files |

---

## Initialisation Files

* `inst/examples/simul_init.txt` – core simulation parameters (species counts, gradients, quadrat settings …)
* `inst/examples/interactions_init.txt` – neighbourhood interaction radius & coefficient matrix/edgelist

---

## Function Index by Source File

### `R/abundance.R`
| Function | Description |
|----------|-------------|
| `generate_fisher_log_series()` | Build rank–abundance vector via Fisher’s log-series with one dominant species. |
| `generate_heterogeneous_distribution()` | Place individuals in the domain using environmental filtering, dominant clustering and local interactions. |

### `R/analysis.R`
| Function | Description |
|----------|-------------|
| `calculate_rank_abundance()` | Compute ranks & abundances table for plotting. |
| `calculate_occupancy_abundance()` | Derive occupancy vs. abundance metrics across quadrats. |
| `calculate_species_area()` | Classic species-area relationship per cumulative quadrats. |
| `calculate_distance_decay()` | Pairwise Bray–Curtis vs. geographic distance. |
| `calculate_rarefaction()` | Rarefaction curve (species richness vs. sites sampled). |
| `plot_rank_abundance()` – `plot_*()` (5) | ggplot helper for each metric above. |

### `R/config.R`
| Function | Description |
|----------|-------------|
| `parse_init_file()` | Low-level parser for `key = value` text configs. |
| `load_config()` | Merge defaults + parsed params; validate and coerce types. |
| `load_interactions()` | Read interaction radius & matrix/edgelist; build defaults if absent. |
| Internal helpers: `.validate_colour()`, `.resolve_param_vector()`, `clamp01()`, etc. |

### `R/domain.R`
| Function | Description |
|----------|-------------|
| `create_sampling_domain()` | Produce irregular polygon representing study area. |

### `R/gradients.R`
| Function | Description |
|----------|-------------|
| `create_environmental_gradients()` | Generate gridded temperature, elevation, rainfall layers with noise. |

### `R/matrices.R`
| Function | Description |
|----------|-------------|
| `create_abundance_matrix()` | Build site × species count table from points & quadrats. |
| `calculate_quadrat_environment()` | Aggregate gradient values per quadrat. |

### `R/panel.R`
| Function | Description |
|----------|-------------|
| `generate_advanced_panel()` | Combine analysis plots into multi-row patchwork panel. |

### `R/plot_spatial_sampling.R`
| Function | Description |
|----------|-------------|
| `plot_spatial_sampling()` | Visualise domain, species points & quadrats; optional gradient overlay. |

### `R/quadrats_random.R`
| Function | Description |
|----------|-------------|
| `place_quadrats()` | Non-overlapping random rectangles within domain. |
| `create_quadrat_from_center()` | Helper constructing sf polygon from centre + size. |

### `R/quadrats_tiled.R`
| Function | Description |
|----------|-------------|
| `place_quadrats_tiled()` | Regular grid fully contained inside domain. |

### `R/quadrats_system.R`
| Function | Description |
|----------|-------------|
| `place_quadrats_systematic()` | Systematic lattice-centred quadrats. |

### `R/quadrats_transect.R`
| Function | Description |
|----------|-------------|
| `place_quadrats_transect()` | Parallel transects with evenly spaced quadrats; supports angle parameter. |

### `R/quadrats_voronoi.R`
| Function | Description |
|----------|-------------|
| `place_quadrats_voronoi()` | Voronoi-based quadrat placement respecting min size. |

### `R/report.R`
| Function | Description |
|----------|-------------|
| `generate_full_report()` | Build plain-text summary of simulation & results. |
| Internal helpers `.opt_to_units()`, `.tol_to_units()`, `.fmt_grad_line()` for formatting. |

### `R/run.R`
| Function | Description |
|----------|-------------|
| `run_spatial_simulation()` | High-level orchestrator: read configs, simulate, sample, analyse & save outputs. |

### `R/utils-internal.R`
Utility infix `%||%`, string/vector parsing helpers, colour validators, etc.

### `R/zzz.R`
Package startup hooks `.onAttach()` and `.onLoad()` – set ggplot theme and global options.

---

## Notes
* Above descriptions derive from roxygen headers and code comments (brevity applied).
* For deeper details, consult each file’s roxygen **@details** section or Rd docs in `man/`.
* Keep this map in sync when adding new functions or files.
