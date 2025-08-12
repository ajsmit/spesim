spesim: Spatial Ecological Simulation in R
================
AJ Smit and contributors
2025-08-12

# spesim <img src="man/figures/logo.png" align="right" width="120" />

**spesim** is an R package for simulating, sampling, and analysing
*spatially heterogeneous ecological communities* in irregular
landscapes.  
It’s built for teaching, methods testing, and exploratory research in
**biogeography** and **community ecology**.

- Generate individuals for multiple species under realistic
  **species–abundance distributions** (e.g., Fisher’s log-series with a
  dominant species).
- Impose **environmental filtering** (per‑species optima/tolerances on
  named gradients).
- Add **spatial structure** (dominant clustering + optional
  interspecific neighbourhood effects).
- Sample with flexible **quadrat schemes** (random, systematic, tiled,
  transect, Voronoi).
- Produce **site × species** matrices, per‑site **environment
  summaries**, and **diversity partitions** (α, β, γ).
- Save publication‑quality **maps** and an optional **advanced analysis
  panel** (rank–abundance, occupancy–abundance, species–area,
  distance–decay, rarefaction).
- Write a human‑readable **text report** that explains what happened.

## Why spesim?

- **Teaching**: make abstract ideas (environmental filtering, distance
  decay, β‑diversity) tangible with quick simulations you can show in
  class.
- **Methods development**: prototype sampling designs (how many
  quadrats? which scheme?) and test metrics before going to the field.
- **Exploration**: stress‑test analyses on known “truth” to see how
  robust your inferences are.

------------------------------------------------------------------------

## Installation

Development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("ajsmit/spesim")
```

No heavy *spatstat* dependency is required for core functionality. The
package relies on common spatial/data/plotting libraries (see
`DESCRIPTION` for imports).

> Building the optional **pkgdown** site locally? Use:
>
> ``` r
> # install.packages("pkgdown")
> pkgdown::build_site()
> ```

------------------------------------------------------------------------

## Vignettes (online)

- [**Interactions (neighbour
  effects)**](https://ajsmit.github.io/spesim/articles/spesim-interactions.html)
- [**Full init-file
  parameters**](https://ajsmit.github.io/spesim/articles/spesim-init-parameters.html)
- [**Environmental
  Gradients**](https://ajsmit.github.io/spesim/articles/spesim-env-gradients.html)
- [**Advanced Analysis
  Panel**](https://ajsmit.github.io/spesim/articles/spesim-advanced-panel.html)
- [**Point
  processes**](https://ajsmit.github.io/spesim/articles/spesim-point-processes.html)

If reading on GitHub before the site is live, you’ll find the sources
under `vignettes/`.

------------------------------------------------------------------------

## Quick start

You can run a simulation in two ways: **programmatically** (in-memory)
or using an **init file** on disk. The latter is more convenient for
larger setups and reproducibility, since access to the init file gives
you a full bird’s eye view of the simulation parameters.

### Quick start (programmatic)

A minimal run that *does not* write files (handy for examples/CI).

``` r
library(spesim)

# Load a complete example init, then tweak a few things in-memory
P <- load_config(system.file("examples/spesim_init_basic.txt", package = "spesim"))
P$N_SPECIES <- 10
P$N_INDIVIDUALS <- 2000
P$SAMPLING_SCHEME <- "random"
P$N_QUADRATS <- 20
P$QUADRAT_SIZE_OPTION <- "medium"

# Run (skip writing to disk to keep examples snappy)
res <- run_spatial_simulation(P = P, write_outputs = FALSE)

# One-map view: individuals + quadrats (no gradient fill)
p <- plot_spatial_sampling(res$domain, res$species_dist, res$quadrats, res$P)
print(p)

# Optional: a 2×2 panel including gradient overlays
p1 <- plot_spatial_sampling(res$domain, res$species_dist, res$quadrats, res$P)
p2 <- plot_spatial_sampling(res$domain, res$species_dist, res$quadrats, res$P,
                            show_gradient = TRUE, env_gradients = res$env_gradients, gradient_type = "temperature_C")
p3 <- plot_spatial_sampling(res$domain, res$species_dist, res$quadrats, res$P,
                            show_gradient = TRUE, env_gradients = res$env_gradients, gradient_type = "elevation_m")
p4 <- plot_spatial_sampling(res$domain, res$species_dist, res$quadrats, res$P,
                            show_gradient = TRUE, env_gradients = res$env_gradients, gradient_type = "rainfall_mm")
(p1 | p2) / (p3 | p4)
```

### Quick start (init file on disk)

Prefer declaring everything in a text file? Point
`run_spatial_simulation()` to it:

``` r
library(spesim)

# Your own config on disk
init <- "inst/examples/spesim_init_complete.txt"   # or a path you created

# Optional: interactions resolved within the init (via INTERACTIONS_EDGELIST)
#          or via a separate interactions file set here:
# interactions <- "inst/examples/interactions_init.txt"

res <- run_spatial_simulation(
  init_file = init,
  interactions_file = NULL,  # use inline settings if present
  output_prefix = "out/demo", 
  write_outputs = TRUE       # write CSVs, figures, report
)

# Outputs written under out/demo_<timestamp>*
# Results also returned in 'res' for programmatic use:
str(res$abund_matrix)
```

------------------------------------------------------------------------

## Typical workflow

1.  **Set parameters**: either with an **init file** (see vignette) or
    in-memory (`load_config()` → edit fields).
2.  **Run the simulator**: `run_spatial_simulation()`.
3.  **Use the outputs**:
    - `res$abund_matrix` (site × species),
    - `res$site_coords` (quadrat centroids),
    - `res$env_gradients` (grid with temp/elev/rain fields),
    - and map(s) via `plot_spatial_sampling()`.
4.  **(Optional) Advanced panel**: enable `P$ADVANCED_ANALYSIS <- TRUE`
    when you want the multi‑plot diagnostic image written to disk.

------------------------------------------------------------------------

## Function overview

| Function | Purpose |
|----|----|
| `run_spatial_simulation()` | End‑to‑end orchestrator; returns a results list and (optionally) writes CSVs/figures/report. |
| `generate_heterogeneous_distribution()` | Place individuals using abundances, gradients, dominant clustering, and interactions. |
| `create_abundance_matrix()` | Build site × species table by intersecting individuals with quadrats. |
| `calculate_quadrat_environment()` | Mean environmental conditions per quadrat from the grid. |
| `plot_spatial_sampling()` | Publication‑ready map of domain, individuals, quadrats, optional gradient fill. |
| `generate_full_report()` | Text report covering gradients, distributions, diversity, and Fisher log‑series validation. |
| `generate_advanced_panel()` | Patchwork panel: rank–abundance, occupancy–abundance, species–area, distance–decay, rarefaction. |

------------------------------------------------------------------------

## Interactions (neighbour effects)

You can specify interspecific effects **inline in the init file** or via
a simple **CSV edgelist**. See the dedicated vignette:

- <https://ajsmit.github.io/spesim/articles/spesim-interactions.html>

Key points:

- Rules are **directed** (`A -> B` may differ from `B -> A`).  
- Default is **neutral** (`1.0`); only list deviations.
- Set the global `INTERACTION_RADIUS` (0 disables the modifier).

------------------------------------------------------------------------

## Environmental gradients

Three synthetic gradients are available out-of-the-box (`temperature`,
`elevation`, `rainfall`), each normalised to \[0,1\] and reported in
common units for plots/reports. You can control per-species **optima**
and **tolerances** as scalars, vectors (by species), or named by
gradient. See the init-file vignette:

- <https://ajsmit.github.io/spesim/articles/spesim-init-parameters.html>

For a quick look at the generated grid:

``` r
domain <- create_sampling_domain()
env <- create_environmental_gradients(domain, resolution = 50, noise_level = 0.05)
head(env)
```

------------------------------------------------------------------------

## Citation

If you use **spesim**, please cite:

    @misc{spesim,
      author = {AJ Smit},
      title  = {spesim: Spatial Ecological Simulation in R},
      year   = {2025},
      url    = {https://github.com/ajsmit/spesim}
    }

------------------------------------------------------------------------

## License

MIT License — see `LICENSE`.
