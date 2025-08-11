spesim: Spatial Ecological Simulation in R
================
Your Name
2025-08-11

# spesim <img src="man/figures/logo.png" align="right" width="120" />

**spesim** is an R package for generating, analysing, and visualising
*spatially heterogeneous ecological communities*.  
It is designed for simulation-based teaching, methodological testing,
and exploratory research.

## Features

- Simulate species distributions across an arbitrary polygonal domain.
- Include **environmental gradients**, **dominance structure**, and
  **spatial autocorrelation**.
- Flexible spatial point generation: CSR, cluster processes, inhibition
  models.
- Automatically produce **species × site abundance matrices**,
  **environmental summaries**, and **diversity metrics**.
- High-quality plotting functions for:
  - Rank–abundance curves
  - Occupancy–abundance relationships
  - Species–area curves
  - Distance–decay patterns
  - Rarefaction curves
- **Advanced analysis panel**: multi-plot patchwork summary with
  optional point–process diagnostics.
- Full-text simulation reports for reproducible summaries.

## Installation

You can install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("yourusername/spesim")
```

The package imports several analysis and plotting libraries; see the
DESCRIPTION file for full details. Optional point-process diagnostics
require:

``` r
install.packages(c("spatstat.geom", "spatstat.core"))
```

## Quick start

``` r
library(spesim)

# Define simulation parameters
P <- list(
  N_SPECIES = 10,
  N_INDIVIDUALS = 1000,
  DOMINANT_FRACTION = 0.3,
  FISHER_ALPHA = 5,
  FISHER_X = 0.6
)

# Generate a simple rectangular domain
domain <- sf::st_as_sf(
  data.frame(id = 1),
  coords = cbind(c(0, 100, 100, 0, 0), c(0, 0, 100, 100, 0)),
  crs = 3857
) |> sf::st_cast("POLYGON")

# Run the spatial simulation
res <- run_spatial_simulation(P, domain)

# Visualise the output
generate_advanced_panel(res)
```

## Function overview

| Function | Purpose |
|----|----|
| `run_spatial_simulation()` | Main simulation engine; returns a list of results containing parameters, species distributions, quadrats, abundance matrices, and coordinates. |
| `generate_heterogeneous_distribution()` | Generate species point locations given a chosen spatial process model (CSR, cluster process, inhibition model, etc.). |
| `create_abundance_matrix()` | Build site × species abundance matrix from species point data and quadrat polygons. |
| `calculate_quadrat_environment()` | Summarise environmental conditions per quadrat from environmental grid data. |
| `generate_full_report()` | Produce a plain-text, multi-section analysis report covering gradients, species distributions, diversity metrics, and model validation. |
| `generate_advanced_panel()` | Assemble a multi-plot patchwork diagnostic panel, optionally with point–process diagnostics if `spatstat` is available. |
| `plot_rank_abundance()` | Plot observed vs. theoretical rank–abundance curves on a log-scaled axis. |
| `plot_occupancy_abundance()` | Plot occupancy–abundance relationships (log–log) with optional trendline. |
| `plot_species_area()` | Plot species–area relationships with mean curve ± 1 SD. |
| `plot_distance_decay()` | Plot distance–decay relationships using Sørensen dissimilarity vs. geographic distance. |
| `plot_rarefaction()` | Plot per-site rarefaction curves showing expected richness vs. sample size. |

## Methodological note: from ad-hoc clustering to point processes

In early versions, spesim generated clusters using simple “hand-rolled”
logic for dominant species. You can now opt for principled point process
models:

- Neyman–Scott / Thomas process: parent–offspring clustering, controlled
  by parent intensity & offspring spread.
- Strauss process: inhibition between points, tunable interaction radius
  & strength.
- Geyer saturation process: combination of clustering and inhibition
  effects.

These models give you a tunable pair correlation function (g(r)) that
you can validate against empirical or simulated data using:

- Ripley’s K function (K(r))
- Transformed (L(r)-r)
- (g(r)) from the pair correlation function

This is integrated in `generate_advanced_panel()` if the **spatstat**
packages are installed.

## Generating environmental gradients

You can attach simulated environmental rasters or grids, then extract
quadrat summaries:

``` r
domain <- create_sampling_domain()
env <- create_environmental_gradients(domain, resolution = 50, noise_level = 0.05)
head(env)
```

## Diagnostics panel

If **spatstat** is available, the advanced panel adds an extra row:

- Ripley’s K (border correction): deviation from CSR.
- L(r) – r: positive values suggest clustering; negative values suggest
  inhibition.
- Pair correlation g(r): \>1 clustering, \<1 inhibition.

## Citation

If you use spesim in your work, please cite it:

@misc{spesim, author = {AJ Smit}, title = {spesim: Spatial Ecological
Simulation in R}, year = {2025}, url =
{<https://github.com/ajsmit/spesim>} }

## License

MIT License — see LICENSE file for details.
