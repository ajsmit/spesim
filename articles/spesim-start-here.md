# Start here

This is the quickest path through the package, depending on what you’re
trying to do.

## 1) First: understand what spesim is (and isn’t)

- **Model card:**
  [`vignette("spesim-model-card")`](https://ajsmit.github.io/spesim/articles/spesim-model-card.md)

If you’re teaching with spesim, start there: it sets expectations and
prevents over-interpretation.

## 2) Run your first simulation

- **Workflow overview:**
  [`vignette("spesim-workflow")`](https://ajsmit.github.io/spesim/articles/spesim-workflow.md)

Minimal, in-memory run:

``` r
library(spesim)

P <- load_config(system.file("examples/spesim_init_basic.txt", package = "spesim"))
res <- spesim_run(P, write_outputs = FALSE, seed = 77)

plot_spatial_sampling(res$domain, res$species_dist, res$quadrats, res$P)
```

## 3) Choose a sampling design

- **Quadrat placement schemes:**
  [`vignette("spesim-quadrat-placement")`](https://ajsmit.github.io/spesim/articles/spesim-quadrat-placement.md)

## 4) Add structure (gradients, interactions, point processes)

- **Environmental gradients:**
  [`vignette("spesim-env-gradients")`](https://ajsmit.github.io/spesim/articles/spesim-env-gradients.md)
- **Interactions:**
  [`vignette("spesim-interactions")`](https://ajsmit.github.io/spesim/articles/spesim-interactions.md)
- **Point processes:**
  [`vignette("spesim-point-processes")`](https://ajsmit.github.io/spesim/articles/spesim-point-processes.md)

## 5) Validate that the patterns match your intent

- **Validation & sanity checks:**
  [`vignette("spesim-validation")`](https://ajsmit.github.io/spesim/articles/spesim-validation.md)

## 6) Stability / teaching materials

If you rely on specific functions in teaching content:

- **Public API (what’s stable):**
  [`vignette("spesim-public-api")`](https://ajsmit.github.io/spesim/articles/spesim-public-api.md)
