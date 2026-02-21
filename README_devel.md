# spesim — development log & plan (README_devel)

This file is a lightweight developer-facing log of what has been done **since Claw started looking into `spesim` (2026-02-09)**, plus a short forward plan.

---

## Current status (as of 2026-02-21)

- Package version: **0.5.2**
- All CRAN checks pass: `0 errors | 0 warnings | 0 notes` (R 4.5.2, macOS Tahoe 26.3, arm64).
- Full vignette rebuild succeeds under `R CMD check`.
- `devtools::document()` runs successfully.
- Changes committed and pushed to GitHub (`main`).

## Profiling history (benchmark: N=1,500, S=15, hybrid model, `profvis` 10 ms interval)

| Version | Change | Wall time | Speedup vs v0.5.1 |
|---------|--------|-----------|-------------------|
| v0.5.1 | C++ PIP engine + vectorised domain checks | 17,690 ms | 1× (baseline) |
| v0.5.2 r1 | Vectorised env-acceptance step | 6,010 ms | 2.9× |
| v0.5.2 r2 | Vectorised candidate displacement sampling | **2,600 ms** | **6.8×** |

## Current profile (v0.5.2, self time)

| Rank | Self | % | Function | Notes |
|------|------|---|----------|-------|
| 1 | 480 ms | 18.5% | `pip_cpp` | C++ PIP engine — already optimal |
| 2 | 310 ms | 11.9% | `spesim_simulate_neutral_recruitment` | R interpreter overhead of the loop body |
| 3 | 170 ms | 6.5% | `.env_accept_prob_vec` | Residual: `findInterval` + matrix lookup |
| 4 | 150 ms | 5.8% | `<GC>` | Heap pressure from small per-step allocations |
| 5 | 120 ms | 4.6% | `cbind` | Candidate matrix construction |

## Next frontier

The hot path is now dominated by `pip_cpp` (the C++ PIP engine, already compiled) and
the R interpreter overhead of the simulation loop itself. Meaningful further gains would
require porting the inner loop to Rcpp/C++, or reducing `batch_size` adaptively when the
inside-domain acceptance rate is high. The `cbind()` call inside `.propose_displacement_batch()`
is a minor candidate (pre-allocating a matrix and filling it would avoid the copy), but
the gain would be small relative to the remaining R-loop overhead.

---

## Refactoring and Hardening (v0.5.0)

A major refactoring of the simulation engine was undertaken to improve correctness, modularity, and performance.

### 1) Code Restructuring
- **`generate_heterogeneous_distribution()` was refactored**: Internal helper functions for point process simulation and other utilities were extracted from `R/abundance.R` and moved into a new file, `R/simulation_helpers.R`.
- **Neutral Model Extracted**: The complex neutral model simulator (`.simulate_neutral_recruitment`) was moved from `R/abundance.R` into its own dedicated file, `R/simulate_neutral_recruitment.R`, and renamed to `spesim_simulate_neutral_recruitment`.

### 2) Simulation Model Improvements
- **Canonical Point Processes**: The simple "surrogate" simulators for the Strauss and Geyer point processes were replaced. The simulation engine now calls `simulate_points_dispatch`, which uses the fast, C++-backed MCMC implementations for these models where available, falling back to the R versions otherwise. This ensures the simulations are more theoretically sound.
- **Neutral Model Convergence**: The neutral model simulation (`spesim_simulate_neutral_recruitment`) was changed from a fixed-step heuristic to a convergence-based process. It now runs until the community composition stabilizes (measured by Bray-Curtis dissimilarity) or a maximum number of steps is reached. New parameters (`NEUTRAL_MAX_STEPS`, `NEUTRAL_CONVERGENCE_INTERVAL`, etc.) were added to `load_config` to control this.
- **Parameter Disambiguation**: Corrected a bug where the `OTHERS_S` parameter was used for both Strauss and Geyer models, which have different constraints. A new parameter, `OTHERS_STRAUSS_GAMMA`, was introduced for the Strauss process to resolve this ambiguity.

### 3) Vignette Performance and Build Process
- **Build Failures**: The `R CMD build` process was consistently timing out due to slow vignette execution.
- **Performance Analysis**: A script was run to time each vignette's build process individually. This identified `spesim-constrained-landscapes.Rmd` and `spesim-recipes-point-processes.Rmd` as the two slowest vignettes, with combined run times exceeding the build timeout.
- **Vignette Optimization**: The slow vignettes were optimized by reducing the number of individuals and quadrats in the examples, significantly decreasing their run time while preserving their instructional value.
- **Successful Build and Check**: After optimization, the package now successfully passes `R CMD build` and `R CMD check` with all vignettes included.

### Recommended Build & Test Workflow

For a full check that mimics CRAN's process, the following is recommended. Note that building vignettes can be slow.

```bash
# Regenerate documentation from Roxygen comments
Rscript -e 'devtools::document()'

# Build the source package, including vignettes
R CMD build .

# Check the resulting tarball.
# The filename will change with version numbers.
# The _R_CHECK_FORCE_SUGGESTS_=false is to avoid dependency errors on suggested packages
# that are not strictly needed.
_R_CHECK_FORCE_SUGGESTS_=false R CMD check spesim_*.tar.gz --no-manual
```
