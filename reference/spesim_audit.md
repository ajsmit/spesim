# Audit whether a spesim run matches the intended qualitative regime

spesim is designed for teaching and method testing: you often want to
know whether the *realised* simulation actually exhibits the qualitative
regime you think you requested (e.g., clustering vs inhibition; strong
vs weak environmental filtering; heavy boundary exclusion in a sampling
scheme).

`spesim_audit()` computes lightweight, dependency-minimal diagnostics
that can be used as a conceptual "tutor": it reports what the simulator
instantiated, and quantifies key patterns.

- **Spatial structure:** nearest-neighbour (NN) diagnostics per species,
  comparing the observed mean NN distance to the CSR expectation.

- **Environmental filtering:** checks whether quadrat abundances for
  gradient-responsive species peak near their configured optima.

- **Sampling design:** reports boundary exclusion / rejection rates for
  the chosen quadrat placement scheme (when available).

For more formal point-process diagnostics (Ripley's K, pair correlation
functions, edge corrections in irregular domains), use `spatstat`
tooling on the exported point patterns.

## Usage

``` r
spesim_audit(res, species = "all", nn_k = 1, diagnostics = c("nn", "spatstat"))
```

## Arguments

- res:

  A `spesim_result` list returned by
  [`spesim_run()`](https://ajsmit.github.io/spesim/reference/spesim_run.md)
  or
  [`run_spatial_simulation()`](https://ajsmit.github.io/spesim/reference/run_spatial_simulation.md).

- species:

  Character vector of species labels to audit for spatial structure.
  Default `"all"` audits all species present.

- nn_k:

  Integer; the neighbour order used for NN distances. Default 1.

- diagnostics:

  Character vector controlling which diagnostics to compute. Default
  `c("nn", "spatstat")` will try to compute both nearest-neighbour and
  spatstat-based diagnostics; if spatstat is not installed, those
  diagnostics are skipped with a message. Use `"nn"` only for a
  lightweight audit.

## Value

A list of class `spesim_audit` with components:

- `spatial`: data.frame of NN diagnostics per species

- `spatial_spatstat`: data.frame of edge-corrected diagnostics (if
  requested)

- `filtering`: data.frame of environment-abundance checks (may be empty)

- `sampling`: a list describing quadrat placement rejection/exclusion
  metrics

- `regime`: green/amber/red classification (see
  [`spesim_regime()`](https://ajsmit.github.io/spesim/reference/spesim_regime.md))
