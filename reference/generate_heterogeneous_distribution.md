# Generate a heterogeneous community (point process + environment + interactions)

Generates individual locations and species identities over an arbitrary
polygon domain by combining:

1.  a spatial point process for baseline locations (clustered,
    inhibited, or Poisson),

2.  environmental filtering for gradient-responsive species, and

3.  local interspecific interactions within a neighbourhood radius.

The dominant species "A" and the pool of non-dominant species can use
different point-process models. Environmental suitability is applied as
a Gaussian preference around a per-species optimum, and local
interactions are incorporated as a multiplicative modifier computed from
nearby already-assigned individuals.

## Usage

``` r
generate_heterogeneous_distribution(domain, P)
```

## Arguments

- domain:

  An sf polygon (or multipolygon) defining the sampling domain; must
  have a valid CRS (projected coordinates recommended).

- P:

  A fully materialized parameter list, typically from
  [`load_config()`](https://ajsmit.github.io/spesim/reference/load_config.md),
  containing at least:

  - *Community size:* `N_SPECIES`, `N_INDIVIDUALS`, `DOMINANT_FRACTION`,
    `FISHER_ALPHA`, `FISHER_X`.

  - *Environment:* `SAMPLING_RESOLUTION`, `ENVIRONMENTAL_NOISE`, and
    `GRADIENT` (tibble with `species`, `gradient` in {temperature,
    elevation, rainfall}, `optimum` in \\\[0,1\]\\, and `tol > 0`). If
    absent, species are treated as neutral.

  - *Point-process selection (strings):* `SPATIAL_PROCESS_A` and
    `SPATIAL_PROCESS_OTHERS` in `"poisson"`, `"thomas"`, `"strauss"`, or
    `"geyer"`.

  - *Thomas (A) params (if used):* `A_PARENT_INTENSITY` (parents per
    area; optional), `A_MEAN_OFFSPRING` (mean children per parent),
    `A_CLUSTER_SCALE` (Gaussian sd of offspring displacement; map
    units).

  - *Strauss/Geyer (others) params (if used):* `OTHERS_R` (interaction
    radius), `OTHERS_GAMMA` (interaction parameter; `< 1` inhibition,
    `> 1` attraction), `OTHERS_S` (Geyer saturation count; ignored by
    Strauss), `OTHERS_BETA` (baseline intensity/multiplier; optional).

  - *Local interactions:* `INTERACTION_RADIUS` (map units) and
    `INTERACTION_MATRIX` (S x S numeric, dimnames = species letters).

## Value

An sf POINT layer with a character column `species` and appended
environmental columns (e.g. `temperature_C`, `elevation_m`,
`rainfall_mm`). Rows correspond to simulated individuals retained after
assignment.

## Details

**Abundances.** Total individuals per species are drawn with
[`generate_fisher_log_series()`](https://ajsmit.github.io/spesim/reference/generate_fisher_log_series.md),
allocating a fixed fraction to species "A" and distributing the
remainder by a log-series across `B, C, ...`.

**Baseline locations (point processes).** Locations are simulated using
the process names in `P`. Supported values (case-insensitive):

- `"poisson"`:

  Homogeneous Poisson inside `domain`.

- `"thomas"`:

  Thomas (Neyman-Scott) cluster process; uses parent intensity (derived
  if missing), mean offspring, and Gaussian cluster scale. Implemented
  internally.

- `"strauss"`:

  Inhibition via a sequential-inhibition surrogate using interaction
  radius `OTHERS_R` and inhibition strength `OTHERS_S` in \\(0,1\]\\.
  Lower `OTHERS_S` increases inhibition.

- `"geyer"`:

  Geyer saturation surrogate: acceptance probability proportional to
  \\\gamma^{\min(m, s)}\\ where \\m\\ neighbours fall within `OTHERS_R`
  and saturation count `OTHERS_S`. Implemented internally.

**Environmental filtering.** For species listed in `P$GRADIENT`, the
assignment probability for each species is proportional to \\\exp\\-(x -
\mu)^2 / (2 \sigma^2)\\\\, where \\x\\ is the normalized environmental
value at the point, \\\mu\\ the species optimum, and \\\sigma\\ the
tolerance. Environmental values are attached by nearest neighbour join
to the grid returned by
[`create_environmental_gradients()`](https://ajsmit.github.io/spesim/reference/create_environmental_gradients.md).

**Local interactions.** For each candidate point and focal species, the
abundance-independent interaction modifier is the geometric mean of the
corresponding coefficients in `P$INTERACTION_MATRIX` for neighbours
found within `P$INTERACTION_RADIUS` (using up to 5 nearest
already-assigned individuals). Coefficients `> 1` favour co-occurrence;
`< 1` penalize it.

**Tie-breaking and robustness.** If all weights for a step are
non-finite or non-positive, a uniform assignment is used for that step
to avoid dead ends. Only points with a non-empty `species` are returned.

## Notes

- Use a projected CRS (e.g. metres) so that process radii and cluster
  scales are in meaningful linear units.

- When `SPATIAL_PROCESS_OTHERS` is inhibitory and `n` is large relative
  to `OTHERS_R` and domain area, you may hit feasibility limits; adjust
  `n`, `r`, or choose a different process.

- Setting `INTERACTION_RADIUS = 0` or an all-ones matrix disables local
  interaction effects.

## See also

[`create_environmental_gradients()`](https://ajsmit.github.io/spesim/reference/create_environmental_gradients.md),
[`generate_fisher_log_series()`](https://ajsmit.github.io/spesim/reference/generate_fisher_log_series.md),
[`load_config()`](https://ajsmit.github.io/spesim/reference/load_config.md),
[`load_interactions()`](https://ajsmit.github.io/spesim/reference/load_interactions.md)

## Examples

``` r
if (FALSE) { # \dontrun{
P <- load_config("simul_init.txt")
domain <- create_sampling_domain()

# Example: A clustered (Thomas), others mildly inhibited (Strauss)
P$SPATIAL_PROCESS_A <- "thomas"
P$A_PARENT_INTENSITY <- NA
P$A_MEAN_OFFSPRING <- 10
P$A_CLUSTER_SCALE <- 0.8

P$SPATIAL_PROCESS_OTHERS <- "strauss"
P$OTHERS_R <- 0.6
P$OTHERS_S <- 0.5

pts <- generate_heterogeneous_distribution(domain, P)
plot(sf::st_geometry(pts), pch = 16, cex = 0.4)
} # }
```
