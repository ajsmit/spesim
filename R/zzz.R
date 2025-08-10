# R/zzz.R

#' @keywords internal
"_PACKAGE"

# ---- Graphics & plotting ----
# Bring in all ggplot2 so helpers like ggplot(), aes(), labs(), ylim(),
# theme_bw(), geom_*(), scale_*(), etc. are available package-wide.
# This is what fixes the repeated "could not find function 'ylim'/'labs'..." errors.
#' @import ggplot2

#' @importFrom patchwork wrap_plots plot_layout plot_annotation

#' @importFrom viridis scale_fill_viridis

# ---- Data wrangling ----
# dplyr verbs used unqualified
#' @importFrom dplyr arrange mutate select group_by summarise n_distinct filter
#' @importFrom dplyr left_join across count row_number
#' @importFrom dplyr desc


#' @importFrom FNN get.knnx
#' @importFrom RANN nn2

# tidyr / tibble
#' @importFrom tidyr pivot_wider
#' @importFrom tibble tibble

# ---- Stats / utils ----
#' @importFrom stats cor cor.test dist
#' @importFrom stats rnorm runif var
#' @importFrom vegan diversity vegdist fisher.alpha
#' @importFrom utils write.csv read.csv globalVariables
#' @importFrom grDevices col2rgb

# utils
#' @importFrom utils globalVariables
#' @importFrom utils capture.output

# ---- sf helpers used unqualified in a few places ----
# (Most sf calls are qualified, but st_drop_geometry() is sometimes unqualified.)
#' @importFrom sf st_drop_geometry

NULL

# Silence R CMD check for NSE column names used in dplyr/ggplot2
utils::globalVariables(c(
  # data columns used across verbs/plots
  "site", "quadrat_id", "species", "x", "y", "X", "Y",
  "Rank", "Abundance", "Source",
  "temperature_C", "elevation_m", "rainfall_mm",
  "richness", "n_ind", "Var1", "Freq",

  # NSE/tidy‑eval columns created or referenced without quotes
  "abundance",
  "TotalAbundance", "Occupancy",
  "Sites", "Richness", "SD",
  "SiteID", "value",
  "Count", "Distance", "Dissimilarity",
  "SampleSize", "RarefiedRichness",

  # columns produced by spatstat data frames used in ggplot
  "r", "border", "trans",

  # helpers used in label/transform contexts (e.g., scales::math_format)
  ".x"
))

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("spesim loaded - try run_spatial_simulation() to generate a simulation.")
}

.onLoad <- function(libname, pkgname) {
  # no special load behavior needed
}
