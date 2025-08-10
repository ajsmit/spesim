#' Rank-Abundance Data
#'
#' @description
#' Construct a tidy table for rank-abundance plotting by combining
#' the **observed** species counts in `species_dist` with the
#' **theoretical** Fisher log-series abundances implied by `P`.
#'
#' @details
#' Observed abundances are computed with base R `table()`, ranked
#' in descending order, and annotated as `Source = "Observed"`.
#' Theoretical abundances are generated via
#' \code{\link{generate_fisher_log_series}} using parameters in `P`
#' (e.g., \code{N_SPECIES}, \code{N_INDIVIDUALS}, \code{FISHER_ALPHA},
#' \code{FISHER_X}, \code{DOMINANT_FRACTION}) and annotated as
#' `Source = "Theoretical"`. The two data frames are row-bound.
#'
#' @param species_dist An `sf` point object with a character column
#'   `species` giving per-individual species identities.
#' @param P A named list as returned by \code{\link{load_config}}
#'   containing at least the fields used by
#'   \code{\link{generate_fisher_log_series}}:
#'   \code{N_SPECIES}, \code{N_INDIVIDUALS}, \code{DOMINANT_FRACTION},
#'   \code{FISHER_ALPHA}, and \code{FISHER_X}.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{`Rank`}{Integer species rank (1 = most abundant).}
#'   \item{`Abundance`}{Non-negative integer abundance.}
#'   \item{`Source`}{`"Observed"` or `"Theoretical"`.}
#' }
#'
#' @seealso \code{\link{plot_rank_abundance}},
#'   \code{\link{generate_fisher_log_series}}
#'
#' @examples
#' \dontrun{
#' ra <- calculate_rank_abundance(species_dist, P)
#' head(ra)
#' }
#' @export
calculate_rank_abundance <- function(species_dist, P) {
  observed_counts <- table(species_dist$species)
  observed_data <- as.data.frame(observed_counts, stringsAsFactors = FALSE)
  names(observed_data) <- c("Species", "Abundance")
  observed_data <- observed_data |>
    arrange(desc(Abundance)) |>
    mutate(Rank = dplyr::row_number(), Source = "Observed") |>
    select(Rank, Abundance, Source)

  theoretical_abundances <- generate_fisher_log_series(P$N_SPECIES, P$N_INDIVIDUALS, P$DOMINANT_FRACTION, P$FISHER_ALPHA, P$FISHER_X)
  theoretical_data <- tibble::tibble(
    Abundance = sort(as.numeric(theoretical_abundances), decreasing = TRUE),
    Rank = seq_along(theoretical_abundances),
    Source = "Theoretical"
  )
  dplyr::bind_rows(observed_data, theoretical_data)
}

#' Occupancy-Abundance Table
#'
#' @description
#' Summarise total abundance and site occupancy (number of sites with
#' non-zero counts) for each species from a site x species matrix.
#'
#' @param abund_matrix A data frame where the first column is `site`
#'   (site / quadrat identifier) and the remaining columns are species
#'   abundances (non-negative integers).
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{`Species`}{Species (column) name.}
#'   \item{`TotalAbundance`}{Column sum across sites.}
#'   \item{`Occupancy`}{Number of sites with abundance \eqn{> 0}.}
#' }
#'
#' @seealso \code{\link{plot_occupancy_abundance}},
#'   \code{\link{create_abundance_matrix}}
#'
#' @examples
#' \dontrun{
#' oa <- calculate_occupancy_abundance(abund_matrix)
#' head(oa)
#' }
#' @export
calculate_occupancy_abundance <- function(abund_matrix) {
  abund_numeric <- abund_matrix[, -which(names(abund_matrix) == "site"), drop = FALSE]
  data.frame(
    Species = names(abund_numeric),
    TotalAbundance = colSums(abund_numeric),
    Occupancy = colSums(abund_numeric > 0),
    row.names = NULL
  )
}

#' Species-Area (Accumulation) Data
#'
#' @description
#' Compute a species accumulation curve (SAR) using
#' \code{vegan::specaccum(method = "random")}.
#'
#' @details
#' Uses 100 random permutations by default to estimate mean richness
#' and its standard deviation as a function of the number of sites.
#'
#' @param abund_matrix A site x species abundance data frame where the first
#'   column is `site` and remaining columns are species abundances.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{`Sites`}{Number of sites sampled.}
#'   \item{`Richness`}{Mean cumulative species richness.}
#'   \item{`SD`}{Standard deviation across permutations.}
#' }
#'
#' @seealso \code{\link{plot_species_area}}, \code{vegan::specaccum}
#'
#' @examples
#' \dontrun{
#' sar <- calculate_species_area(abund_matrix)
#' head(sar)
#' }
#' @export
calculate_species_area <- function(abund_matrix) {
  abund_numeric <- abund_matrix[, -which(names(abund_matrix) == "site"), drop = FALSE]
  sar_curve <- vegan::specaccum(abund_numeric, method = "random", permutations = 100)
  data.frame(Sites = sar_curve$sites, Richness = sar_curve$richness, SD = sar_curve$sd)
}

#' Distance-Decay Data
#'
#' @description
#' Pair geographic distances between sites with community dissimilarities
#' (Sorensen index computed as binary Bray-Curtis) for distance-decay plots.
#'
#' @param abund_matrix A site x species abundance data frame where the first
#'   column is `site` and the remaining columns are abundances.
#' @param site_coords A data frame with numeric columns `x` and `y`
#'   giving site coordinates (in a projected CRS) in the same row order
#'   as `abund_matrix`.
#'
#' @return A data frame with two numeric columns:
#' \describe{
#'   \item{`Distance`}{Euclidean distance between site pairs.}
#'   \item{`Dissimilarity`}{Sorensen dissimilarity (0-1).}
#' }
#'
#' @seealso \code{\link{plot_distance_decay}}, \code{vegan::vegdist}
#'
#' @examples
#' \dontrun{
#' dd <- calculate_distance_decay(abund_matrix, site_coords)
#' head(dd)
#' }
#' @export
calculate_distance_decay <- function(abund_matrix, site_coords) {
  coords <- site_coords[, c("x", "y")]
  abund_numeric <- abund_matrix[, -which(names(abund_matrix) == "site"), drop = FALSE]
  geo_dist <- stats::dist(coords, method = "euclidean")
  comm_dissim <- vegan::vegdist(abund_numeric, method = "bray", binary = TRUE)
  data.frame(Distance = as.vector(geo_dist), Dissimilarity = as.vector(comm_dissim))
}

#' Rarefaction Curves Data
#'
#' @description
#' Generate per-site rarefaction curves (expected richness vs. sample size)
#' using \code{vegan::rarecurve}.
#'
#' @details
#' The result is returned in long (tidy) format with one row per site x
#' sample size point. The `SampleSize` values come from the `Subsample`
#' attribute provided by \code{vegan::rarecurve}.
#'
#' @param abund_matrix A site x species abundance data frame where the first
#'   column is `site` and remaining columns are species counts.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{`SiteID`}{Factor identifying the site (from `abund_matrix$site`).}
#'   \item{`SampleSize`}{Number of individuals subsampled.}
#'   \item{`RarefiedRichness`}{Expected species richness at that sample size.}
#' }
#'
#' @seealso \code{\link{plot_rarefaction}}, \code{vegan::rarecurve}
#'
#' @examples
#' \dontrun{
#' rr <- calculate_rarefaction(abund_matrix)
#' head(rr)
#' }
#' @export
calculate_rarefaction <- function(abund_matrix) {
  abund_numeric <- abund_matrix[, -which(names(abund_matrix) == "site"), drop = FALSE]
  site_ids <- abund_matrix$site
  rarefaction_list <- vegan::rarecurve(abund_numeric, step = 1)
  output_list <- vector("list", length(rarefaction_list))
  for (i in seq_along(rarefaction_list)) {
    richness_values <- rarefaction_list[[i]]
    sample_sizes <- attr(richness_values, "Subsample")
    output_list[[i]] <- data.frame(
      SiteID = as.factor(site_ids[i]),
      SampleSize = sample_sizes,
      RarefiedRichness = richness_values
    )
  }
  do.call(rbind, output_list)
}

#' Plot Rank-Abundance Curve
#'
#' @description
#' Draws observed and theoretical rank-abundance (SAD) curves on a
#' log-scaled abundance axis. The input should be the combined table
#' returned by \code{\link{calculate_rank_abundance}()}, which contains
#' one row per rank with a \code{Source} column distinguishing
#' \emph{Observed} vs \emph{Theoretical}.
#'
#' @param rank_abundance_data A data frame with (at minimum) the columns
#'   \code{Rank} (integer rank starting at 1), \code{Abundance}
#'   (non-negative numeric counts), and \code{Source} (factor/character
#'   with levels such as \code{"Observed"} and \code{"Theoretical"}).
#'
#' @details
#' The y-axis is displayed on a base-10 logarithmic scale. Lines and points
#' are styled by \code{Source} (solid vs dashed; different shapes and colors).
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{calculate_rank_abundance}}
#'
#' @examples
#' \dontrun{
#' ra <- calculate_rank_abundance(species_dist, P)
#' plot_rank_abundance(ra)
#' }
#' @export
plot_rank_abundance <- function(rank_abundance_data) {
  ggplot(rank_abundance_data, aes(x = Rank, y = Abundance, color = Source)) +
    geom_line(aes(linetype = Source), linewidth = 1.1) +
    geom_point(aes(shape = Source), size = 3, fill = "white", stroke = 1.2) +
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_color_manual(values = c("Observed" = "black", "Theoretical" = "#e41a1c")) +
    scale_linetype_manual(values = c("Observed" = "solid", "Theoretical" = "dashed")) +
    scale_shape_manual(values = c("Observed" = 21, "Theoretical" = 22)) +
    labs(
      title = "Species-Abundance Distribution (SAD)",
      subtitle = "Observed vs. Fisher's log-series",
      x = "Species Rank",
      y = "Abundance (log10)",
      color = "Distribution",
      linetype = "Distribution",
      shape = "Distribution"
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
}

#' Plot Occupancy-Abundance Relationship
#'
#' @description
#' Scatterplot of total abundance (per species) versus site occupancy
#' (number of quadrats in which the species occurs) with both axes on a
#' log scale. Adds an optional least-squares trend line for visual guidance.
#'
#' @param oa_data A data frame as returned by
#'   \code{\link{calculate_occupancy_abundance}} with columns
#'   \code{Species}, \code{TotalAbundance} (non-negative numeric),
#'   and \code{Occupancy} (non-negative integer).
#'
#' @details
#' If the input contains no observations (or all totals are zero), a
#' minimal placeholder plot is returned indicating that there is nothing
#' to draw. Otherwise both axes are shown on base-10 logarithmic scales.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{calculate_occupancy_abundance}}
#'
#' @examples
#' \dontrun{
#' oa <- calculate_occupancy_abundance(abund_matrix)
#' plot_occupancy_abundance(oa)
#' }
#' @export
plot_occupancy_abundance <- function(oa_data) {
  if (nrow(oa_data) == 0 || all(oa_data$TotalAbundance == 0)) {
    return(ggplot() +
      labs(title = "Occupancy-Abundance Relationship", subtitle = "No data to plot.") +
      theme_void())
  }
  ggplot(oa_data, aes(x = TotalAbundance, y = Occupancy)) +
    geom_point(alpha = 0.6, size = 3) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.8) +
    scale_x_log10(labels = scales::label_log()) +
    scale_y_log10(labels = scales::label_log()) +
    labs(
      title = "Occupancy-Abundance Relationship",
      x = "Total Abundance (log10)",
      y = "Sites Occupied (log10)"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
}

#' Plot Species-Area Relationship (SAR)
#'
#' @description
#' Plots the mean species accumulation (species-area) curve with a ribbon
#' denoting +/-1 standard deviation across permutations, based on the summary
#' returned by \code{\link{calculate_species_area}}.
#'
#' @param sar_data A data frame with columns \code{Sites} (number of
#'   quadrats sampled), \code{Richness} (mean cumulative richness),
#'   and \code{SD} (standard deviation across permutations).
#'
#' @details
#' The curve represents expected richness as sampling effort increases;
#' uncertainty is displayed as a semi-transparent ribbon.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{calculate_species_area}}
#'
#' @examples
#' \dontrun{
#' sar <- calculate_species_area(abund_matrix)
#' plot_species_area(sar)
#' }
#' @export
plot_species_area <- function(sar_data) {
  ggplot(sar_data, aes(x = Sites, y = Richness)) +
    geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.4) +
    geom_line(linewidth = 1.2) +
    labs(
      title = "Species-Area Relationship (SAR)",
      x = "Number of Quadrats",
      y = "Cumulative Species Richness"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
}

#' Plot Distance-Decay Relationship
#'
#' @description
#' Displays community dissimilarity versus geographic distance with a
#' smooth loess trend. The y-axis is constrained to \[0, 1\] to reflect
#' the range of Sorensen (binary Bray-Curtis) dissimilarity.
#'
#' @param decay_data A data frame as returned by
#'   \code{\link{calculate_distance_decay}} with columns
#'   \code{Distance} (pairwise Euclidean distances among quadrat centroids)
#'   and \code{Dissimilarity} (pairwise Sorensen dissimilarities).
#'
#' @details
#' Points show pairwise site comparisons; the loess smoother summarizes
#' the overall distance-decay pattern. Input should already exclude
#' self-comparisons.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{calculate_distance_decay}}
#'
#' @examples
#' \dontrun{
#' dd <- calculate_distance_decay(abund_matrix, site_coords)
#' plot_distance_decay(dd)
#' }
#' @export
plot_distance_decay <- function(decay_data) {
  ggplot(decay_data, aes(x = Distance, y = Dissimilarity)) +
    geom_point(alpha = 0.3, shape = 16) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1.1) +
    ylim(0, 1) +
    labs(
      title = "Distance-Decay of Community Similarity",
      x = "Geographic Distance",
      y = "Community Dissimilarity (Sorensen)"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
}

#' Plot Rarefaction Curves
#'
#' @description
#' Draws per-site rarefaction curves (expected richness as a function of
#' the number of individuals sampled) using the long-format output from
#' \code{\link{calculate_rarefaction}}.
#'
#' @param rarefaction_data A data frame with columns \code{SiteID}
#'   (factor/character site label), \code{SampleSize} (non-negative
#'   integer), and \code{RarefiedRichness} (expected species count).
#'
#' @details
#' Each site/quadrat is drawn as a separate line. Colors are mapped to
#' \code{SiteID} using a discrete viridis palette for readability.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{calculate_rarefaction}}
#'
#' @examples
#' \dontrun{
#' rr <- calculate_rarefaction(abund_matrix)
#' plot_rarefaction(rr)
#' }
#' @export
plot_rarefaction <- function(rarefaction_data) {
  ggplot(rarefaction_data, aes(x = SampleSize, y = RarefiedRichness, group = SiteID, color = SiteID)) +
    geom_line(linewidth = 0.8) +
    scale_color_viridis_d(option = "plasma") +
    labs(
      title = "Rarefaction Curves by Quadrat",
      x = "Individuals Sampled",
      y = "Expected Species (Rarefied)",
      color = "Quadrat ID"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
}
