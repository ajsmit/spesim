#' Rank-Abundance Data
#'
#' @description
#' Construct a tidy table for rank-abundance plotting by combining
#' the **observed** species counts in `species_dist` with the
#' **theoretical** SAD abundances implied by `P`.
#'
#' @details
#' Observed abundances are computed with base R `table()`, ranked
#' in descending order, and annotated as `Source = "Observed"`.
#' Theoretical abundances are generated via
#' \code{\link{generate_sad}} using parameters in `P`
#' (e.g., \code{N_SPECIES}, \code{N_INDIVIDUALS}, and \code{SAD_MODEL}; plus
#' model-specific parameters such as \code{FISHER_ALPHA}, \code{FISHER_X},
#' \code{DOMINANT_FRACTION}) and annotated as
#' `Source = "Theoretical"`. The two data frames are row-bound.
#'
#' @param species_dist An `sf` point object with a character column
#'   `species` giving per-individual species identities.
#' @param P A named list as returned by \code{\link{load_config}}
#'   containing at least the fields used by
#'   \code{\link{generate_sad}}:
#'   \code{N_SPECIES}, \code{N_INDIVIDUALS}, and \code{SAD_MODEL} (plus any
#'   model-specific parameters).
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{`Rank`}{Integer species rank (1 = most abundant).}
#'   \item{`Abundance`}{Non-negative integer abundance.}
#'   \item{`Source`}{`"Observed"` or `"Theoretical"`.}
#' }
#'
#' @seealso \code{\link{plot_rank_abundance}},
#'   \code{\link{generate_sad}},
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

  theoretical_abundances <- generate_sad(
    n_species = P$N_SPECIES,
    n_individuals = P$N_INDIVIDUALS,
    model = P$SAD_MODEL %||% "fisher",
    dominant_fraction = P$DOMINANT_FRACTION,
    alpha = P$FISHER_ALPHA,
    x = P$FISHER_X,
    k = P$GEOMETRIC_K,
    exponent = P$ZIPF_EXPONENT,
    q = P$ZIPF_Q,
    meanlog = P$LOGNORMAL_MEANLOG,
    sdlog = P$LOGNORMAL_SDLOG,
    shape = P$POIGAMMA_SHAPE,
    rate = P$POIGAMMA_RATE,
    theta = P$ZSM_THETA,
    m = P$ZSM_M,
    sad = P$SAD_VECTOR
  )
  theoretical_abundances <- theoretical_abundances[theoretical_abundances > 0]
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

  # vegan::specaccum() expects >= 2 sites for permutation-based accumulation.
  if (nrow(abund_numeric) < 2) {
    richness_1 <- sum(colSums(abund_numeric > 0, na.rm = TRUE) > 0)
    return(data.frame(Sites = 1, Richness = richness_1, SD = NA_real_))
  }

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
#' @param metric Character scalar. One of:
#'   \itemize{
#'     \item \code{"auto"} (default): use along-path distance if
#'       a graph distance matrix or \code{site_coords$linear_pos} exists, else Euclidean.
#'     \item \code{"euclidean"}: force Euclidean distance on \code{x,y}.
#'     \item \code{"along_path"}: use absolute differences in
#'       \code{site_coords$linear_pos}; if \code{linear_wrap} attribute is TRUE,
#'       use wrapped circular distance. If attribute
#'       \code{attr(site_coords, "graph_dist_matrix")} is present, that matrix is
#'       used instead.
#'   }
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
calculate_distance_decay <- function(abund_matrix, site_coords, metric = c("auto", "euclidean", "along_path")) {
  metric <- match.arg(metric)
  # Ensure coordinate rows align with abundance rows by site id
  if ("site" %in% names(site_coords) && "site" %in% names(abund_matrix)) {
    ord <- match(abund_matrix$site, site_coords$site)
    if (any(is.na(ord))) {
      stop("site_coords$site must contain all sites in abund_matrix$site")
    }
    site_coords <- site_coords[ord, , drop = FALSE]
  }

  coords <- site_coords[, c("x", "y"), drop = FALSE]
  abund_numeric <- abund_matrix[, -which(names(abund_matrix) == "site"), drop = FALSE]

  # Note: we keep empty sites (all-zero rows). For distance–decay, we compute a
  # Sorensen dissimilarity matrix that defines empty–empty dissimilarity as 0 and
  # empty–nonempty as 1, avoiding NA outputs.

  if (nrow(coords) < 2 || nrow(abund_numeric) < 2) {
    return(data.frame(Distance = numeric(0), Dissimilarity = numeric(0)))
  }

  use_metric <- metric
  if (identical(use_metric, "auto")) {
    has_graph <- is.matrix(attr(site_coords, "graph_dist_matrix"))
    use_metric <- if (has_graph || "linear_pos" %in% names(site_coords)) "along_path" else "euclidean"
  }

  if (identical(use_metric, "along_path")) {
    gdist <- attr(site_coords, "graph_dist_matrix")
    if (is.matrix(gdist)) {
      if (nrow(gdist) != nrow(coords) || ncol(gdist) != nrow(coords)) {
        stop("graph_dist_matrix must be square and match number of sites.")
      }
      if (isTRUE(attr(site_coords, "graph_directed"))) {
        gdist <- pmin(gdist, t(gdist))
      }
      geo_dist <- stats::as.dist(gdist)
    } else {
      if (!"linear_pos" %in% names(site_coords)) {
        stop("metric='along_path' requires site_coords$linear_pos or graph_dist_matrix.")
      }
      lp <- as.numeric(site_coords$linear_pos)
      dmat <- abs(outer(lp, lp, "-"))
      if (isTRUE(attr(site_coords, "linear_wrap"))) {
        dmat <- pmin(dmat, 1 - dmat)
      }
      geo_dist <- stats::as.dist(dmat)
    }
  } else {
    geo_dist <- stats::dist(coords, method = "euclidean")
  }
  # Sorensen dissimilarity on presence/absence (binary), with well-defined
  # behaviour for empty rows.
  pa <- (as.matrix(abund_numeric) > 0) * 1
  tot <- rowSums(pa)
  shared <- pa %*% t(pa)
  denom <- outer(tot, tot, "+")
  sim <- (2 * shared) / denom
  sim[denom == 0] <- 1 # empty-empty -> similarity 1
  dissim <- 1 - sim
  diag(dissim) <- 0
  comm_dissim <- stats::as.dist(dissim)

  out <- data.frame(Distance = as.vector(geo_dist), Dissimilarity = as.vector(comm_dissim))
  out[is.finite(out$Distance) & is.finite(out$Dissimilarity), , drop = FALSE]
}

#' Rarefaction Curves Data
#'
#' @description
#' Generate per-site rarefaction curves (expected richness vs. sample size)
#' using \code{vegan::rarefy} (i.e., the same computation as
#' \code{vegan::rarecurve}, but without any plotting side-effects).
#'
#' @details
#' The result is returned in long (tidy) format with one row per site x
#' sample size point. The `SampleSize` values come from the `Subsample`
#' attribute attached by \code{vegan::rarefy}.
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
#' @seealso \code{\link{plot_rarefaction}}, \code{vegan::rarefy}
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

  # NOTE: We intentionally do NOT call vegan::rarecurve() here.
  # rarecurve() opens/uses a graphics device even when you only want the data,
  # which causes unwanted standalone plots when used inside other plotting
  # workflows (e.g., generate_advanced_panel()).
  x <- as.matrix(abund_numeric, rownames.force = TRUE)

  if (!isTRUE(all.equal(x, round(x)))) {
    stop("calculate_rarefaction() requires integer count data (site x species).")
  }
  x <- round(x)

  tot <- rowSums(x)

  output_list <- vector("list", nrow(x))
  for (i in seq_len(nrow(x))) {
    if (tot[i] <= 0) {
      # Return a 0-row slice with the correct column types.
      output_list[[i]] <- data.frame(
        SiteID = factor(character(0), levels = as.character(site_ids)),
        SampleSize = integer(0),
        RarefiedRichness = numeric(0)
      )
      next
    }

    n <- seq(1, tot[i], by = 1)
    if (n[length(n)] != tot[i]) n <- c(n, tot[i])

    richness_values <- drop(suppressWarnings(vegan::rarefy(x[i, ], n)))
    sample_sizes <- attr(richness_values, "Subsample")

    output_list[[i]] <- data.frame(
      SiteID = as.factor(site_ids[i]),
      SampleSize = sample_sizes,
      RarefiedRichness = as.numeric(richness_values)
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
