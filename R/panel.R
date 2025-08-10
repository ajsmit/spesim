#' Advanced Analysis Panel (rank/occupancy/SAR/decay/rarefaction + optional PPC)
#'
#' @description
#' Build a multi‑panel summary of key ecological diagnostics from a simulation
#' run: rank–abundance (SAD), occupancy–abundance, species–area (accumulation),
#' distance–decay, and per‑quadrat rarefaction curves. If the point‑process
#' packages \pkg{spatstat.geom} and \pkg{spatstat.explore} are available, an
#' additional row of **point–process diagnostics** is appended, showing
#' Ripley’s K (border correction), \eqn{L(r)-r}, and the pair‑correlation
#' function \eqn{g(r)}.
#'
#' @param res A list produced by the main simulator containing at least:
#' \describe{
#'   \item{\code{P}}{Parameter list (used for labels/themes if needed).}
#'   \item{\code{species_dist}}{An \code{sf} POINT layer of individuals with a \code{species} column.}
#'   \item{\code{abund_matrix}}{Site \eqn{\times} species abundance data frame (first column \code{site}).}
#'   \item{\code{site_coords}}{Data frame with \code{site, x, y} for quadrat centroids.}
#'   \item{\code{domain}}{An \code{sf} polygon/multipolygon of the study area (used for PPC window).}
#' }
#'
#' @return A \code{patchwork} object (a \code{ggplot} layout). You can print it
#'   or save it with \code{ggplot2::ggsave()}.
#'
#' @details
#' The optional point–process diagnostics are computed only when both
#' \pkg{spatstat.geom} (for spatial windows and point patterns) and
#' \pkg{spatstat.explore} (for summary functions such as \code{Kest} and
#' \code{pcf}) are installed. The domain is converted to a window (\code{owin})
#' using \code{spatstat.geom::as.owin()} on the \code{sf} polygon; the
#' individuals are converted to a \code{ppp} with that window. We then compute:
#' \itemize{
#'   \item \code{Kest(ppp, correction = "border")} and plot the border‑corrected K.
#'   \item \eqn{L(r) = \sqrt{K(r)/\pi}} and \eqn{L(r)-r} (\eqn{>0} suggests clustering).
#'   \item \code{pcf(ppp)} as an estimate of the pair‑correlation \eqn{g(r)}
#'         (\eqn{>1} suggests clustering; \eqn{<1} inhibition).
#' }
#' If conversion to \code{owin}/\code{ppp} fails, the PPC row is omitted gracefully.
#'
#' @seealso
#' \code{\link{calculate_rank_abundance}}, \code{\link{plot_rank_abundance}},
#' \code{\link{calculate_occupancy_abundance}}, \code{\link{plot_occupancy_abundance}},
#' \code{\link{calculate_species_area}}, \code{\link{plot_species_area}},
#' \code{\link{calculate_distance_decay}}, \code{\link{plot_distance_decay}},
#' \code{\link{calculate_rarefaction}}, \code{\link{plot_rarefaction}}
#'
#' @examples
#' \dontrun{
#' panel <- generate_advanced_panel(results_list)
#' ggplot2::ggsave("advanced_panel.png", panel, width = 12, height = 14, dpi = 300)
#' }
#' @export
generate_advanced_panel <- function(res) {
  # --- Common theme tweaks for smaller panels
  theme_panel <-
    ggplot2::theme(
      text = ggplot2::element_text(size = 11),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9)
    )

  # --- Core five plots -------------------------------------------------
  p_rank <- plot_rank_abundance(calculate_rank_abundance(res$species_dist, res$P)) +
    ggplot2::labs(subtitle = NULL) + theme_panel + ggplot2::theme(legend.position = "bottom")

  p_oa <- plot_occupancy_abundance(calculate_occupancy_abundance(res$abund_matrix)) +
    ggplot2::labs(subtitle = NULL) + theme_panel

  p_sar <- plot_species_area(calculate_species_area(res$abund_matrix)) +
    ggplot2::labs(subtitle = NULL) + theme_panel

  p_decay <- plot_distance_decay(calculate_distance_decay(res$abund_matrix, res$site_coords)) +
    ggplot2::labs(subtitle = NULL) + theme_panel

  p_rare <- plot_rarefaction(calculate_rarefaction(res$abund_matrix)) +
    ggplot2::labs(subtitle = NULL) + theme_panel + ggplot2::guides(color = "none")

  # --- Optional: Point-process diagnostics (K, L-r, g) ----------------
  p_ppc <- NULL
  if (all(vapply(c("spatstat.geom", "spatstat.explore"),
    requireNamespace, logical(1),
    quietly = TRUE
  ))) {
    # Convert domain -> owin (prefer exact polygon if possible)
    W <- try(spatstat.geom::as.owin(sf::st_as_sf(res$domain)), silent = TRUE)
    if (!inherits(W, "try-error")) {
      # Extract point coordinates and build ppp with that window
      xy <- sf::st_coordinates(res$species_dist)
      ppp <- spatstat.geom::ppp(x = xy[, 1], y = xy[, 2], window = W, check = FALSE)

      # Kest and pcf from spatstat.explore (CRAN)
      K <- spatstat.explore::Kest(ppp, correction = "border")
      L <- with(as.data.frame(K), data.frame(r = r, L = sqrt(border / pi)))
      g <- spatstat.explore::pcf(ppp)

      # Quick ggplots for diagnostics
      pK <- ggplot2::ggplot(as.data.frame(K), ggplot2::aes(r, border)) +
        ggplot2::geom_line() +
        ggplot2::labs(title = "Ripley's K (border)", x = "r", y = "K(r)") +
        ggplot2::theme_bw(11)

      pL <- ggplot2::ggplot(L, ggplot2::aes(r, L - r)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::geom_line() +
        ggplot2::labs(title = "L(r) - r", x = "r", y = "L(r) - r") +
        ggplot2::theme_bw(11)

      pg <- ggplot2::ggplot(as.data.frame(g), ggplot2::aes(r, trans)) +
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
        ggplot2::geom_line() +
        ggplot2::labs(title = "Pair correlation g(r)", x = "r", y = "g(r)") +
        ggplot2::theme_bw(11)

      p_ppc <- (pK | pL | pg)
    }
  }

  # --- Assemble layout -------------------------------------------------
  layout <- if (!is.null(p_ppc)) {
    (p_rank | p_oa) / (p_sar | p_decay) / (p_rare | p_ppc)
  } else {
    (p_rank | p_oa) / (p_sar | p_decay) / (p_rare | patchwork::plot_spacer())
  }

  layout + patchwork::plot_annotation(
    title = "Advanced Ecological Analysis Panel",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 18, hjust = 0.5))
  )
}
