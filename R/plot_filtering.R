#' Plot environmental filtering responses (sanity check)
#'
#' @description
#' A teaching/method-testing helper that visualises whether a gradient-responsive
#' species shows higher abundance in quadrats whose environments are near the
#' configured optimum (with tolerance band).
#'
#' This plot is designed to pair with [audit_environmental_filtering()]. It is
#' intentionally simple: it uses per-quadrat mean environments (as computed by
#' [calculate_quadrat_environment()]) and the site x species abundance matrix.
#'
#' @param res A `spesim_result` list returned by [spesim_run()].
#' @param species Species labels to plot. Default `"all"` uses all
#'   gradient-responsive species in `res$P$GRADIENT`.
#' @param add_smooth Logical; add a loess smooth (teaching-friendly). Default
#'   `TRUE`.
#'
#' @return A ggplot object.
#' @export
plot_filtering_response <- function(res, species = "all", add_smooth = TRUE) {
  .data <- NULL
  if (is.null(res$P$GRADIENT) || nrow(res$P$GRADIENT) == 0) {
    stop("No gradient-responsive species found in res$P$GRADIENT")
  }
  if (is.null(res$abund_matrix) || is.null(res$quadrats) || is.null(res$env_gradients)) {
    stop("res must contain abund_matrix, quadrats, and env_gradients")
  }

  G <- res$P$GRADIENT
  spp_all <- as.character(G$species)
  spp <- if (identical(species, "all")) spp_all else intersect(as.character(species), spp_all)
  if (length(spp) == 0) stop("No requested species found in res$P$GRADIENT")

  # per-quadrat env in natural units
  E <- calculate_quadrat_environment(res$env_gradients, res$quadrats, sf::st_crs(res$domain))
  A <- res$abund_matrix

  .opt_to_units <- function(opt, gname) {
    switch(gname,
      "temperature" = opt * 30 - 2,
      "elevation"   = opt * 2000,
      "rainfall"    = opt * 700 + 200,
      NA_real_
    )
  }
  .tol_to_units <- function(tol, gname) {
    switch(gname,
      "temperature" = tol * 30,
      "elevation"   = tol * 2000,
      "rainfall"    = tol * 700,
      NA_real_
    )
  }
  .env_col <- function(gname) {
    switch(gname,
      "temperature" = "temperature_C",
      "elevation"   = "elevation_m",
      "rainfall"    = "rainfall_mm",
      NA_character_
    )
  }
  .env_label <- function(gname) {
    switch(gname,
      "temperature" = "Temperature (deg C)",
      "elevation"   = "Elevation (m)",
      "rainfall"    = "Rainfall (mm)",
      "Environment"
    )
  }

  # Build long data frame: one row per site x species to plot
  rows <- lapply(spp, function(s) {
    gi <- G[G$species == s, , drop = FALSE][1, ]
    gname <- as.character(gi$gradient)
    env_col <- .env_col(gname)
    if (is.na(env_col) || !(s %in% colnames(A))) return(NULL)

    df <- dplyr::left_join(E, A[, c("site", s), drop = FALSE], by = c("site" = "site"))
    names(df)[names(df) == s] <- "abundance"

    df$species <- s
    df$gradient <- gname
    df$env <- df[[env_col]]

    opt_u <- .opt_to_units(as.numeric(gi$optimum), gname)
    tol_u <- .tol_to_units(as.numeric(gi$tol), gname)

    df$optimum <- opt_u
    df$tol <- tol_u
    df$env_label <- .env_label(gname)

    df[, c("site", "species", "gradient", "env_label", "env", "abundance", "optimum", "tol")]
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) stop("No plot-able species found")

  dat <- do.call(rbind, rows)
  dat <- dat[is.finite(dat$env) & is.finite(dat$abundance), , drop = FALSE]

  # Plot: abundance vs environment, with optimum line and +/- tol band
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$env, y = .data$abundance)) +
    ggplot2::geom_point(alpha = 0.7, size = 2) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = .data$optimum), linetype = "dashed", linewidth = 0.6) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = .data$optimum - .data$tol, xmax = .data$optimum + .data$tol, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey85",
      alpha = 0.25
    ) +
    ggplot2::facet_wrap(~ species + gradient, scales = "free_x") +
    ggplot2::labs(
      x = "Environment (per-quadrat mean)",
      y = "Quadrat abundance",
      title = "Environmental filtering sanity check",
      subtitle = "Dashed line = configured optimum; grey band = +/- tolerance (in natural units)"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  if (isTRUE(add_smooth)) {
    p <- p + ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 0.7)
  }

  p
}
