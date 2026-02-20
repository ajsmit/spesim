#' Spatially-explicit neutral community simulation with convergence checks
#'
#' @description
#' Simulates a spatially-explicit neutral or hybrid community model, running
#' until the community composition converges or a maximum number of steps is
#' reached. This is an internal function called by
#' `generate_heterogeneous_distribution` when the model family is set to
#' `"neutral_hubbell_like"` or `"hybrid"`.
#'
#' @details
#' The simulation follows a Moran-style process of death, followed by recruitment.
#' Recruitment can be from a local parent (with dispersal limitation) or from an
#' external metacommunity.
#'
#' **Convergence:**
#' Instead of running for a fixed number of steps, this implementation checks
#' for convergence. At intervals specified by `NEUTRAL_CONVERGENCE_INTERVAL`,
#' it calculates the Bray-Curtis dissimilarity between the current community's
#' species abundance vector and the vector from the previous check. If this
#' dissimilarity is below `NEUTRAL_CONVERGENCE_THRESHOLD` for a number of
#' consecutive checks equal to `NEUTRAL_CONVERGENCE_PATIENCE`, the community
#' is considered to have reached a stationary state, and the simulation stops.
#' A `NEUTRAL_MAX_STEPS` parameter prevents infinite loops.
#'
#' @param poly An `sf` polygon defining the sampling domain.
#' @param P A fully materialized parameter list from `load_config()`. It must
#'   also contain convergence parameters: `NEUTRAL_MAX_STEPS`,
#'   `NEUTRAL_CONVERGENCE_INTERVAL`, `NEUTRAL_CONVERGENCE_THRESHOLD`,
#'   `NEUTRAL_CONVERGENCE_PATIENCE`.
#' @return An `sf` POINT layer with species and environmental data.
#'
#' @keywords internal
spesim_simulate_neutral_recruitment <- function(poly, P) {
  # This function uses helpers from R/simulation_helpers.R
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  J <- as.integer(P$N_INDIVIDUALS %||% 0L)
  S <- as.integer(P$N_SPECIES %||% 1L)
  spp <- LETTERS[seq_len(max(1L, S))]
  crs_dom <- sf::st_crs(poly)

  if (!is.finite(J) || J < 0L) J <- 0L
  if (J == 0L) {
    return(sf::st_sf(species = character(0), geometry = sf::st_sfc(crs = crs_dom)))
  }

  family <- tolower(as.character(P$MODEL_FAMILY %||% "manual"))
  family <- gsub("-", "_", family)
  hybrid_mode <- identical(family, "hybrid")
  linear_st <- .linear_domain_state(poly, P)

  # --- Model parameters ---
  m <- as.numeric(P$NEUTRAL_M %||% 0.1)
  nu <- as.numeric(P$NEUTRAL_NU %||% 0.0)
  m <- min(1, max(0, ifelse(is.finite(m), m, 0.1)))
  nu <- min(1, max(0, ifelse(is.finite(nu), nu, 0.0)))

  kernel <- tolower(as.character(P$DISPERSAL_KERNEL %||% "gaussian"))
  if (!kernel %in% c("gaussian", "exponential", "power_law")) kernel <- "gaussian"
  scale <- as.numeric(P$DISPERSAL_SCALE %||% 0.5)
  if (!is.finite(scale) || scale <= 0) scale <- 0.5
  alpha <- as.numeric(P$DISPERSAL_ALPHA %||% 2)
  if (!is.finite(alpha) || alpha <= 0) alpha <- 2
  dir_bias <- as.numeric(P$DISPERSAL_DIRECTION_BIAS %||% 0)
  dir_bias <- max(-1, min(1, ifelse(is.finite(dir_bias), dir_bias, 0)))
  lambda_env <- as.numeric(P$HYBRID_ENV_WEIGHT %||% 1)
  if (!is.finite(lambda_env) || lambda_env < 0) lambda_env <- 1

  # --- Convergence parameters ---
  max_steps <- as.integer(P$NEUTRAL_MAX_STEPS %||% (J * 50))
  conv_interval <- as.integer(P$NEUTRAL_CONVERGENCE_INTERVAL %||% J)
  conv_threshold <- as.numeric(P$NEUTRAL_CONVERGENCE_THRESHOLD %||% 0.01)
  conv_patience <- as.integer(P$NEUTRAL_CONVERGENCE_PATIENCE %||% 3)

  # --- Metacommunity ---
  meta_model <- as.character(P$NEUTRAL_META_MODEL %||% P$SAD_MODEL %||% "zsm")
  meta_counts <- generate_sad(
    n_species = length(spp), n_individuals = max(J, length(spp)),
    model = meta_model, dominant_fraction = P$DOMINANT_FRACTION,
    alpha = P$FISHER_ALPHA, x = P$FISHER_X, k = P$GEOMETRIC_K,
    exponent = P$ZIPF_EXPONENT, q = P$ZIPF_Q, meanlog = P$LOGNORMAL_MEANLOG,
    sdlog = P$LOGNORMAL_SDLOG, shape = P$POIGAMMA_SHAPE,
    rate = P$POIGAMMA_RATE, theta = P$ZSM_THETA, m = P$ZSM_M, sad = P$SAD_VECTOR
  )
  meta_prob <- as.numeric(meta_counts)
  if (!all(is.finite(meta_prob)) || sum(meta_prob) <= 0) {
    meta_prob <- rep(1 / length(spp), length(spp))
  } else {
    meta_prob <- meta_prob / sum(meta_prob)
  }

  # --- Environment setup ---
  drv <- unique(c(as.character(P$ENV_DRIVERS %||% character(0)),
                  as.character(P$GRADIENT$gradient %||% character(0))))
  env_grid <- create_environmental_gradients(
    poly, P$SAMPLING_RESOLUTION, P$ENVIRONMENTAL_NOISE,
    covariates = P$ENV_COVARIATES %||% NULL, drivers = drv
  )
  env_sf <- sf::st_as_sf(env_grid, coords = c("x", "y"), crs = crs_dom)

  # Grid helpers for fast lookup
  gx <- sort(unique(env_grid$x))
  gy <- sort(unique(env_grid$y))
  ix <- match(env_grid$x, gx)
  iy <- match(env_grid$y, gy)
  .to_grid_index <- function(v, grid_vals) {
    if (length(grid_vals) <= 1L) return(1L)
    mids <- (grid_vals[-1L] + grid_vals[-length(grid_vals)]) / 2
    findInterval(v, mids) + 1L
  }
  .build_grid <- function(colname) {
    mat <- matrix(NA_real_, nrow = length(gy), ncol = length(gx))
    if(colname %in% names(env_grid)) {
        mat[cbind(iy, ix)] <- as.numeric(env_grid[[colname]])
    }
    mat
  }
  temp_mat <- .build_grid("temperature")
  elev_mat <- .build_grid("elevation")
  rain_mat <- .build_grid("rainfall")

  grad_lookup <- list()
  if (hybrid_mode && !is.null(P$GRADIENT) && NROW(P$GRADIENT) > 0) {
    for (i in seq_len(NROW(P$GRADIENT))) {
      row <- P$GRADIENT[i, , drop = FALSE]
      grad_lookup[[as.character(row$species)]] <- list(
        gradient = as.character(row$gradient),
        optimum = as.numeric(row$optimum),
        tol = as.numeric(row$tol)
      )
    }
  }

  # --- Dispersal & Environment helpers ---
  .sample_immigrant_species <- function() sample(spp, size = 1L, prob = meta_prob)
  .step_magnitude <- function() {
    if (kernel == "gaussian") {
      abs(stats::rnorm(1L, 0, scale))
    } else if (kernel == "exponential") {
      stats::rexp(1L, rate = 1 / scale)
    } else { # power_law
      u <- stats::runif(1L)
      scale * ((1 - u)^(-1 / alpha) - 1)
    }
  }
  .propose_displacement <- function(parent_xy) {
    ang <- stats::runif(1L, 0, 2 * pi)
    r <- .step_magnitude()
    dx <- r * cos(ang); dy <- r * sin(ang)
    c(parent_xy[1] + dx, parent_xy[2] + dy)
  }
  .propose_linear_pos <- function(parent_pos) {
    mag <- .step_magnitude()
    axis_range <- if (identical(linear_st$axis, "x")) linear_st$xmax - linear_st$xmin else linear_st$ymax - linear_st$ymin
    if (!is.finite(axis_range) || axis_range <= 0) axis_range <- 1
    step <- mag / axis_range
    p_plus <- (1 + dir_bias) / 2
    sgn <- if (stats::runif(1) < p_plus) 1 else -1
    out <- parent_pos + sgn * step
    if (isTRUE(linear_st$wrap)) out %% 1 else max(0, min(1, out))
  }
  .env_accept_prob <- function(sp, xy) {
    if (!hybrid_mode || lambda_env == 0 || is.null(grad_lookup[[sp]])) return(1)
    g <- grad_lookup[[sp]]
    env_col <- .resolve_env_col(g$gradient, names(env_grid))
    if (is.na(env_col)) return(1)
    ixx <- .to_grid_index(as.numeric(xy[1]), gx)
    iyy <- .to_grid_index(as.numeric(xy[2]), gy)
    val <- if (env_col == "temperature") {
      temp_mat[iyy, ixx]
    } else if (env_col == "elevation") {
      elev_mat[iyy, ixx]
    } else if (env_col == "rainfall") {
      rain_mat[iyy, ixx]
    } else {
      pt <- sf::st_as_sf(data.frame(x=xy[1], y=xy[2]), coords=c("x","y"), crs=crs_dom)
      j <- sf::st_nearest_feature(pt, env_sf)
      as.numeric(env_sf[[env_col]][j])
    }
    if (!is.finite(val) || !is.finite(g$tol) || g$tol <= 0) return(1)
    w <- exp(-((val - g$optimum)^2) / (2 * g$tol^2))
    max(0, min(1, w^lambda_env))
  }

  # --- Simulation State ---
  xy <- matrix(NA_real_, nrow = J, ncol = 2)
  linear_pos <- rep(NA_real_, J)
  spp_out <- rep(NA_character_, J)

  # Initial community from metacommunity
  for(i in 1:J) {
      spp_out[i] <- .sample_immigrant_species()
      if(isTRUE(linear_st$enabled)) {
          linear_pos[i] <- stats::runif(1)
          xy[i,] <- .linear_pos_to_xy(linear_pos[i], linear_st, poly)[1,]
      } else {
          xy[i,] <- .sample_uniform_in_domain(1L, poly)[1,]
      }
  }

  patience_counter <- 0
  last_community_abund <- table(factor(spp_out, levels = spp))

  # --- Main Simulation Loop ---
  for (t in seq_len(max_steps)) {
    # 1. Death of one random individual
    dead_idx <- sample.int(J, 1L)
    
    # 2. Birth to replace it
    is_immigrant <- (stats::runif(1) < m)
    if (is_immigrant) {
      parent_idx <- NA_integer_
      birth_sp <- .sample_immigrant_species()
    } else {
      parent_idx <- sample.int(J, 1L) 
      birth_sp <- spp_out[parent_idx]
      if (stats::runif(1L) < nu) { # speciation
        birth_sp <- .sample_immigrant_species()
      }
    }

    # Propose new location and accept/reject
    accepted <- FALSE
    for (att in seq_len(40L)) {
      if(is_immigrant || is.na(parent_idx)) {
        new_xy <- if(isTRUE(linear_st$enabled)) .linear_pos_to_xy(stats::runif(1), linear_st, poly) else .sample_uniform_in_domain(1L, poly)
      } else {
        new_xy <- if(isTRUE(linear_st$enabled)) .linear_pos_to_xy(.propose_linear_pos(linear_pos[parent_idx]), linear_st, poly) else matrix(.propose_displacement(xy[parent_idx,]), ncol=2)
      }
      
      if (NROW(new_xy) > 0 && .in_domain(new_xy, poly)) {
        p_acc <- .env_accept_prob(birth_sp, new_xy)
        if (stats::runif(1L) <= p_acc) {
            xy[dead_idx, ] <- new_xy[1,]
            spp_out[dead_idx] <- birth_sp
            if(isTRUE(linear_st$enabled)) linear_pos[dead_idx] <- .xy_to_linear_pos(new_xy, linear_st)
            accepted <- TRUE
            break
        }
      }
    }
    if (!accepted) { # Fallback: if no valid spot found, keep species but location is unchanged.
      spp_out[dead_idx] <- birth_sp
    }

    # 3. Convergence Check
    if (t %% conv_interval == 0) {
      current_abund <- table(factor(spp_out, levels = spp))
      diss <- vegan::vegdist(rbind(current_abund, last_community_abund), method = "bray")[1]
      if (is.finite(diss) && diss < conv_threshold) {
        patience_counter <- patience_counter + 1
      } else {
        patience_counter <- 0
      }
      
      if (patience_counter >= conv_patience) {
        if(!isTRUE(P$QUIET)) message(sprintf("Neutral simulation converged after %d steps.", t))
        break
      }
      last_community_abund <- current_abund
    }
  }
  
  if (t == max_steps && !isTRUE(P$QUIET)) {
    warning("Neutral simulation reached max_steps without converging.", call. = FALSE)
  }

  # --- Finalize and return sf object ---
  out <- sf::st_sf(species = spp_out, geometry = .as_sfc_points(xy, crs_dom))
  if (isTRUE(linear_st$enabled)) {
    out$linear_pos <- linear_pos
  }
  sf::st_join(out, env_sf, join = sf::st_nearest_feature)
}
