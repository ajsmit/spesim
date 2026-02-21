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
#' **Performance:**
#' All polygon domain checks use a compiled C++ ray-casting point-in-polygon
#' test (see `pip_cpp`), with the polygon ring coordinates extracted once and
#' cached for the duration of the simulation. The initial community placement
#' is performed with a single vectorised rejection-sample call rather than
#' J individual calls, and the same pre-built tester is threaded through every
#' candidate-location filter in the main loop.
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
  # Precompute midpoints once; reused by the vectorised env-acceptance step.
  mids_gx <- if (length(gx) <= 1L) numeric(0L) else (gx[-1L] + gx[-length(gx)]) / 2
  mids_gy <- if (length(gy) <= 1L) numeric(0L) else (gy[-1L] + gy[-length(gy)]) / 2
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
      g_name <- as.character(row$gradient)
      grad_lookup[[as.character(row$species)]] <- list(
        gradient = g_name,
        env_col  = .resolve_env_col(g_name, names(env_grid)),
        optimum  = as.numeric(row$optimum),
        tol      = as.numeric(row$tol)
      )
    }
  }

  # --- Dispersal & Environment helpers ---
  .sample_immigrant_species <- function() sample(spp, size = 1L, prob = meta_prob)
  # Vectorised step magnitudes: one RNG call for n candidates instead of n×1.
  .step_magnitude_batch <- function(n) {
    if (kernel == "gaussian") {
      abs(stats::rnorm(n, 0, scale))
    } else if (kernel == "exponential") {
      stats::rexp(n, rate = 1 / scale)
    } else { # power_law
      u <- stats::runif(n)
      scale * ((1 - u)^(-1 / alpha) - 1)
    }
  }
  # Returns an n×2 matrix of proposed (x, y) positions from a single parent.
  .propose_displacement_batch <- function(parent_xy, n) {
    angs <- stats::runif(n, 0, 2 * pi)
    r    <- .step_magnitude_batch(n)
    cbind(parent_xy[1] + r * cos(angs),
          parent_xy[2] + r * sin(angs))
  }
  # Returns a length-n vector of proposed linear positions from a single parent.
  .propose_linear_pos_batch <- function(parent_pos, n) {
    mag        <- .step_magnitude_batch(n)
    axis_range <- if (identical(linear_st$axis, "x")) linear_st$xmax - linear_st$xmin else linear_st$ymax - linear_st$ymin
    if (!is.finite(axis_range) || axis_range <= 0) axis_range <- 1
    step  <- mag / axis_range
    p_plus <- (1 + dir_bias) / 2
    sgn   <- ifelse(stats::runif(n) < p_plus, 1, -1)
    out   <- parent_pos + sgn * step
    if (isTRUE(linear_st$wrap)) out %% 1 else pmax(0, pmin(1, out))
  }
  # Vectorised over all candidate points at once: avoids per-point closure
  # overhead, repeated .resolve_env_col calls, and repeated mids recomputation.
  .env_accept_prob_vec <- function(sp, xy_mat) {
    n <- nrow(xy_mat)
    if (!hybrid_mode || lambda_env == 0 || is.null(grad_lookup[[sp]])) return(rep(1, n))
    g <- grad_lookup[[sp]]
    env_col <- g$env_col   # resolved once at setup, not per-call
    if (is.na(env_col)) return(rep(1, n))
    ixx <- if (length(mids_gx) == 0L) rep(1L, n) else findInterval(xy_mat[, 1], mids_gx) + 1L
    iyy <- if (length(mids_gy) == 0L) rep(1L, n) else findInterval(xy_mat[, 2], mids_gy) + 1L
    val <- if (env_col == "temperature") {
      temp_mat[cbind(iyy, ixx)]
    } else if (env_col == "elevation") {
      elev_mat[cbind(iyy, ixx)]
    } else if (env_col == "rainfall") {
      rain_mat[cbind(iyy, ixx)]
    } else {
      pts <- sf::st_as_sf(data.frame(x = xy_mat[, 1], y = xy_mat[, 2]),
                          coords = c("x", "y"), crs = crs_dom)
      j <- sf::st_nearest_feature(pts, env_sf)
      as.numeric(env_sf[[env_col]][j])
    }
    w <- rep(1, n)
    fin <- is.finite(val) & is.finite(g$tol) & g$tol > 0
    if (any(fin)) w[fin] <- exp(-((val[fin] - g$optimum)^2) / (2 * g$tol^2))
    pmax(0, pmin(1, w^lambda_env))
  }

  # Pre-build PIP tester once from the fixed polygon; shared by all domain
  # checks in the initialisation and simulation loop.
  pip_fn <- .build_pip_tester(poly)

  # --- Simulation State ---
  # Vectorised initialisation: sample all J positions in one call rather than
  # J individual calls to .sample_uniform_in_domain(1, poly).
  spp_out <- sample(spp, J, replace = TRUE, prob = meta_prob)
  linear_pos <- rep(NA_real_, J)
  if (isTRUE(linear_st$enabled)) {
    linear_pos <- stats::runif(J)
    xy <- .linear_pos_to_xy(linear_pos, linear_st, poly, pip_fn = pip_fn)
  } else {
    xy <- .sample_uniform_in_domain(J, poly, pip_fn = pip_fn)
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
    batch_size <- 40 # Propose 40 candidates at once

    # Generate a batch of candidate locations
    if (is_immigrant || is.na(parent_idx)) {
      cand_xy <- if(isTRUE(linear_st$enabled)) .linear_pos_to_xy(stats::runif(batch_size), linear_st, poly, pip_fn = pip_fn) else .sample_uniform_in_domain(batch_size, poly, pip_fn = pip_fn)
    } else {
      if(isTRUE(linear_st$enabled)) {
        cand_pos <- .propose_linear_pos_batch(linear_pos[parent_idx], batch_size)
        cand_xy <- .linear_pos_to_xy(cand_pos, linear_st, poly, pip_fn = pip_fn)
      } else {
        cand_xy <- .propose_displacement_batch(xy[parent_idx,], batch_size)
      }
    }

    # Filter candidates to those inside the domain
    if (NROW(cand_xy) > 0) {
      in_domain_idx <- which(pip_fn(cand_xy))
      cand_xy <- cand_xy[in_domain_idx, , drop = FALSE]

      if (NROW(cand_xy) > 0) {
        # Check environmental acceptance for valid candidates
        p_acc_vec <- .env_accept_prob_vec(birth_sp, cand_xy)
        
        successful_idx <- which(stats::runif(nrow(cand_xy)) <= p_acc_vec)
        
        if (length(successful_idx) > 0) {
            # Select one of the successful candidates
            chosen_idx <- if (length(successful_idx) == 1) successful_idx else sample(successful_idx, 1)
            new_xy <- cand_xy[chosen_idx, , drop = FALSE]

            xy[dead_idx, ] <- new_xy[1,]
            spp_out[dead_idx] <- birth_sp
            if(isTRUE(linear_st$enabled)) linear_pos[dead_idx] <- .xy_to_linear_pos(new_xy, linear_st)
            accepted <- TRUE
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
