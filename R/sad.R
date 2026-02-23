#' Species abundance distribution (SAD) generators
#'
#' @description
#' Utilities to generate **species abundance distributions** (SADs) for a
#' community of `n_species` and `n_individuals`.
#'
#' In spesim, SADs are used to create a per-species abundance vector (named
#' `A`, `B`, ...), after which individuals are placed in space.
#'
#' @details
#' Most classic SADs are defined as a probability vector
#' \eqn{p_1,\ldots,p_S} over species ranks, then sampled to integer counts via
#' `rmultinom(1, n_individuals, prob = p)`.
#'
#' Some sampling-model SADs (e.g. Poisson--lognormal, Poisson--gamma) are
#' implemented by first simulating latent rates and then drawing counts. Because
#' spesim's downstream spatial simulation typically expects exactly
#' `n_individuals` points, these generators rescale/adjust to sum exactly to
#' `n_individuals`.
#'
#' @section Models:
#' `generate_sad()` supports (case-insensitive):
#' \describe{
#'   \item{`"fisher"`}{Fisher log-series with an explicit dominant species via
#'     [generate_fisher_log_series()].}
#'   \item{`"geometric"`}{Geometric series (niche pre-emption / Motomura).}
#'   \item{`"brokenstick"`}{Broken stick (MacArthur).}
#'   \item{`"zipf"`}{Zipf rank-abundance: \eqn{p_i \propto i^{-a}}.}
#'   \item{`"zipf-mandelbrot"`}{Zipf--Mandelbrot:
#'     \eqn{p_i \propto (i + q)^{-a}}.}
#'   \item{`"lognormal"`}{Lognormal weights + multinomial sampling.}
#'   \item{`"poisson-lognormal"`}{Poisson--lognormal sampling model (rates are
#'     lognormal; counts are Poisson).}
#'   \item{`"poisson-gamma"`}{Poisson--gamma sampling model (rates are Gamma;
#'     counts are Poisson). Closely related to negative-binomial mixtures.}
#'   \item{`"nbd"`}{Negative Binomial Distribution. A flexible model for count
#'     data that can handle high variance (overdispersion).}
#'   \item{`"zsm"`}{Neutral-theory SAD helper. With \code{theta} only, this uses a
#'     theta-only Ewens sampler (Chinese restaurant) to control the
#'     \strong{rank--abundance} curve. If an immigration probability
#'     \eqn{m \in (0,1)} is supplied, spesim uses a simple neutral
#'     \emph{death--birth with immigration} heuristic (Moran style) to make the
#'     SAD more uneven as \code{m} decreases. This is a pragmatic SAD helper and
#'     is not a full neutral dynamics engine (spatial patterns are controlled
#'     separately by the point-process settings). \strong{Important notes:}
#'     (1) The immigration-modulated model runs until convergence is detected
#'     (monitoring the abundance vector) or `max_steps` is reached.
#'     (2) The simulation may produce more or fewer species than `n_species`;
#'     the resulting list is truncated or padded with zeros, and the total
#'     number of individuals is adjusted to match `n_individuals` using a
#'     proportional allocation method, which avoids lumping tail-end abundances
#'     into a single species.}
#'   \item{`"custom"`}{User-supplied SAD via a numeric vector or a function
#'     (`sad` argument).}
#' }
#'
#' @name sad_generators
NULL

# --- internal helpers -------------------------------------------------------

.sad_require_species_labels <- function(n_species) {
  n_species <- as.integer(n_species)
  if (!is.finite(n_species) || n_species < 1L) {
    stop("n_species must be an integer >= 1")
  }
  if (n_species > length(LETTERS)) {
    stop("n_species > 26 is not supported (uses LETTERS for species labels).")
  }
  LETTERS[seq_len(n_species)]
}

.sad_counts_adjust_to_n <- function(counts, n) {
  # Robustly coerce a numeric vector to non-negative integer counts summing to n.
  n <- as.integer(n)
  if (!is.finite(n) || n < 0L) stop("n must be a non-negative integer")

  x <- as.numeric(counts)
  x[!is.finite(x)] <- 0
  x[x < 0] <- 0

  if (length(x) == 0L) {
    return(integer(0))
  }

  if (sum(x) == 0) {
    # uniform seed to avoid all-zeros
    x <- rep(1, length(x))
  }

  # proportional allocation with remainder distribution
  scaled <- x / sum(x) * n
  base <- floor(scaled)
  rem <- n - sum(base)

  if (rem > 0L) {
    # add remainder to species with largest fractional parts
    frac <- scaled - base
    o <- order(frac, decreasing = TRUE)
    base[o[seq_len(rem)]] <- base[o[seq_len(rem)]] + 1L
  } else if (rem < 0L) {
    # remove excess from most abundant species to be less biased
    k <- -rem
    while (k > 0) {
      # find max abundances, break ties by picking first one
      idx_max <- which.max(base)
      if (base[idx_max] > 0) {
        base[idx_max] <- base[idx_max] - 1L
        k <- k - 1L
      } else {
        # This should not happen if sum(base) > n and there are items.
        # If it does, we must take from another cell.
        other_positives <- which(base > 0)
        if (length(other_positives) > 0) {
          base[other_positives[1]] <- base[other_positives[1]] - 1L
          k <- k - 1L
        } else {
          # All are zero, can't remove any more.
          break
        }
      }
    }
  }

  base <- as.integer(base)
  # final guard: if there is still a mismatch, adjust the most abundant species.
  # This is less biased than always adjusting the first.
  diff <- n - sum(base)
  if (diff != 0 && length(base) >= 1L) {
    idx_max <- which.max(base)
    base[idx_max] <- base[idx_max] + diff
  }
  base
}

.sad_prob_normalize <- function(p) {
  p <- as.numeric(p)
  p[!is.finite(p)] <- 0
  p[p < 0] <- 0
  if (length(p) == 0L) return(p)
  s <- sum(p)
  if (!is.finite(s) || s <= 0) {
    return(rep(1 / length(p), length(p)))
  }
  p / s
}

.sad_counts_from_prob <- function(n, prob) {
  prob <- .sad_prob_normalize(prob)
  as.integer(stats::rmultinom(1, size = as.integer(n), prob = prob)[, 1])
}

# --- probability generators ------------------------------------------------

#' Geometric-series (niche pre-emption) SAD probabilities
#'
#' @param n_species Integer (>= 1). Number of species.
#' @param k Numeric in (0, 1). Pre-emption parameter. Larger `k` increases
#'   dominance (more individuals in the top-ranked species).
#'
#' @return Numeric vector of probabilities of length `n_species` summing to 1.
#' @export
generate_geometric_probabilities <- function(n_species, k = 0.5) {
  spp <- .sad_require_species_labels(n_species)
  k <- as.numeric(k)
  if (!is.finite(k) || k <= 0 || k >= 1) stop("k must be in (0, 1)")
  i <- seq_along(spp)
  p <- k * (1 - k)^(i - 1)
  .sad_prob_normalize(p)
}

#' Broken-stick SAD probabilities
#'
#' @description
#' Simulates a broken-stick partition of the unit interval using `runif()`.
#'
#' @param n_species Integer (>= 1). Number of species.
#'
#' @return Numeric vector of probabilities of length `n_species` summing to 1.
#' @export
generate_brokenstick_probabilities <- function(n_species) {
  spp <- .sad_require_species_labels(n_species)
  if (length(spp) == 1L) return(1)
  cuts <- sort(stats::runif(length(spp) - 1L))
  seg <- diff(c(0, cuts, 1))
  p <- sort(seg, decreasing = TRUE)
  .sad_prob_normalize(p)
}

#' Zipf / Zipf--Mandelbrot SAD probabilities
#'
#' @param n_species Integer (>= 1). Number of species.
#' @param exponent Numeric (> 0). Tail exponent `a`.
#' @param q Numeric (>= 0). Mandelbrot offset. `q = 0` gives Zipf.
#'
#' @return Numeric vector of probabilities of length `n_species` summing to 1.
#' @export
generate_zipf_probabilities <- function(n_species, exponent = 1, q = 0) {
  spp <- .sad_require_species_labels(n_species)
  exponent <- as.numeric(exponent)
  q <- as.numeric(q)
  if (!is.finite(exponent) || exponent <= 0) stop("exponent must be > 0")
  if (!is.finite(q) || q < 0) stop("q must be >= 0")
  i <- seq_along(spp)
  p <- (i + q)^(-exponent)
  .sad_prob_normalize(p)
}

#' Lognormal SAD probabilities (weights)
#'
#' @param n_species Integer (>= 1). Number of species.
#' @param meanlog,sdlog Passed to [stats::rlnorm()].
#'
#' @return Numeric vector of probabilities of length `n_species` summing to 1.
#' @export
generate_lognormal_probabilities <- function(n_species, meanlog = 0, sdlog = 1) {
  spp <- .sad_require_species_labels(n_species)
  w <- stats::rlnorm(length(spp), meanlog = meanlog, sdlog = sdlog)
  .sad_prob_normalize(w)
}

# --- count generators -------------------------------------------------------

#' Generate a species abundance distribution (SAD)
#'
#' @description
#' Returns a named vector of integer abundances for `n_species` that sums to
#' `n_individuals`.
#'
#' @param n_species Integer (>= 1). Number of species.
#' @param n_individuals Integer (>= 0). Number of individuals.
#' @param model Character scalar. One of the models described under
#'   \dQuote{Models} in `?sad_generators`.
#' @param sad For `model = "custom"`: either a numeric vector (probabilities or
#'   counts) or a function `function(n_species, n_individuals, ...)` returning a
#'   numeric vector.
#' @param ... Model-specific parameters. Common arguments include:
#'   `k` (geometric), `exponent` and `q` (Zipf / Mandelbrot), `meanlog`, `sdlog`
#'   (lognormal), `shape`, `rate` (Poisson--gamma), `nbd_mu`, `nbd_size` (NBD),
#'   and neutral-theory `theta`, `m`.
#'
#' For `model = "fisher"`, pass `dominant_fraction`, `alpha`, and `x` (see
#' [generate_fisher_log_series()]).
#'
#' @return A named integer vector of length `n_species` (some entries may be 0)
#'   whose names are `LETTERS[1:n_species]`.
#' @export
generate_sad <- function(n_species, n_individuals, model = "fisher", sad = NULL, ...) {
  spp <- .sad_require_species_labels(n_species)
  n_individuals <- as.integer(n_individuals)
  if (!is.finite(n_individuals) || n_individuals < 0L) {
    stop("n_individuals must be an integer >= 0")
  }

  model <- tolower(as.character(model %||% "fisher"))
  model <- gsub("_", "-", model)
  if (model %in% c("broken-stick", "brokenstick")) model <- "brokenstick"
  if (model %in% c("zipfmandelbrot", "zipf-mandelbrot")) model <- "zipf-mandelbrot"
  if (model %in% c("poissonlognormal", "poisson-lognormal")) model <- "poisson-lognormal"
  if (model %in% c("poissongamma", "poisson-gamma")) model <- "poisson-gamma"
  if (model %in% c("negbin", "nbd")) model <- "nbd"

  dots <- list(...)

  out_counts <- switch(model,
    fisher = {
      dominant_fraction <- dots$dominant_fraction %||% dots$DOMINANT_FRACTION
      alpha <- dots$alpha %||% dots$FISHER_ALPHA
      x <- dots$x %||% dots$FISHER_X
      generate_fisher_log_series(
        n_species = length(spp),
        n_individuals = n_individuals,
        dominant_fraction = as.numeric(dominant_fraction %||% 0.3),
        alpha = as.numeric(alpha %||% 3),
        x = as.numeric(x %||% 0.95)
      )
    },

    geometric = {
      k <- dots$k %||% 0.5
      p <- generate_geometric_probabilities(length(spp), k = k)
      counts <- .sad_counts_from_prob(n_individuals, p)
      stats::setNames(counts, spp)
    },

    brokenstick = {
      p <- generate_brokenstick_probabilities(length(spp))
      counts <- .sad_counts_from_prob(n_individuals, p)
      stats::setNames(counts, spp)
    },

    zipf = {
      exponent <- dots$exponent %||% dots$a %||% 1
      p <- generate_zipf_probabilities(length(spp), exponent = exponent, q = 0)
      counts <- .sad_counts_from_prob(n_individuals, p)
      stats::setNames(counts, spp)
    },

    `zipf-mandelbrot` =,
    zipfmandelbrot =,
    zipf_mandelbrot = {
      exponent <- dots$exponent %||% dots$a %||% 1
      q <- dots$q %||% 1
      p <- generate_zipf_probabilities(length(spp), exponent = exponent, q = q)
      counts <- .sad_counts_from_prob(n_individuals, p)
      stats::setNames(counts, spp)
    },

    lognormal = {
      meanlog <- dots$meanlog %||% 0
      sdlog <- dots$sdlog %||% 1
      p <- generate_lognormal_probabilities(length(spp), meanlog = meanlog, sdlog = sdlog)
      counts <- .sad_counts_from_prob(n_individuals, p)
      stats::setNames(counts, spp)
    },

    `poisson-lognormal` =,
    poisson_lognormal =,
    poilog = {
      meanlog <- dots$meanlog %||% 0
      sdlog <- dots$sdlog %||% 1
      lambda <- stats::rlnorm(length(spp), meanlog = meanlog, sdlog = sdlog)
      lambda <- lambda / sum(lambda) * n_individuals
      raw <- stats::rpois(length(spp), lambda)
      counts <- .sad_counts_adjust_to_n(raw, n_individuals)
      stats::setNames(counts, spp)
    },

    `poisson-gamma` =,
    poisson_gamma =,
    poigamma = {
      shape <- as.numeric(dots$shape %||% dots$k %||% 1)
      rate <- as.numeric(dots$rate %||% 1)
      if (!is.finite(shape) || shape <= 0) stop("shape must be > 0")
      if (!is.finite(rate) || rate <= 0) stop("rate must be > 0")
      lambda <- stats::rgamma(length(spp), shape = shape, rate = rate)
      lambda <- lambda / sum(lambda) * n_individuals
      raw <- stats::rpois(length(spp), lambda)
      counts <- .sad_counts_adjust_to_n(raw, n_individuals)
      stats::setNames(counts, spp)
    },

    nbd = {
      mu <- as.numeric(dots$nbd_mu %||% (n_individuals / length(spp)))
      size <- as.numeric(dots$nbd_size %||% 1)
      if (!is.finite(mu) || mu <= 0) stop("nbd_mu must be > 0")
      if (!is.finite(size) || size <= 0) stop("nbd_size must be > 0")
      raw <- stats::rnbinom(length(spp), size = size, mu = mu)
      counts <- .sad_counts_adjust_to_n(raw, n_individuals)
      stats::setNames(counts, spp)
    },

    zsm = {
      theta <- as.numeric(dots$theta %||% 10)
      m <- as.numeric(dots$m %||% NA_real_)

      # NOTE: we keep the "zsm" label for user familiarity, but spesim treats this
      # as a SAD helper.
      # - With theta only: Ewens sampling formula (theta-only), via Chinese restaurant.
      # - With m in (0,1): a simple neutral death--birth with immigration heuristic
      #   to introduce a dispersal-limitation knob on SAD unevenness.
      if (!is.finite(theta) || theta <= 0) stop("theta must be > 0")

      if (is.na(m) || !is.finite(m) || m >= 1) {
        # --- Ewens sampling formula (theta-only) -------------------------------
        counts <- integer(0)
        for (i in seq_len(n_individuals)) {
          if (length(counts) == 0L) {
            counts <- 1L
          } else {
            p_new <- theta / (theta + i - 1)
            if (stats::runif(1) < p_new) {
              counts <- c(counts, 1L)
            } else {
              j <- sample.int(length(counts), size = 1L, prob = counts)
              counts[j] <- counts[j] + 1L
            }
          }
        }
        counts <- sort(counts, decreasing = TRUE)
      } else {
        # --- Immigration-modulated heuristic with convergence check ---
        if (!is.finite(m) || m <= 0 || m > 1) stop("m must be in (0,1] or NA")

        S <- length(spp)
        J <- as.integer(n_individuals)

        # Metacommunity weights
        alpha <- rep(theta / S, S)
        w <- stats::rgamma(S, shape = alpha, rate = 1)
        p <- w / sum(w)

        # Initial local community
        counts <- as.integer(stats::rmultinom(1, size = J, prob = p)[, 1])

        # Moran-style death--birth with immigration until convergence
        max_steps <- as.integer(dots$max_steps %||% 1e6)
        burn_in <- as.integer(dots$burn_in %||% max(1000L, 20L * J))
        check_every <- as.integer(dots$check_every %||% 100L)
        tolerance <- as.numeric(dots$tolerance %||% 1e-4)

        last_counts <- counts

        for (t in 1:max_steps) {
          # death
          d_prob <- counts / sum(counts)
          if (any(is.na(d_prob)) || sum(d_prob) == 0) d_prob <- rep(1 / S, S)
          d <- sample.int(S, size = 1L, prob = d_prob)
          counts[d] <- counts[d] - 1L

          # birth
          if (stats::runif(1) < m) {
            b <- sample.int(S, size = 1L, prob = p)
          } else {
            b_prob <- counts / sum(counts)
            if (any(is.na(b_prob)) || sum(b_prob) == 0) {
              b <- sample.int(S, size = 1L, prob = p)
            } else {
              b <- sample.int(S, size = 1L, prob = b_prob)
            }
          }
          counts[b] <- counts[b] + 1L

          if (t > burn_in && t %% check_every == 0) {
            # Check for convergence: relative change in abundance vector norm
            norm_diff <- sqrt(sum((counts - last_counts)^2)) / J
            if (norm_diff < tolerance) {
              break
            }
            last_counts <- counts
          }
        }

        if (t == max_steps) {
          warning("ZSM simulation reached max_steps without converging.")
        }

        counts <- sort(counts, decreasing = TRUE)
      }

      # Map to a fixed label set of size n_species by truncating.
      # The previous method of lumping tail abundances created artifacts.
      # Truncating and then adjusting to N is a less biased approach.
      if (length(counts) > length(spp)) {
        counts <- counts[seq_len(length(spp))]
      } else if (length(counts) < length(spp)) {
        counts <- c(counts, rep(0L, length(spp) - length(counts)))
      }
      counts <- .sad_counts_adjust_to_n(counts, n_individuals)
      stats::setNames(as.integer(counts), spp)
    },

    custom = {
      if (is.null(sad)) stop("For model='custom', provide `sad` as a numeric vector or function.")

      v <- NULL
      if (is.function(sad)) {
        v <- sad(n_species = length(spp), n_individuals = n_individuals, ...)
      } else {
        v <- sad
      }
      if (!is.numeric(v)) stop("Custom SAD must return/provide a numeric vector.")

      v <- as.numeric(v)
      if (length(v) == length(spp)) {
        # ok
      } else if (length(v) > length(spp)) {
        v <- v[seq_len(length(spp))]
      } else {
        v <- c(v, rep(0, length(spp) - length(v)))
      }

      # Interpret as probabilities if it approximately sums to 1.
      if (isTRUE(abs(sum(v) - 1) < 1e-6)) {
        counts <- .sad_counts_from_prob(n_individuals, v)
      } else {
        # Interpret as weights/counts and adjust to N.
        counts <- .sad_counts_adjust_to_n(v, n_individuals)
      }
      stats::setNames(as.integer(counts), spp)
    },

    stop("Unknown SAD model: '", model, "'.")
  )

  # Ensure full-length output (including zeros) and correct naming.
  out <- integer(length(spp))
  names(out) <- spp
  out[names(out_counts)] <- as.integer(out_counts)
  out <- .sad_counts_adjust_to_n(out, n_individuals)
  names(out) <- spp
  out
}
