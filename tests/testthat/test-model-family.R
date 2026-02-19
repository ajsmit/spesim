test_that("model-family neutral_hubbell_like runs and returns expected size", {
  init <- system.file("examples/spesim_init_basic.txt", package = "spesim")
  P <- spesim::load_config(init)
  P$N_INDIVIDUALS <- 400
  P$N_SPECIES <- 12
  P$MODEL_FAMILY <- "neutral_hubbell_like"
  P$NEUTRAL_M <- 0.12
  P$NEUTRAL_NU <- 0.01
  P$DISPERSAL_KERNEL <- "gaussian"
  P$DISPERSAL_SCALE <- 0.5
  P$INTERACTION_RADIUS <- 0
  P$INTERACTIONS_FILE <- NULL
  P$INTERACTIONS_EDGELIST <- NULL

  res <- suppressWarnings(spesim::spesim_run(
    P,
    seed = 123,
    write_outputs = FALSE,
    quiet = TRUE,
    interactions_print = FALSE
  ))

  expect_s3_class(res$species_dist, "sf")
  expect_equal(nrow(res$species_dist), P$N_INDIVIDUALS)
  expect_true(all(res$species_dist$species %in% LETTERS[seq_len(P$N_SPECIES)]))
})

test_that("model-family hybrid executes with gradient sorting enabled", {
  init <- system.file("examples/spesim_init_basic.txt", package = "spesim")
  P <- spesim::load_config(init)
  P$N_INDIVIDUALS <- 350
  P$N_SPECIES <- 10
  P$MODEL_FAMILY <- "hybrid"
  P$NEUTRAL_M <- 0.08
  P$DISPERSAL_KERNEL <- "exponential"
  P$DISPERSAL_SCALE <- 0.4
  P$HYBRID_ENV_WEIGHT <- 1.2
  P$GRADIENT_SPECIES <- c("A", "B", "C")
  P$GRADIENT_ASSIGNMENTS <- c("temperature", "elevation", "rainfall")
  P$GRADIENT_OPTIMA <- c(A = 0.2, B = 0.5, C = 0.75)
  P$GRADIENT_TOLERANCE <- c(A = 0.15, B = 0.12, C = 0.18)
  P$INTERACTION_RADIUS <- 0
  P$INTERACTIONS_FILE <- NULL
  P$INTERACTIONS_EDGELIST <- NULL

  res <- suppressWarnings(spesim::spesim_run(
    P,
    seed = 321,
    write_outputs = FALSE,
    quiet = TRUE,
    interactions_print = FALSE
  ))

  expect_equal(nrow(res$species_dist), P$N_INDIVIDUALS)
  expect_true(all(c("temperature", "elevation", "rainfall") %in% names(res$species_dist)))
  expect_true(nrow(res$abund_matrix) > 0)
})
