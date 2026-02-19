test_that("coastline domain type emits linear positions and along-path decay runs", {
  init <- system.file("examples/spesim_init_basic.txt", package = "spesim")
  P <- spesim::load_config(init)
  P$N_INDIVIDUALS <- 300
  P$N_SPECIES <- 8
  P$DOMAIN_TYPE <- "coastline"
  P$LINEAR_AXIS <- "x"
  P$LINEAR_WRAP <- TRUE
  P$LINEAR_JITTER_SD <- 0.05
  P$MODEL_FAMILY <- "neutral_hubbell_like"
  P$DISPERSAL_DIRECTION_BIAS <- 0.5
  P$INTERACTION_RADIUS <- 0
  P$INTERACTIONS_FILE <- NULL
  P$INTERACTIONS_EDGELIST <- NULL

  res <- suppressWarnings(spesim::spesim_run(
    P,
    seed = 99,
    write_outputs = FALSE,
    quiet = TRUE,
    interactions_print = FALSE
  ))

  expect_true("linear_pos" %in% names(res$species_dist))
  expect_true("linear_pos" %in% names(res$site_coords))
  expect_equal(attr(res$site_coords, "distance_metric"), "along_path")

  dd <- suppressWarnings(spesim::calculate_distance_decay(res$abund_matrix, res$site_coords))
  expect_true(nrow(dd) > 0)
  expect_true(all(dd$Distance >= 0))
})

test_that("load_config reads external ENV_COVARIATES_FILE", {
  td <- tempdir()
  cov_file <- file.path(td, "covs.csv")
  write.csv(data.frame(
    x = c(0, 1, 0, 1),
    y = c(0, 0, 1, 1),
    temperature_C = c(10, 12, 14, 16),
    elevation_m = c(100, 200, 300, 400),
    rainfall_mm = c(500, 550, 600, 650)
  ), cov_file, row.names = FALSE)

  init_file <- file.path(td, "init_cov.txt")
  writeLines(c(
    "SEED = 7",
    "N_INDIVIDUALS = 100",
    "N_SPECIES = 5",
    paste0('ENV_COVARIATES_FILE = "', cov_file, '"')
  ), init_file)

  P <- spesim::load_config(init_file)
  expect_true(is.data.frame(P$ENV_COVARIATES))
  expect_true(all(c("x", "y", "temperature_C") %in% names(P$ENV_COVARIATES)))
})
