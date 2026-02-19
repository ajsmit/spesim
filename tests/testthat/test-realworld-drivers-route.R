test_that("arbitrary ENV_DRIVERS and GRADIENT_ASSIGNMENTS run", {
  init <- tempfile(fileext = ".txt")
  writeLines(c(
    "SEED = 77",
    "N_INDIVIDUALS = 200",
    "N_SPECIES = 8",
    "SAMPLING_SCHEME = random",
    "N_QUADRATS = 20",
    "QUADRAT_SIZE_OPTION = medium"
  ), init)
  P <- spesim::load_config(init)
  P$N_INDIVIDUALS <- 250
  P$N_SPECIES <- 8
  P$ENV_DRIVERS <- c("sst", "salinity", "nitrate")
  P$GRADIENT_SPECIES <- c("A", "B", "C")
  P$GRADIENT_ASSIGNMENTS <- c("sst", "salinity", "nitrate")
  P$GRADIENT_OPTIMA <- c(sst = 0.2, salinity = 0.6, nitrate = 0.4)
  P$GRADIENT_TOLERANCE <- c(sst = 0.2, salinity = 0.25, nitrate = 0.2)
  P$INTERACTION_RADIUS <- 0
  P$INTERACTIONS_FILE <- NULL
  P$INTERACTIONS_EDGELIST <- NULL

  res <- suppressWarnings(spesim::spesim_run(P, write_outputs = FALSE, seed = 111, quiet = TRUE, interactions_print = FALSE))

  expect_true(all(c("sst", "salinity", "nitrate") %in% names(res$env_gradients)))
  expect_true(all(c("sst", "salinity", "nitrate") %in% names(res$species_dist)))
  expect_true(all(c("sst", "salinity", "nitrate") %in% names(res$site_env)))
})

test_that("route sampling places quadrats for linearized domain", {
  init <- tempfile(fileext = ".txt")
  writeLines(c(
    "SEED = 77",
    "N_INDIVIDUALS = 180",
    "N_SPECIES = 6",
    "SAMPLING_SCHEME = route",
    "N_QUADRATS = 6",
    "QUADRAT_SIZE_OPTION = medium"
  ), init)
  P <- spesim::load_config(init)
  P$N_INDIVIDUALS <- 200
  P$N_SPECIES <- 6
  P$DOMAIN_TYPE <- "coastline"
  P$SAMPLING_SCHEME <- "route"
  P$ROUTE_QUADRAT_MODE <- "specified"
  P$ROUTE_POSITIONS <- c(0.1, 0.2, 0.4, 0.6, 0.8, 0.9)
  P$N_QUADRATS <- length(P$ROUTE_POSITIONS)
  P$INTERACTION_RADIUS <- 0
  P$INTERACTIONS_FILE <- NULL
  P$INTERACTIONS_EDGELIST <- NULL

  res <- suppressWarnings(spesim::spesim_run(P, write_outputs = FALSE, seed = 222, quiet = TRUE, interactions_print = FALSE))

  expect_s3_class(res$quadrats, "sf")
  expect_true(nrow(res$quadrats) >= 1)
  expect_equal(nrow(res$site_coords), nrow(res$quadrats))
  expect_true("linear_pos" %in% names(res$site_coords))
})
