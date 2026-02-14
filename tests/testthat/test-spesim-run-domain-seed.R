test_that("spesim_run(seed=) makes the whole run reproducible (domain + outputs)", {
  skip_if_not_installed("sf")

  init <- system.file("examples/spesim_init_basic.txt", package = "spesim")
  P <- spesim::load_config(init)

  # keep it fast
  P$N_INDIVIDUALS <- 400L
  P$N_SPECIES <- 6L
  P$N_QUADRATS <- 10L
  P$SAMPLING_RESOLUTION <- 20L
  P$ADVANCED_ANALYSIS <- FALSE

  r1 <- suppressWarnings(spesim::spesim_run(P, seed = 202, write_outputs = FALSE, quiet = TRUE,
    interactions_print = FALSE))
  r2 <- suppressWarnings(spesim::spesim_run(P, seed = 202, write_outputs = FALSE, quiet = TRUE,
    interactions_print = FALSE))

  expect_equal(sf::st_coordinates(r1$domain), sf::st_coordinates(r2$domain))
  expect_equal(r1$abund_matrix, r2$abund_matrix)
})

test_that("spesim_run() with an init-file respects seed= override", {
  skip_if_not_installed("sf")

  init <- system.file("examples/spesim_init_basic.txt", package = "spesim")

  r <- suppressWarnings(spesim::spesim_run(init, seed = 555, write_outputs = FALSE, quiet = TRUE,
    interactions_print = FALSE))

  expect_equal(r$P$SEED, 555L)
})
