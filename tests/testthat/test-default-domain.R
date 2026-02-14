test_that("spesim_run() default domain varies across runs when seed is NULL", {
  skip_if_not_installed("sf")

  init <- system.file("examples/spesim_init_basic.txt", package = "spesim")
  P <- spesim::load_config(init)

  # keep it fast
  P$N_INDIVIDUALS <- 300L
  P$N_SPECIES <- 6L
  P$N_QUADRATS <- 8L
  P$SAMPLING_RESOLUTION <- 20L
  P$ADVANCED_ANALYSIS <- FALSE

  res1 <- spesim::spesim_run(P, write_outputs = FALSE, quiet = TRUE)
  res2 <- spesim::spesim_run(P, write_outputs = FALSE, quiet = TRUE)

  c1 <- sf::st_coordinates(res1$domain)[, c("X", "Y"), drop = FALSE]
  c2 <- sf::st_coordinates(res2$domain)[, c("X", "Y"), drop = FALSE]

  expect_false(isTRUE(all.equal(c1, c2)))
})

test_that("spesim_run() default domain is reproducible when seed is supplied", {
  skip_if_not_installed("sf")

  init <- system.file("examples/spesim_init_basic.txt", package = "spesim")
  P <- spesim::load_config(init)

  # keep it fast
  P$N_INDIVIDUALS <- 300L
  P$N_SPECIES <- 6L
  P$N_QUADRATS <- 8L
  P$SAMPLING_RESOLUTION <- 20L
  P$ADVANCED_ANALYSIS <- FALSE

  res1 <- spesim::spesim_run(P, write_outputs = FALSE, quiet = TRUE, seed = 123)
  res2 <- spesim::spesim_run(P, write_outputs = FALSE, quiet = TRUE, seed = 123)

  c1 <- sf::st_coordinates(res1$domain)[, c("X", "Y"), drop = FALSE]
  c2 <- sf::st_coordinates(res2$domain)[, c("X", "Y"), drop = FALSE]

  expect_equal(c1, c2)
})
