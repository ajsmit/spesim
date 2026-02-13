test_that("place_quadrats_transect respects compass bearing", {
  skip_if_not_installed("sf")

  set.seed(1)
  dom <- spesim::create_sampling_domain()
  qs <- c(1.2, 1.2)

  angle <- 30
  quads <- spesim::place_quadrats_transect(
    domain = dom,
    n_transects = 1,
    n_quadrats_per_transect = 10,
    quadrat_size = qs,
    angle = angle
  )

  expect_true(inherits(quads, "sf"))
  expect_gt(nrow(quads), 1)

  ctr <- sf::st_coordinates(sf::st_centroid(quads))

  # Estimate transect direction from the principal axis of quadrat centroids.
  pc <- stats::prcomp(ctr, center = TRUE, scale. = FALSE)
  v <- pc$rotation[, 1]

  bearing <- (atan2(v[1], v[2]) * 180 / pi) %% 360
  bearing2 <- (bearing + 180) %% 360

  # Compare modulo 180 degrees (a line has two directions).
  diff1 <- abs(((bearing - angle + 180) %% 360) - 180)
  diff2 <- abs(((bearing2 - angle + 180) %% 360) - 180)
  diff <- min(diff1, diff2)

  expect_lt(diff, 1) # within 1 degree
})
