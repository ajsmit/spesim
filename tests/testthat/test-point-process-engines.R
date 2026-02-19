test_that("rthomas_fast returns points inside polygon", {
  skip_if_not_installed("sf")

  set.seed(1)
  dom <- create_sampling_domain(shape = "concave")
  pts <- rthomas_fast(dom, n_target = 80, mu = 8, sigma = 0.5)

  expect_s3_class(pts, "sf")
  expect_equal(nrow(pts), 80)

  inside <- as.logical(sf::st_within(pts, sf::st_union(dom), sparse = FALSE)[, 1])
  expect_true(all(inside))
})

# simulate_points_dispatch() is intentionally internal; public API is via
# rthomas_fast()/simulate_points_*_fast() and the main simulation runner.

test_that("simulate_points_strauss_fast clips to polygon and returns n_target", {
  skip_if_not_installed("sf")

  set.seed(3)
  dom <- create_sampling_domain(shape = "concave")
  pts <- simulate_points_strauss_fast(dom, n_target = 40, r = 0.6, gamma = 0.6, sweeps = 20, burnin = 0)

  expect_s3_class(pts, "sf")
  expect_equal(nrow(pts), 40)

  inside <- as.logical(sf::st_within(pts, sf::st_union(dom), sparse = FALSE)[, 1])
  expect_true(all(inside))
})
