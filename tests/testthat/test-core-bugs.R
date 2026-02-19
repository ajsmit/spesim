test_that("generate_fisher_log_series sums correctly when n_species == 1", {
  abund <- generate_fisher_log_series(
    n_species = 1,
    n_individuals = 100,
    dominant_fraction = 0.1,
    alpha = 3,
    x = 0.95
  )
  expect_equal(sum(abund), 100)
  expect_equal(names(abund), "A")
})

test_that("calculate_distance_decay aligns site order by site id", {
  abund_matrix <- data.frame(
    site = c(2, 1, 3),
    A = c(1, 0, 1),
    B = c(0, 1, 1)
  )
  site_coords <- data.frame(
    site = c(1, 2, 3),
    x = c(0, 10, 20),
    y = c(0, 0, 0)
  )

  # If alignment works, distance for pair (site1, site2) is 10 and dissimilarity finite.
  dd <- calculate_distance_decay(abund_matrix, site_coords, metric = "euclidean")
  expect_true(all(is.finite(dd$Distance)))
  expect_true(all(is.finite(dd$Dissimilarity)))
  expect_true(any(abs(dd$Distance - 10) < 1e-8))
})
