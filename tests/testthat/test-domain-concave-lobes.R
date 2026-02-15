test_that("create_sampling_domain(concave) accepts concave_lobes in 2..5", {
  skip_if_not_installed("sf")

  for (k in 2:5) {
    set.seed(1)
    dom <- spesim::create_sampling_domain(shape = "concave", concave_lobes = k, n_vertices = 80, size = 10)
    expect_true(isTRUE(as.logical(sf::st_is_valid(dom)[1])))
  }
})
