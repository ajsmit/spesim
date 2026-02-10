test_that("plot_quadrats returns a ggplot and handles empty quadrats", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("sf")

  dom <- spesim::create_sampling_domain()

  # empty sf
  empty <- sf::st_sf(quadrat_id = integer(0), geometry = sf::st_sfc(crs = sf::st_crs(dom)))

  p0 <- spesim::plot_quadrats(dom, empty, title = "Empty")
  expect_true(inherits(p0, "ggplot"))

  qs <- spesim::place_quadrats(dom, n_quadrats = 5, quadrat_size = c(1.2, 1.2))
  p1 <- spesim::plot_quadrats(dom, qs, title = "Non-empty", label_ids = "auto", label_threshold = 3)
  expect_true(inherits(p1, "ggplot"))
})

test_that("plot_quadrats can overlay Voronoi cells when present", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("sf")

  dom <- spesim::create_sampling_domain()
  set.seed(1)
  qs <- spesim::place_quadrats_voronoi(dom, 10, c(1.2, 1.2), 4, show_voronoi = TRUE)

  p <- spesim::plot_quadrats(dom, qs, show_voronoi = TRUE, show_seeds = TRUE)
  expect_true(inherits(p, "ggplot"))
})
