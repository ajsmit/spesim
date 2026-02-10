test_that("place_quadrats_voronoi(show_voronoi=TRUE) returns attributes", {
  skip_if_not_installed("sf")

  dom <- spesim::create_sampling_domain()
  set.seed(1)
  qs <- spesim::place_quadrats_voronoi(
    dom,
    n_quadrats = 10,
    quadrat_size = c(1.2, 1.2),
    voronoi_seed_factor = 4,
    show_voronoi = TRUE
  )

  expect_s3_class(qs, "sf")
  expect_true(nrow(qs) > 0)

  vc <- attr(qs, "voronoi_cells")
  vs <- attr(qs, "voronoi_seeds")
  ic <- attr(qs, "inscribed_circles")

  expect_true(!is.null(vc))
  expect_true(!is.null(vs))
  expect_true(!is.null(ic))

  expect_s3_class(vc, "sf")
  # seeds are an sfc
  expect_true(inherits(vs, "sfc"))
  expect_true(inherits(ic, "sfc"))
})
