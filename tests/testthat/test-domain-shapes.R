test_that("create_sampling_domain() supports multiple shape families and returns valid polygons", {
  skip_if_not_installed("sf")

  shapes <- c("concave", "convex", "ellipse", "rectangle")

  for (sh in shapes) {
    set.seed(123)
    if (sh == "concave") {
      dom <- spesim::create_sampling_domain(shape = sh, n_vertices = 60, size = 10, aspect = 1.3, concave_lobes = 3)
    } else {
      dom <- spesim::create_sampling_domain(shape = sh, n_vertices = 60, size = 10, aspect = 1.3)
    }

    expect_s3_class(dom, "sf")
    expect_equal(nrow(dom), 1L)

    gt <- as.character(sf::st_geometry_type(dom$geometry))
    expect_true(gt %in% c("POLYGON", "MULTIPOLYGON"))

    expect_true(isTRUE(as.logical(sf::st_is_valid(dom)[1])))
  }
})
