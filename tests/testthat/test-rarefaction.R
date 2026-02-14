test_that("calculate_rarefaction() has no graphics side-effects", {
  skip_if_not_installed("vegan")

  # Ensure we start from the null device. If calculate_rarefaction() ever
  # triggers base plotting (as vegan::rarecurve() does by default), R will try
  # to open a graphics device.
  while (grDevices::dev.cur() > 1) grDevices::dev.off()

  old <- options(device = function(...) stop("Graphics device opened"))
  on.exit(options(old), add = TRUE)

  abund <- data.frame(
    site = paste0("Q", 1:3),
    A = c(5, 2, 0),
    B = c(1, 0, 3),
    C = c(0, 1, 1)
  )

  expect_silent({
    rr <- spesim::calculate_rarefaction(abund)
  })

  expect_true(is.data.frame(rr))
  expect_named(rr, c("SiteID", "SampleSize", "RarefiedRichness"))
  expect_equal(unname(grDevices::dev.cur()), 1)
})
