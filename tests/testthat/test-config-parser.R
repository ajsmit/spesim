test_that("maximal init parses quoted c(...) vectors", {
  f <- system.file("examples/spesim_init_maximal.txt", package = "spesim")
  skip_if(f == "")
  P <- spesim::load_config(f)
  expect_true(is.list(P))
  expect_true("GRADIENT" %in% names(P))
})

test_that("basic init exists and loads", {
  f <- system.file("examples/spesim_init_basic.txt", package = "spesim")
  expect_true(nzchar(f))
  P <- spesim::load_config(f)
  expect_true(is.list(P))
})
