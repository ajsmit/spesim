test_that("generate_sad returns correct length, names, and sum", {
  set.seed(1)
  S <- 10
  N <- 500

  models <- c(
    "fisher",
    "geometric",
    "brokenstick",
    "zipf",
    "zipf-mandelbrot",
    "lognormal",
    "poisson-lognormal",
    "poisson-gamma",
    "custom"
  )

  for (m in models) {
    if (m == "custom") {
      sad <- rep(1 / S, S)
      x <- generate_sad(S, N, model = m, sad = sad)
    } else if (m == "fisher") {
      x <- generate_sad(S, N, model = m, dominant_fraction = 0.3, alpha = 3, x = 0.95)
    } else if (m == "geometric") {
      x <- generate_sad(S, N, model = m, k = 0.4)
    } else if (m == "zipf") {
      x <- generate_sad(S, N, model = m, exponent = 1.2)
    } else if (m == "zipf-mandelbrot") {
      x <- generate_sad(S, N, model = m, exponent = 1.2, q = 2)
    } else if (m %in% c("lognormal", "poisson-lognormal")) {
      x <- generate_sad(S, N, model = m, meanlog = 0, sdlog = 1)
    } else if (m == "poisson-gamma") {
      x <- generate_sad(S, N, model = m, shape = 1.5, rate = 1)
    } else {
      x <- generate_sad(S, N, model = m)
    }

    expect_type(x, "integer")
    expect_length(x, S)
    expect_identical(names(x), LETTERS[1:S])
    expect_true(all(x >= 0))
    expect_equal(sum(x), N)
  }
})

test_that("custom generator function works", {
  set.seed(1)
  S <- 8
  N <- 200

  f <- function(n_species, n_individuals, ...) {
    # return unnormalised weights
    seq_len(n_species)
  }

  x <- generate_sad(S, N, model = "custom", sad = f)
  expect_length(x, S)
  expect_equal(sum(x), N)
})
