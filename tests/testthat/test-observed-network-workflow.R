test_that("spesim_from_observed builds result and graph-aware distance decay runs", {
  spa <- read.csv(system.file("extdata", "DoubsSpa.csv", package = "spesim"), check.names = FALSE)
  env <- read.csv(system.file("extdata", "DoubsEnv.csv", package = "spesim"), check.names = FALSE)
  spe <- read.csv(system.file("extdata", "DoubsSpe.csv", package = "spesim"), check.names = FALSE)

  names(spa)[1] <- "site"
  names(env)[1] <- "site"
  names(spe)[1] <- "site"
  names(spa)[names(spa) == "X"] <- "x"
  names(spa)[names(spa) == "Y"] <- "y"
  dat <- merge(spa, env, by = "site")
  dat$linear_pos <- (dat$dfs - min(dat$dfs)) / (max(dat$dfs) - min(dat$dfs))

  A <- spe
  C <- dat[, c("site", "x", "y", "linear_pos")]
  E <- dat[, c("site", "dfs", "alt", "flo", "pH", "nit", "oxy")]
  edges <- spesim::build_network_edges(C, order_by = "linear_pos", directed = TRUE)

  res <- suppressWarnings(spesim::spesim_from_observed(
    abund_matrix = A,
    site_coords = C,
    site_env = E,
    edges = edges,
    directed = TRUE
  ))

  expect_s3_class(res, "spesim_result")
  expect_true(is.matrix(attr(res$site_coords, "graph_dist_matrix")))

  dd <- suppressWarnings(spesim::calculate_distance_decay(
    res$abund_matrix,
    res$site_coords,
    metric = "along_path"
  ))
  expect_true(nrow(dd) > 0)
  expect_true(all(c("Distance", "Dissimilarity") %in% names(dd)))
})

test_that("simulate_between_observed returns interpolated sites", {
  # Minimal synthetic observed object
  A <- data.frame(
    site = c("s1", "s2", "s3"),
    A = c(3, 2, 1),
    B = c(0, 1, 2),
    C = c(1, 0, 1),
    check.names = FALSE
  )
  C <- data.frame(
    site = c("s1", "s2", "s3"),
    x = c(0, 1, 2),
    y = c(0, 0, 0),
    linear_pos = c(0, 0.5, 1)
  )
  E <- data.frame(site = c("s1", "s2", "s3"), temp = c(10, 11, 12))
  edges <- spesim::build_network_edges(C, order_by = "linear_pos", directed = TRUE)
  res <- spesim::spesim_from_observed(A, C, site_env = E, edges = edges, directed = TRUE)

  out <- spesim::simulate_between_observed(res, n_new_sites = 4, direction_bias = 0.2)
  expect_equal(nrow(out$site_coords), 4)
  expect_equal(nrow(out$abund_matrix), 4)
  expect_true(all(c("site", "x", "y", "linear_pos") %in% names(out$site_coords)))
})
