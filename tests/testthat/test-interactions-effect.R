test_that("local interactions influence species assignment as documented", {
  skip_if_not_installed("sf")
  skip_if_not_installed("FNN")

  # Build a simple scenario where B is attracted to A.
  P <- load_config(system.file("examples/spesim_init_complete.txt", package = "spesim"))
  P$SEED <- 101
  P$N_SPECIES <- 3
  P$N_INDIVIDUALS <- 1200
  P$GRADIENT_SPECIES <- character(0)
  P$GRADIENT_ASSIGNMENTS <- character(0)
  P$GRADIENT_OPTIMA <- NULL
  P$GRADIENT_TOLERANCE <- NULL

  # Make A spatially structured to create local neighbourhoods.
  P$SPATIAL_PROCESS_A <- "thomas"
  P$A_MEAN_OFFSPRING <- 25
  P$A_CLUSTER_SCALE <- 0.6

  P$SPATIAL_PROCESS_OTHERS <- "poisson"

  dom <- create_sampling_domain(shape = "concave")

  # Neutral interactions
  P0 <- P
  P0$INTERACTION_RADIUS <- 0
  P0$INTERACTIONS_FILE <- NULL
  P0$INTERACTIONS_EDGELIST <- NULL

  res0 <- spesim_run(P0, domain = dom, write_outputs = FALSE, seed = P0$SEED, interactions_print = FALSE)

  # Attraction: B favours A neighbours
  P1 <- P
  P1$INTERACTION_RADIUS <- 1.5
  IM <- matrix(1, 3, 3, dimnames = list(c("A","B","C"), c("A","B","C")))
  IM["B", "A"] <- 5
  P1$INTERACTION_MATRIX <- IM
  P1$INTERACTIONS_FILE <- NULL
  P1$INTERACTIONS_EDGELIST <- NULL

  res1 <- spesim_run(P1, domain = dom, write_outputs = FALSE, seed = P1$SEED, interactions_print = FALSE)

  # Compare mean distance from B points to nearest A point.
  xy0 <- sf::st_coordinates(res0$species_dist)
  xy1 <- sf::st_coordinates(res1$species_dist)

  sp0 <- res0$species_dist$species
  sp1 <- res1$species_dist$species

  A0 <- xy0[sp0 == "A", , drop = FALSE]
  B0 <- xy0[sp0 == "B", , drop = FALSE]
  A1 <- xy1[sp1 == "A", , drop = FALSE]
  B1 <- xy1[sp1 == "B", , drop = FALSE]

  skip_if(nrow(A0) < 5 || nrow(B0) < 5 || nrow(A1) < 5 || nrow(B1) < 5)

  d0 <- FNN::get.knnx(A0, B0, k = 1)$nn.dist[, 1]
  d1 <- FNN::get.knnx(A1, B1, k = 1)$nn.dist[, 1]

  expect_true(mean(d1) < mean(d0))
})
