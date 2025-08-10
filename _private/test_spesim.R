library(spesim)

# Quick smoke test:
run_spatial_simulation(
  init_file = system.file("examples/simul_init.txt", package = "spesim"),
  interactions_file = system.file("examples/interactions_init.txt", package = "spesim"),
  output_prefix = file.path(tempdir(), "simulation_output")
)



P <- load_config("inst/examples/spesim_init_complete.txt")
