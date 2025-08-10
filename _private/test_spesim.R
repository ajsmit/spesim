library(spesim)

# Quick smoke test:
run_spatial_simulation(
  init_file = "_test/spesim_init.txt",
  # interactions_file = "_test/interactions_init.txt",
  output_prefix = "_test/simulation_output"
)


P <- load_config("inst/examples/spesim_init_complete.txt")


tools::showNonASCIIfile("R/report.R")
