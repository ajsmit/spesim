library(spesim)

# Quick smoke test:
run_spatial_simulation(
  init_file = "_test/spesim_init.txt",
  # interactions_file = "_test/interactions_init.txt",
  output_prefix = "_test/simulation_output"
)

P <- load_config("inst/examples/spesim_init_complete.txt")

tools::showNonASCIIfile("R/report.R")




library(spesim)

## A) File-driven
run_spatial_simulation(
  # init_file = system.file("examples", "spesim_init_complete.txt", package = "spesim"),
  init_file = "_test/spesim_init.txt",
  interactions_file = NULL,
  output_prefix = "_test/ajs1"
)

## B) Programmatic (as in README)
dom1 <- create_sampling_domain()
dom2 <- create_sampling_domain()
init <- system.file("examples", "spesim_init_complete.txt", package = "spesim")
P1 <- load_config(init)
P2 <- load_config("_test/spesim_init.txt")
res <- run_spatial_simulation(
  domain = dom1,
  P = P2,
  interactions_file = NULL,
  write_outputs = TRUE,
  output_prefix = "_test/ajs4"
)
str(res$abund_matrix)

plot(dom$geometry, col = "grey95", border = "grey60")




# Step 1: Parameters & domain
P <- list(
  N_SPECIES = 6,
  N_INDIVIDUALS = 500,
  DOMINANT_FRACTION = 0.4,
  FISHER_ALPHA = 4,
  FISHER_X = 0.7
)
domain <- sf::st_as_sf(data.frame(id = 1),
  wkt = "POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))",
  crs = 4326
)

plot(domain)

# Step 2: Run simulation
run_spatial_simulation(
  P = P1,
  domain = dom1,
  interactions_file = NULL,
  write_outputs = TRUE,
  output_prefix = "_test/ajs5"
)

# Step 3: Analysis panel
panel <- generate_advanced_panel(res)
print(panel)

# Step 4: Text report
cat(generate_full_report(res)) # needs TLC







# Step 1: Parameters & domain
# using full example parameters
init <- system.file("examples", "spesim_init_complete.txt", package = "spesim")
P <- load_config(init)
domain <- create_sampling_domain()

# Step 2: Run simulation
res <- run_spatial_simulation(P, domain)

# Step 3: Analysis panel
generate_advanced_panel(res)

# Step 4: Text report
generate_full_report(res) # no effect

rpt <- generate_full_report(res)
writeLines(rpt, "spesim_report.txt") # <-- nothing...
cat(generate_advanced_panel(res)) # <-- nothing


rpt <- generate_full_report(res)

## show in console
cat(rpt)

## or write to a specific file and confirm
outfile <- file.path(getwd(), "_test/spesim_report.txt")
writeLines(rpt, outfile)
file.exists(outfile)
file.info(outfile)$size
