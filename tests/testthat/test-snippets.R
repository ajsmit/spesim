test_that("All vignette code chunks are syntactically valid", {
  skip_if_not_installed("knitr")

  # In a development checkout, vignette sources are in ./vignettes.
  # Under R CMD check, the working directory can vary, so we search a few
  # plausible locations and skip if vignette sources aren't present.
  cand_dirs <- c(
    "vignettes",
    file.path("..", "vignettes"),
    file.path("..", "..", "vignettes")
  )
  vign_dir <- cand_dirs[dir.exists(cand_dirs)][1]
  if (is.na(vign_dir) || !dir.exists(vign_dir)) {
    skip("Vignette sources not available in this test context")
  }

  rmd <- list.files(vign_dir, pattern = "\\.Rmd$", full.names = TRUE)
  if (length(rmd) == 0) {
    skip("No vignette .Rmd files found")
  }

  for (f in rmd) {
    tmp <- tempfile(fileext = ".R")

    # purl extracts *all* code chunks, including eval=FALSE examples.
    knitr::purl(f, output = tmp, documentation = 0L, quiet = TRUE)

    # Ensure the extracted code is parseable.
    expect_error(parse(tmp), NA, info = paste("Failed to parse purl from", basename(f)))
  }
})
