# Open spesim documentation pages in your browser

Opens pkgdown documentation in a web browser. This is a small usability
helper for teaching settings.

## Usage

``` r
spesim_open_docs(
  page = c("home", "workflow", "quadrat-placement", "public-api", "recipes"),
  url_base = "https://ajsmit.github.io/spesim/"
)
```

## Arguments

- page:

  One of: `"home"`, `"workflow"`, `"quadrat-placement"`, `"public-api"`,
  or `"recipes"`.

- url_base:

  Optional base URL. Defaults to the package URL in DESCRIPTION.
