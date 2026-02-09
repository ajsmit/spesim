# Validate an interaction specification

Checks that the resolved interaction list (radius + matrix) is well
formed for the given species set. Warns or errors on
shape/NA/finite/name issues.

## Usage

``` r
validate_interactions(I, spp_names, stop_on_error = FALSE)
```

## Arguments

- I:

  A list with elements `radius` (numeric scalar) and `matrix` (S×S
  numeric).

- spp_names:

  Character vector of expected species names (e.g., `LETTERS[1:S]`).

- stop_on_error:

  Logical; stop on validation failure (`TRUE`) or just warn (`FALSE`).

## Value

Invisibly returns `TRUE` on success.

## Examples

``` r
if (FALSE) { # \dontrun{
I <- list(radius = 2, matrix = diag(10))
validate_interactions(I, LETTERS[1:10])
} # }
```
