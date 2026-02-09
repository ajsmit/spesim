# Pretty-print a compact interaction matrix summary

Prints the radius and a sorted list of non-1.0 coefficients from the
interaction matrix.

## Usage

``` r
print_interactions(I, digits = 3, top_n = NULL)
```

## Arguments

- I:

  A list with `radius` and `matrix`.

- digits:

  Integer; number of digits to print.

- top_n:

  Integer; show at most this many non-1.0 entries.

## Value

Invisibly returns `NULL`.

## Examples

``` r
if (FALSE) { # \dontrun{
print_interactions(I, digits = 3, top_n = 20)
} # }
```
