# Load interspecific interaction settings (radius + matrix)

Reads a separate interactions config file and returns a
local-interaction radius along with an \\S \times S\\ interaction matrix
for \\S\\ species (rows = focal, columns = neighbour). If the file is
missing or incomplete, a neutral matrix of ones and radius 0 are
returned.

## Usage

``` r
load_interactions(interactions_file, n_species)
```

## Arguments

- interactions_file:

  Character scalar or `NULL`. Path to the interactions configuration
  file. If `NULL` or not found, the returned radius is 0 and all
  coefficients are `1`.

- n_species:

  Integer. Number of species; defines the size of the returned matrix
  and the expected label set `LETTERS[1:n_species]`.

## Value

A list with components:

- `radius`:

  Numeric length-1; the interaction radius (0 disables interactions).

- `matrix`:

  Numeric matrix of size \\S \times S\\ with dimnames `LETTERS[1:S]`.
  Non-finite values are rejected.

## Details

The interactions file is parsed as simple `KEY = value` pairs. The
following keys are recognized:

- `INTERACTION_RADIUS`: non-negative numeric; distance threshold within
  which neighbours influence assignment probabilities.

- `MATRIX_CSV`: path to a CSV containing a (sub)matrix with row names
  and column names corresponding to species labels (e.g. A..Z). Any
  missing rows/columns are filled with `1`.

- `EDGELIST_CSV`: path to a CSV with columns `focal, neighbor, value`;
  entries are inserted into the appropriate cells of the full matrix.

- `AUTO`: optional flag (`TRUE`/`FALSE`) enabling a simple built-in
  example structure if no CSV is provided.

If both `MATRIX_CSV` and `EDGELIST_CSV` are supplied, `MATRIX_CSV` takes
precedence.

## See also

[`load_config()`](https://ajsmit.github.io/spesim/reference/load_config.md)

## Examples

``` r
if (FALSE) { # \dontrun{
I <- load_interactions("config/interactions.txt", n_species = 10)
I$radius
I$matrix[1:3, 1:3]
} # }
```
