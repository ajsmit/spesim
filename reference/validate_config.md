# Validate a spesim configuration (without running the simulation)

Checks a configuration file or parameter list for common problems and
returns a structured report (errors + warnings) without stopping.

This is aimed at teaching workflows where students iterate on init files
and need fast feedback before running a full simulation.

## Usage

``` r
validate_config(config, interactions_file = NULL, strict = FALSE)
```

## Arguments

- config:

  Either a path to an init file (character scalar) or an in-memory
  parameter list `P`.

- interactions_file:

  Optional interactions file path; if provided, it is validated as well.

- strict:

  Logical; if `TRUE`, treat warnings as errors.

## Value

A list with components `ok` (logical), `errors` (character), `warnings`
(character), and (when possible) `P` (resolved config).
