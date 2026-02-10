# Copy a shipped example init file into your working directory

Convenience helper for workshops: copy an example init file from the
package into the current directory (or a specified location) so students
can edit it.

## Usage

``` r
spesim_use_example(
  example = "spesim_init_basic.txt",
  to = example,
  overwrite = FALSE
)
```

## Arguments

- example:

  Basename of an example file (e.g. `"spesim_init_basic.txt"`).

- to:

  Destination path. Default is `example` in the current working
  directory.

- overwrite:

  Logical; overwrite an existing file?

## Value

The destination path (invisibly).
