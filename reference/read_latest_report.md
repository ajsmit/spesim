# Read the most recent simulation report from an output directory

Convenience helper to locate and read the latest `*_report.txt` written
by `run_spatial_simulation(write_outputs = TRUE)` under a given
`output_dir`. Useful for users who want the report text without calling
[`generate_full_report()`](https://ajsmit.github.io/spesim/reference/generate_full_report.md)
directly.

## Usage

``` r
read_latest_report(output_dir = "out")
```

## Arguments

- output_dir:

  Directory to search. If you set `OUTPUT_PREFIX` in your init file to
  something like `"out/run"`, pass `"out"` here (default). The function
  searches recursively for files matching `_report.txt`.

## Value

A single character scalar containing the report contents.

## Examples

``` r
if (FALSE) { # \dontrun{
txt <- read_latest_report("out")
cat(txt)
} # }
```
