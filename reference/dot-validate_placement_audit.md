# Validate a placement_audit attribute

Quadrat placement functions in spesim attach a `placement_audit`
attribute to the returned `sf` object. This is used by
[`audit_sampling_scheme()`](https://ajsmit.github.io/spesim/reference/audit_sampling_scheme.md)
and the report's conceptual audit section.

This validator is intentionally permissive: it checks for a minimal
common schema and normalises field types where possible, but it does not
error on scheme-specific fields.

## Usage

``` r
.validate_placement_audit(x, stop_on_error = FALSE)
```

## Arguments

- x:

  An object with (optional) `attr(x, "placement_audit")`.

- stop_on_error:

  Logical; stop on validation errors? Default `FALSE`.

## Value

A normalised `placement_audit` list (possibly empty).
