# Check for (optional) spatstat availability

These simulations no longer require spatstat. This helper remains as a
no-op stub to avoid breaking any internal calls that previously relied
on it. If you later re-enable spatstat-powered diagnostics, you can
repurpose this function to warn or fail accordingly.

## Usage

``` r
.require_spatstat()
```
