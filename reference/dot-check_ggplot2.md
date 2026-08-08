# Check ggplot2 availability

Verifies that ggplot2 is installed. Called at the top of every
visualization function since ggplot2 is in Suggests, not Imports.

## Usage

``` r
.check_ggplot2()
```

## Value

Invisible NULL. Throws an error if ggplot2 is not installed.
