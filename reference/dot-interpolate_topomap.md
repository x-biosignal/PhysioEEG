# Inverse distance weighted interpolation for topographic maps

Interpolates channel values onto a regular grid within a unit circle
using inverse distance weighting with power 2. Points outside the head
circle are set to NA.

## Usage

``` r
.interpolate_topomap(values, x_pos, y_pos, resolution = 100)
```

## Arguments

- values:

  Numeric vector of channel values.

- x_pos:

  Numeric vector of channel x positions.

- y_pos:

  Numeric vector of channel y positions.

- resolution:

  Integer grid resolution (default: 100).

## Value

A data.frame with columns x, y, value suitable for ggplot2::geom_raster.
