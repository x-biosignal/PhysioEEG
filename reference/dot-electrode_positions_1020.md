# Standard 10-20 electrode positions

Returns a data.frame with label and (pos_x, pos_y, pos_z) for the 19
standard 10-20 system electrodes. Coordinates are based on the standard
spherical projection to 2D (azimuthal equidistant) with z computed to
lie on the unit sphere.

## Usage

``` r
.electrode_positions_1020()
```

## Value

data.frame with columns: label, pos_x, pos_y, pos_z.
