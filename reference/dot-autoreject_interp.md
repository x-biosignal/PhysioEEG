# Interpolate given channels of one epoch slice via spherical spline (internal)

Repairs the `bad_idx` channels of a 2D epoch slice (time x channels)
using the same spherical-spline weights as
[`eegInterpolate()`](https://x-biosignal.github.io/PhysioEEG/reference/eegInterpolate.md),
caching each good/bad configuration's weight matrix.

## Usage

``` r
.autoreject_interp(slice, bad_idx, pos, cache = NULL)
```

## Arguments

- slice:

  A time x channels matrix.

- bad_idx:

  Integer indices of channels to interpolate.

- pos:

  Electrode positions matrix (channels x 3).

- cache:

  Optional environment used to memoise weight matrices.

## Value

The slice with `bad_idx` channels replaced by interpolation.
