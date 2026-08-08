# Compute spherical spline interpolation weights

Implements the spherical spline method of Perrin et al. (1989) for
interpolating EEG channel data. Uses Legendre polynomials to compute the
interpolation matrix G and the weight vector g_x for each bad channel.

## Usage

``` r
.spherical_spline_weights(good_pos, bad_pos, order = 4)
```

## Arguments

- good_pos:

  Matrix with columns (pos_x, pos_y, pos_z) for good channels.

- bad_pos:

  Matrix with columns (pos_x, pos_y, pos_z) for bad channels.

- order:

  Maximum Legendre polynomial order (default: 4).

## Value

Matrix of weights (n_bad x n_good) for interpolation.
