# Helper function for multivariate Bayesian bootstrap KL calibration

Helper function for multivariate Bayesian bootstrap KL calibration

## Usage

``` r
.ppmn_kl_bboot(
  resid1,
  resid0,
  method = "multikernel",
  n_grid = 5000,
  n_boot = 1000,
  shrink = 0.1
)
```

## Arguments

- resid1:

  residual matrix for best particle.

- resid0:

  residual array from prior.

- n_grid:

  size of the grid in grid approximation.

- n_boot:

  number of Bayesian bootstrap replicates.

- shrink:

  covariance shrinkage toward a diagonal matrix.
