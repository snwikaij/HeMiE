# Helper function for Adapted Multiple Importance Sampling

This is based on the article "Adaptive Multiple Importance Sampling" by
Cornuet et al. (2012). I keep lambda fixed at 1. and the method is still
experimental but looks promising. CORNUET, J. M., MARIN, J. M., Mira,
A., & Robert, C. P. (2012). Adaptive multiple importance sampling.
Scandinavian Journal of Statistics, 39(4), 798-812.

## Usage

``` r
.gbmn_amis_easy(
  object,
  preds,
  loss,
  new_data,
  first_mad,
  obs_nodes,
  obs_mat,
  niter,
  lambda,
  nsim
)
```

## Arguments

- object:

  An foundational or updated gbmn object.

- preds:

  Predictions of that object.

- loss:

  Loss based on those predictions.

- new_data:

  The new data used for that loss.

- first_mad:

  The Maximum Absolute Deviation of that prediction

- obs_nodes:

  The observed vertices (still inconsistent with terminology)

- niter:

  Number of adaptive iterations (T).

- lambda:

  Lambda which is by default 1.

- nsim:

  Number of simulations per iteration.

- obs_mats:

  Observed matrix after cleaning new data.
