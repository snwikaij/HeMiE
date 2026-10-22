# Bootstrap coverage of generalized posterior expected predictions

Very experimental it concentrates the generalized posterior too fast. It
checks whether the average of all expectations is covered x% by the
credibility under bootstrapping.

## Usage

``` r
.gbmn_coverage_fast(
  loss_full,
  preds,
  obs_nodes,
  obs_mat,
  max_lambda,
  level,
  nboot,
  nsim,
  shrink = 3
)
```

## Arguments

- loss_full:

  Full loss array (N x J x S).

- preds:

  Predictions

- obs_nodes:

  Observed children

- obs_mat:

  Observed new data

- max_lambda:

  Maximum value for lambda.

- level:

  Credibility level.

- nboot:

  Number of bootstraps

- nsim:

  Number of simulations

- shrink:

  number of contractions around lambda
