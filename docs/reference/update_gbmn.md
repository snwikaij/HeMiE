# Update GBMN function

The updating function uses loss-based Generalized Bayesian Updating
(GBU) to update the prior to a generalized posterior (often called GIBBS
posterior). In simplicity it simulates parameters (particles) from prior
model These particles are used in a function to predict. Each prediction
the mismatch with the data (loss) is calculated. To updated the prior
these losses each particles is then weighted inversely relative to the
losses.

## Usage

``` r
update_gbmn(
  object,
  new_data,
  nsim = 3000,
  level = 0.9,
  method = "hold-amis",
  up_character = NULL,
  max_lambda = 200,
  diagnostics = T,
  covar = T,
  seed = 123
)
```

## Arguments

- object:

  A foundational or updated GBMN.

- new_data:

  New data set with all root vertices and least one child vertex
  present.

- nsim:

  Number of simulations used for updating.

- level:

  Generalized credibility level (default level = 0.9).

- method:

  The method being used to set the updating strength of the model
  (default = = 'hold-amis'). The method 'manual' concentrates by
  selecting the best x% of the simulations. The method coverage
  calculates the coverage of the expected value using bootstrapped
  subset. The method 'holdout' is a version that looks similar to
  SafeBayes. The method 'amis' uses Adaptive Multiple Importance
  Sampling. The method 'hold-amis' uses a fixed 5 repetitions for the
  'holdout' method and uses the derived learning rate as a fixed value
  for amis.

- up_character:

  A parameter controlling the updating method. For method "hold-amis",
  it specifies the number adaptive iterations (default = 5). For method
  "amis", it specifies the number adaptive iterations (default = 5). For
  method "holdout", it specifies the number of repetitions (default =
  20). For method "manual", it specifies the fraction of best-performing
  simulations (default = NULL). For method "coverage", it specifies the
  number of bootstraps (default = 50).

- max_lambda:

  Maximum lambda that can be set for Generalized Bayesian Updating
  (default max_lambda = 200).

- diagnostics:

  If residual diagnostics should be returned (default diagnostics = T).

- covar:

  Covariance matrix that is calculated between model parameters (default
  covar = T).

- seed:

  Seed value 123.
