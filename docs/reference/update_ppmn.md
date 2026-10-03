# Update PPMN function

The updating function uses loss-based Generalized Bayesian Updating
(GBU) to update the prior to a generalized posterior (often called GIBBS
posterior). In simplicity it simulates parameters (particles) from prior
model These particles are used in a function to predict. Each prediction
the mismatch with the data (loss) is calculated. To updated the prior
these losses each particles is then weighted inversely relative to the
losses.

## Usage

``` r
update_ppmn(
  object,
  new_data,
  nsim = 3000,
  level = 0.9,
  method = "kernel",
  up_strength = 0.95,
  max_lambda = 1000,
  diagnostics = T,
  covar = T,
  seed = 123
)
```

## Arguments

- object:

  A foundational or updated PPMN.

- new_data:

  New data set with all root vertices and least one child vertex
  present.

- nsim:

  Number of simulations used for updating.

- level:

  Generalized credibility level (default level = 0.9).

- method:

  The method being used to set the updating strength of the model. The
  default method is 'kernel', The other method 'manual' uses the
  Kullenback-Leibler divergence applied over the number of simulations
  to express the concentration strength relative to no change relative
  to the prior. Here the default of the update strength is 95%.

- up_strength:

  The fraction of information that is allowed to concentrate in the
  generalized posterior when using the method 'manual'.

- max_lambda:

  Maximum lambda that can be set for Generalized Bayesian Updating
  (default max_lambda = 1000).

- diagnostics:

  If residual diagnostics should be returned (default diagnostics = T).

- covar:

  Covariance matrix that is calculated between model parameters (default
  covar = T).

- seed:

  Seed value 123.
