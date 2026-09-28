# Helper functions for empirical Bayesian meta analysis

Empirical Bayesian meta-analysis with random-effect (RE) or fixed-effect
(FE). For the estimation of tau^2 an Empirical Bayesian method was used
utilizing DerSimonian and Laird (1986) or DSL estimator. This approach
is easy but intervals are slightly smaller (O Bai et al. 2016). That
being said the between study variance is therefore treated as having an
improper uniform prior, which simplifies most calculations.

The heuristic methods was (naive) attempt to formulate an approximation
to a closed form solution. I started with a simple method of moment
estimator but did not progress further. The HE is still rather empirical
and results in similar results to the DSL and REML approach. \\\tau^2 =
\max\left(0, \frac{1}{n} \sum\_{i=1}^{n} (\hat{\theta_i} -
\theta\_{\text{pooled}})^2 - \frac{\sum\_{i=1}^{n} (w_i \cdot
se_i)}{\sum\_{i=1}^{n} w_i} \right)\\ However, HE will create wider
intervals over the estimate under smaller sample sizes. Yet, by default
DSL is still preferred.

## Usage

``` r
.ppmn_ebmeta(
  estimate,
  stderr,
  prior_mu = 0,
  prior_mu_se = 1000,
  prior_weights = NULL,
  a = -Inf,
  b = Inf,
  tau_2 = "DSL",
  interval = 0.9,
  RE = T,
  warnings = F
)
```

## Arguments

- estimate:

  A vector containing the effect-size

- stderr:

  A vector containing the standard error

- prior_mu:

  Prior for the mean

- prior_mu_se:

  Prior for the se

- prior_weights:

  Weights for prior odds (default=1/number of priors)

- a:

  Lower truncation bound

- b:

  Upper truncation bound

- tau_2:

  The method used to estimate tau^2 either a heuristic method 'HE0' or
  'DSL'=DerSimonian and Laird method (default = DSL).

- interval:

  Credibility intervals for the summary (default=0.9)

- RE:

  An argument indicating if RE or FE should be used (default RE=TRUE)

- warnings:

  An argument that if True returns a warning is the number of studies is
  1 (default warnings = FALSE)
