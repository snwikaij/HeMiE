# Root estimation function

A comparatively heuristic method that is strongly dependent on the
choice of the priors. Proper priors appoint neglectable weights to
completely unrealistic values. For example, it is highly implausible
global average surface water temperatures \<0 or \>50 degrees Celsius.
At 0 it is not water but ice, and \>50 it is either a hot-spring, desert
stream, soon to evaporate or climate change has run its course and there
is no one to describe and analyse these patterns anymore.

## Usage

``` r
root_est_ppmn(
  object,
  child = NULL,
  target_value = NULL,
  prior = NULL,
  level = 0.9,
  line_col_post = "dodgerblue4",
  line_col_prior = "orange2",
  lwd = 0.9,
  plot_lab = F,
  lab_x = 0.75,
  lab_y = 0.75,
  round_lab = 1,
  lab_size = 3.5,
  kde_method = "Gaussian",
  kde_adjust = 1,
  kl_frac = 0.95,
  seed = 123
)
```

## Arguments

- object:

  A foundational or updated PPMN.

- child:

  The name of the child vertex observed in the network.

- target_value:

  The value of the child vertex observed in the network you want to
  infer the root values associated to.

- prior:

  A data frame of values describing the prior values in the root
  vertices. These need to be generated manual so the user is completely
  free what shape the prior can have and is not restricted to a
  conjugated prior.

- level:

  Credibility level (default level = 0.9).

- line_col_post:

  The colour of the generalized posterior line.

- line_col_prior:

  The colour of the prior line.

- lwd:

  Line width of the generalized posterior and prior

- plot_lab:

  Plots the point estimate (median) and standard error (se) in the
  figure. It is by default false, perhaps it is less aesthetically
  appealing (?).

- lab_x:

  Position of the label in percentage of the x-axis.

- lab_y:

  Position of the label in percentage of the y-axis.

- round_lab:

  Rounding the numbers in the label. Provide a number so it rounds the
  digits behind the comma.

- lab_size:

  Size of the symbols within the label.

- kde_method:

  The type of kernel density method (default kde_method = Gaussian) can
  be Gaussian or Gamma. The Gaussian kernal moves \<0 when this is
  technically not possible. This is not an issue for the estimated
  results, but the visual display of it. Therefore a heuristically
  Created kernel method using the Gamm distribution was created.

- kde_adjust:

  Adjustment value of the kernel bandwidth (default kde_adjust = 1; also
  see ?density).

- kl_frac:

  The strength of the update can be derived from the Kullenback-Leiber
  divergence applied over the number of simulations (default kl_frac =
  0.9) smaller values mean less precision and broader intervals.

- seed:

  Seed value 123.
