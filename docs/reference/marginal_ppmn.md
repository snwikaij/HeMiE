# Plots the marginal expectations of each edge function with credibility intervals

This is essentially an added variable plot. It displays the marginal
relation of the dependent variable and all independent variables for
each edge-function belonging to the selected child. If multiple
edge-functions predict the same child, a separate marginal plot is
returned for every edge-function and parent combination.

## Usage

``` r
marginal_ppmn(
  object,
  child,
  data,
  level = 0.9,
  nsim = 1000,
  n_grid = 100,
  smooth_f = 0.5,
  ylim = NULL,
  ylab = NULL,
  lwd_exp = 0.8,
  lwd_col = "dodgerblue4",
  size_title = 8,
  size_text = 8,
  int_col = "dodgerblue3",
  int_alpha = 0.3,
  pt_size = 2,
  pt_alpha = 0.25,
  seed = 123
)
```

## Arguments

- object:

  A foundational or updated PPMN which will be used for predictions.

- child:

  The dependent (child) variable to display.

- data:

  Data that contains x and y coordinates (points) to plot.

- level:

  Credibility level (default level = 0.9).

- nsim:

  Number of simulations to construct the intervals.

- n_grid:

  Grid size of the added variable plot.

- smooth_f:

  A value for the smoother span only used to display specific non-linear
  models.

- ylim:

  Limit of the y-axis.

- ylab:

  Title of the y-axis.

- lwd_exp:

  Line width of the expected value.

- lwd_col:

  Colour of the expected value.

- size_title:

  Title size of the y-axis label.

- size_text:

  Text size of the title.

- int_col:

  Colour of the credibility intervals around the expected value.

- int_alpha:

  Transparency of the credibility intervals around the expected value.

- pt_size:

  Size of the points in the plot.

- pt_alpha:

  Transparency of the points in the plot.

- seed:

  Seed value 123.
