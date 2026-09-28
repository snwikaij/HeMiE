# Predictive fit function to assess performance

Predictive fit function to assess performance

## Usage

``` r
pred_fit_ppmn(
  object,
  new_data,
  nsim = 1000,
  level = 0.9,
  rotate_x = 0,
  hjust = 0,
  seed = 123
)
```

## Arguments

- object:

  A foundational or updated PPMN.

- new_data:

  New dataset with all vertices that are roots and least one child.

- nsim:

  The number of simulations.

- level:

  Credibility level (default level = 0.9)

- rotate_x:

  The angle to rotate the variable names on x-axis.

- hjust:

  A number determining the height of the variable names on the x-axis.

- seed:

  Seed value 123.
