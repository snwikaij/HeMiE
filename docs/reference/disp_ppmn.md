# Displacement test for PPMN

Displacement test for PPMN

## Usage

``` r
disp_ppmn(
  object,
  new_data,
  nsim = 3000,
  level = 0.9,
  gd_display = F,
  seed = 123
)
```

## Arguments

- object:

  A foundational or updated PPMN

- new_data:

  New dataset with all root nodes and least one child node present

- nsim:

  Number of simulations needed for updating

- level:

  Credibility level (default level = 0.9)

- gd_display:

  Display the Generalized Displacement (GD) in the name

- seed:

  Seed value 123.
