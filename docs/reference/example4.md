# Example meta-analytic data for a single PPMN edge

Example meta-analytic dataset containing parameter estimates and
associated uncertainty for constructing a simple Posterior Predictive
Meta-analytic Network (PPMN).

## Usage

``` r
example4
```

## Format

A data frame containing meta-analytic information, including:

- Source:

  Identifier for the source or fitted model.

- Authors:

  Authors associated with the source.

- Comment:

  Additional information about the source or model.

- function:

  Functional form used for the edge-function.

- transformation:

  Transformation applied to the predictor variable.

- edge:

  Directed relationship represented as a model formula.

- parameter:

  Model parameter associated with the estimate.

- estimate:

  Estimated value of the model parameter.

- error:

  Standard error associated with the parameter estimate.

## Details

The dataset contains estimates from fitted models describing the
relationship between total phosphorus (`TP`) and chlorophyll-a (`Chl`).
It is used in the HMiE examples to demonstrate the construction of a
foundational PPMN.

## Examples

``` r
data(example4)
head(example4)
#> # A tibble: 6 × 9
#>   Source Authors Comment `function` transformation edge   parameter estimate
#>   <chr>  <chr>   <lgl>   <chr>      <chr>          <chr>  <chr>        <dbl>
#> 1 A      A       NA      log        log            Chl~TP b0           -0.8 
#> 2 A      A       NA      log        log            Chl~TP b1            1.1 
#> 3 B      A       NA      log        log            Chl~TP b0           -0.2 
#> 4 B      A       NA      log        log            Chl~TP b1            0.93
#> 5 C      C       NA      log        log            Chl~TP b0           -1.1 
#> 6 C      C       NA      log        log            Chl~TP b1            1.05
#> # ℹ 1 more variable: error <dbl>
```
