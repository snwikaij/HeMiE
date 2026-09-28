# Example meta-analytic data for a multi-edge PPMN

Example meta-analytic dataset used to demonstrate construction of a
Posterior Predictive Meta-analytic Network (PPMN) containing multiple
ecological relationships.

## Usage

``` r
example5
```

## Format

A data frame containing meta-analytic parameter estimates, standard
errors, edge definitions, functional forms, and associated information
used by
[`build_ppmn()`](https://snwikaij.github.io/HeMiE/reference/build_ppmn.md).

## Details

The example includes relationships involving total nitrogen (`TN`),
total phosphorus (`TP`), chlorophyll-a (`Chl`), and phytoplankton taxa
richness (`Taxa`).

## Examples

``` r
data(example5)
head(example5)
#> # A tibble: 6 × 9
#>   Source Authors Comment `function` transformation edge      parameter estimate
#>   <chr>  <chr>   <lgl>   <chr>      <chr>          <chr>     <chr>        <dbl>
#> 1 A      A       NA      log        log            Chl~TN+TP b0           -0.8 
#> 2 A      A       NA      log        log            Chl~TN+TP b1            0.34
#> 3 A      A       NA      log        log            Chl~TN+TP b2            1.1 
#> 4 B      B       NA      log        log            Chl~TN+TP b0           -0.2 
#> 5 B      B       NA      log        log            Chl~TN+TP b1            0.7 
#> 6 B      B       NA      log        log            Chl~TN+TP b2            0.93
#> # ℹ 1 more variable: error <dbl>
```
