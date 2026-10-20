# Example meta-analytic data for a multi-edge GBMN

Example meta-analytic dataset used to demonstrate construction of a
Generalized Bayesian Meta-analytic Network (GBMN) containing multiple
ecological relationships.

## Usage

``` r
example2
```

## Format

A data frame containing meta-analytic parameter estimates, standard
errors, edge definitions, functional forms, and associated information
used by
[`build_gbmn()`](https://snwikaij.github.io/HeMiE/reference/build_gbmn.md).

## Details

The example includes relationships involving total nitrogen (`TN`),
total phosphorus (`TP`), chlorophyll-a (`Chl`), and phytoplankton taxa
richness (`Taxa`).

## Examples

``` r
data(example2)
head(example2)
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
