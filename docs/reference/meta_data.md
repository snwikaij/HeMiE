# Example meta-analytic data for a larger GBMN

Example meta-analytic dataset used to construct a larger ecological
Generalized Bayesian Meta-analytic Network (GBMN).

## Usage

``` r
meta_data
```

## Format

A data frame containing meta-analytic parameter estimates, standard
errors, edge definitions, functional forms, transformations, and
associated information used by
[`build_gbmn()`](https://snwikaij.github.io/HeMiE/reference/build_gbmn.md).

## Details

The corresponding example network contains several ecological variables
and edge-functions and is used to demonstrate construction and
subsequent updating of a GBMN.

## Examples

``` r
data(meta_data)
head(meta_data)
#> # A tibble: 6 × 10
#>   Source         Authors Type  Comment `function` transformation edge  parameter
#>   <chr>          <chr>   <chr> <chr>   <chr>      <chr>          <chr> <chr>    
#> 1 Macrophyte gr… Kaijse… Raw … NA      class      NA             Char… b0       
#> 2 Macrophyte gr… Kaijse… Raw … NA      class      NA             Char… b1       
#> 3 Macrophyte gr… Kaijse… Raw … NA      class      NA             Char… b2       
#> 4 Macrophyte gr… Kaijse… Raw … NA      class      NA             Char… b3       
#> 5 Macrophyte gr… Kaijse… Raw … NA      class      NA             Char… b4       
#> 6 Charophyte va… Kolada… Supp… NA      class      NA             Char… b0       
#> # ℹ 2 more variables: estimate <dbl>, error <dbl>
```
