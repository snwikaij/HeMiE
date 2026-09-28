# Example meta-analytic data for a larger PPMN

Example meta-analytic dataset used to construct a larger ecological
Posterior Predictive Meta-analytic Network (PPMN).

## Usage

``` r
example6
```

## Format

A data frame containing meta-analytic parameter estimates, standard
errors, edge definitions, functional forms, transformations, and
associated information used by
[`build_ppmn()`](https://snwikaij.github.io/HeMiE/reference/build_ppmn.md).

## Details

The corresponding example network contains several ecological variables
and edge-functions and is used to demonstrate construction and
subsequent updating of a PPMN.

## Examples

``` r
data(example6)
head(example6)
#> # A tibble: 6 × 9
#>   Source      Authors Comment `function` transformation edge  parameter estimate
#>   <chr>       <chr>   <chr>   <chr>      <chr>          <chr> <chr>        <dbl>
#> 1 Charophyte… Kolada… NA      class      NA             Char… b0           9.5  
#> 2 Charophyte… Kolada… NA      class      NA             Char… b1          13.4  
#> 3 Charophyte… Kolada… NA      class      NA             Char… b2          37.7  
#> 4 Charophyte… Kolada… NA      class      NA             Char… b3          38.7  
#> 5 Charophyte… Kolada… NA      class      NA             Char… b4           0.530
#> 6 Macrophyte… Kaijse… NA      class      NA             Char… b0          14.3  
#> # ℹ 1 more variable: error <dbl>
```
