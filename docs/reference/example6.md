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
#> # A tibble: 6 × 10
#>   Source         Authors Type  Comment `function` transformation edge  parameter
#>   <chr>          <chr>   <chr> <chr>   <chr>      <chr>          <chr> <chr>    
#> 1 Contribution … Algest… Tabl… NA      log        log            DOC~… b0       
#> 2 Contribution … Algest… Tabl… NA      log        log            DOC~… b1       
#> 3 Short-term su… Arthau… Figu… NA      log        log            Macr… b0       
#> 4 Short-term su… Arthau… Figu… NA      log        log            Macr… b1       
#> 5 Controls ofd1… Bade e… Tabl… HCO3 m… sigmoidal  NA             HCO3… b0       
#> 6 Controls ofd1… Bade e… Tabl… HCO3 m… sigmoidal  NA             HCO3… b1       
#> # ℹ 2 more variables: estimate <dbl>, error <dbl>
```
