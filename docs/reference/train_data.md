# Example training data for updating a PPMN

Multivariate ecological data used to demonstrate sequential updating of
a foundational Posterior Predictive Meta-analytic Network (PPMN).

## Usage

``` r
train_data
```

## Format

A data frame in which rows represent observations and columns represent
ecological variables corresponding to vertices in the example PPMN.

## Details

The dataset provides observations for multiple vertices of the example
network and is used by
[`update_ppmn()`](https://snwikaij.github.io/HeMiE/reference/update_ppmn.md)
to update the foundational PPMN using multivariate information.

## Examples

``` r
data(train_data)
head(train_data)
#> # A tibble: 6 × 10
#>      TP  Temp   Chl    pH  HCO3    DO Characeae Macrophytes CO2onlyuser   DOC
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>     <dbl>       <dbl>       <dbl> <dbl>
#> 1   6    9.30  7.05  5.94   NA   5.50         0           6       0.833 NA   
#> 2  46   16.7   7.9   8.2   146.  6.04         1          16       0.375 NA   
#> 3  86.6 14.4  11     7.84   NA   6.96         1          24       0.292  8.35
#> 4 195.   7.8  22.4   8.79   NA  12.5          0           7       0.286 NA   
#> 5  44.6 19    13     7.85   NA   8.68         1          14       0.286 NA   
#> 6  12.0 17.0   2.6   5.86   NA  NA            0           6       0.833 NA   
```
