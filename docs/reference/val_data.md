# Example validation data for a PPMN

Multivariate ecological data used to demonstrate predictive assessment
of an updated Posterior Predictive Meta-analytic Network (PPMN).

## Usage

``` r
val_data
```

## Format

A data frame in which rows represent observations and columns represent
ecological variables corresponding to vertices in the example PPMN.

## Details

The dataset is intended for evaluating predictions after the example
network has been updated using `train_data`.

## Examples

``` r
data(val_data)
head(val_data)
#> # A tibble: 6 × 10
#>      TP  Temp   Chl    pH  HCO3    DO Characeae Macrophytes CO2onlyuser   DOC
#>   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>     <dbl>       <dbl>       <dbl> <dbl>
#> 1  860.  6.3  14.9   9.23   NA   12           0           5       0        NA
#> 2  142.  7.9  48.7   8.66   NA    8.1         1           6       0        NA
#> 3   78   7.4  22.1   8.94   NA   11.3         1          10       0.1      NA
#> 4  158.  6.57 73.3   8.75  116.  16.7         0           8       0.125    NA
#> 5  131.  5.6   6.29  8.3    NA   15.7         1           6       0        NA
#> 6  212. 13.6  19.0   8.5    NA   14           0           5       0.2      NA
```
