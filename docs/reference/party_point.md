# Assemblage partitioning point estimation

Assemblage partitioning point estimation

## Usage

``` r
party_point(x, y, nboot = 1000, level = 0.9, method = "weighted_x")
```

## Arguments

- x:

  A single numeric vector of the environmental gradient.

- y:

  A community matrix with presence/absence or abundance.

- nboot:

  Number of bootstraps (default = 1000).

- level:

  Credibility level (default = 0.9).

- method:

  The method deriving the intervals for the partitioning point.
