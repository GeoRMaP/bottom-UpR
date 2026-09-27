# bottom-UpR

`bottom-UpR` is a generic R package for Bayesian bottom-up population estimation
from sparse enumeration data, mapped buildings and geospatial covariates.

## Core idea

The user can supply only the response column and let the package infer its
support:

- non-negative integer response → **COUNT**
- strictly positive continuous response → **PPB**

It then automatically compares the compatible likelihoods:

| Response | Candidate likelihoods |
|---|---|
| COUNT | Poisson, Negative Binomial |
| PPB | Gamma, Lognormal |

Automatic likelihood selection uses INLA conditional predictive ordinates (CPO):
the candidate with the highest mean leave-one-out log predictive density is
selected. Lognormal PPB is implemented as Gaussian modelling of `log(PPB)` and
the predictive density is Jacobian-adjusted back to the PPB scale before model
comparison.

Because a rounded PPB variable can look integer-valued, automatic response
classification is necessarily heuristic. Users can always set
`response_type = "ppb"` or `"count"` explicitly.

## Installation

INLA is distributed from the R-INLA repository rather than ordinary CRAN:

```r
install.packages(
  "INLA",
  repos = c(
    getOption("repos"),
    INLA = "https://inla.r-inla-download.org/R/stable"
  ),
  dep = TRUE
)
```

Then install from GitHub:

```r
install.packages("remotes")
remotes::install_github("GeoRMaP/bottom-UpR")
```

## Minimal examples

```r
fit_count <- bottom_up(
  data = ea,
  response = "population",
  buildings = "buildings",
  covariates = c("roads", "water", "market"),
  coords = c("x", "y"),
  random_effects = c("settlement", "source")
)

fit_ppb <- bottom_up(
  data = ea,
  response = "ppb",
  buildings = "buildings",
  covariates = c("roads", "water", "market"),
  coords = c("x", "y")
)

fit_count
summary(fit_count)
```

## Explicit override

```r
fit <- bottom_up(
  ea,
  response = "population",
  buildings = "buildings",
  response_type = "count",
  likelihood = "negative binomial"
)
```

## Design principle

Automatic likelihood selection is based on predictive evidence, not only a
variance-to-mean rule. The package first restricts candidates to likelihoods
whose support matches the response and then compares their leave-one-out
predictive densities.

For publication analyses, report the selected likelihood, all candidate scores,
the response-type rule/override, spatial priors, mesh, covariates and validation
design.

## Development status

Version 0.1.0 is an early development release and has not yet been validated with a full `R CMD check` in an environment with R-INLA installed.
