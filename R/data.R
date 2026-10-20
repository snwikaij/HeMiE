#' Example meta-analytic data for a single GBMN edge
#'
#' Example meta-analytic dataset containing parameter estimates and
#' associated uncertainty for constructing a simple Generalized Bayesian
#' Meta-analytic Network (GBMN).
#'
#' The dataset contains estimates from fitted models describing the
#' relationship between total phosphorus (`TP`) and chlorophyll-a (`Chl`).
#' It is used in the HMiE examples to demonstrate the construction of a
#' foundational GBMN.
#'
#' @format A data frame containing meta-analytic information, including:
#' \describe{
#'   \item{Source}{Identifier for the source or fitted model.}
#'   \item{Authors}{Authors associated with the source.}
#'   \item{Comment}{Additional information about the source or model.}
#'   \item{function}{Functional form used for the edge-function.}
#'   \item{transformation}{Transformation applied to the predictor variable.}
#'   \item{edge}{Directed relationship represented as a model formula.}
#'   \item{parameter}{Model parameter associated with the estimate.}
#'   \item{estimate}{Estimated value of the model parameter.}
#'   \item{error}{Standard error associated with the parameter estimate.}
#' }
#'
#' @examples
#' data(example1)
#' head(example1)
#'
#' @keywords datasets
"example1"


#' Example meta-analytic data for a multi-edge GBMN
#'
#' Example meta-analytic dataset used to demonstrate construction of a
#' Generalized Bayesian Meta-analytic Network (GBMN) containing multiple
#' ecological relationships.
#'
#' The example includes relationships involving total nitrogen (`TN`),
#' total phosphorus (`TP`), chlorophyll-a (`Chl`), and phytoplankton
#' taxa richness (`Taxa`).
#'
#' @format A data frame containing meta-analytic parameter estimates,
#' standard errors, edge definitions, functional forms, and associated
#' information used by [build_gbmn()].
#'
#' @examples
#' data(example2)
#' head(example2)
#'
#' @keywords datasets
"example2"


#' Example meta-analytic data for a larger GBMN
#'
#' Example meta-analytic dataset used to construct a larger ecological
#' Generalized Bayesian Meta-analytic Network (GBMN).
#'
#' The corresponding example network contains several ecological variables
#' and edge-functions and is used to demonstrate construction and subsequent
#' updating of a GBMN.
#'
#' @format A data frame containing meta-analytic parameter estimates,
#' standard errors, edge definitions, functional forms, transformations,
#' and associated information used by [build_gbmn()].
#'
#' @examples
#' data(meta_data)
#' head(meta_data)
#'
#' @keywords datasets
"meta_data"


#' Example training data for updating a GBMN
#'
#' Multivariate ecological data used to demonstrate sequential updating
#' of a foundational Generalized Bayesian Meta-analytic Network (GBMN).
#'
#' The dataset provides observations for multiple vertices of the example
#' network and is used by [update_gbmn()] to update the foundational GBMN
#' using multivariate information.
#'
#' @format A data frame in which rows represent observations and columns
#' represent ecological variables corresponding to vertices in the example
#' GBMN.
#'
#' @examples
#' data(train_data)
#' head(train_data)
#'
#' @keywords datasets
"train_data"


#' Example validation data for a GBMN
#'
#' Multivariate ecological data used to demonstrate predictive assessment
#' of an updated Generalized Bayesian Meta-analytic Network (GBMN).
#'
#' The dataset is intended for evaluating predictions after the example
#' network has been updated using `train_data`.
#'
#' @format A data frame in which rows represent observations and columns
#' represent ecological variables corresponding to vertices in the example
#' GBMN.
#'
#' @examples
#' data(val_data)
#' head(val_data)
#'
#' @keywords datasets
"val_data"
