#' Example meta-analytic data for a single PPMN edge
#'
#' Example meta-analytic dataset containing parameter estimates and
#' associated uncertainty for constructing a simple Posterior Predictive
#' Meta-analytic Network (PPMN).
#'
#' The dataset contains estimates from fitted models describing the
#' relationship between total phosphorus (`TP`) and chlorophyll-a (`Chl`).
#' It is used in the HMiE examples to demonstrate the construction of a
#' foundational PPMN.
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
#' data(example4)
#' head(example4)
#'
#' @keywords datasets
"example4"


#' Example meta-analytic data for a multi-edge PPMN
#'
#' Example meta-analytic dataset used to demonstrate construction of a
#' Posterior Predictive Meta-analytic Network (PPMN) containing multiple
#' ecological relationships.
#'
#' The example includes relationships involving total nitrogen (`TN`),
#' total phosphorus (`TP`), chlorophyll-a (`Chl`), and phytoplankton
#' taxa richness (`Taxa`).
#'
#' @format A data frame containing meta-analytic parameter estimates,
#' standard errors, edge definitions, functional forms, and associated
#' information used by [build_ppmn()].
#'
#' @examples
#' data(example5)
#' head(example5)
#'
#' @keywords datasets
"example5"


#' Example meta-analytic data for a larger PPMN
#'
#' Example meta-analytic dataset used to construct a larger ecological
#' Posterior Predictive Meta-analytic Network (PPMN).
#'
#' The corresponding example network contains several ecological variables
#' and edge-functions and is used to demonstrate construction and subsequent
#' updating of a PPMN.
#'
#' @format A data frame containing meta-analytic parameter estimates,
#' standard errors, edge definitions, functional forms, transformations,
#' and associated information used by [build_ppmn()].
#'
#' @examples
#' data(example6)
#' head(example6)
#'
#' @keywords datasets
"example6"


#' Example training data for updating a PPMN
#'
#' Multivariate ecological data used to demonstrate sequential updating
#' of a foundational Posterior Predictive Meta-analytic Network (PPMN).
#'
#' The dataset provides observations for multiple vertices of the example
#' network and is used by [update_ppmn()] to update the foundational PPMN
#' using multivariate information.
#'
#' @format A data frame in which rows represent observations and columns
#' represent ecological variables corresponding to vertices in the example
#' PPMN.
#'
#' @examples
#' data(train_data)
#' head(train_data)
#'
#' @keywords datasets
"train_data"


#' Example validation data for a PPMN
#'
#' Multivariate ecological data used to demonstrate predictive assessment
#' of an updated Posterior Predictive Meta-analytic Network (PPMN).
#'
#' The dataset is intended for evaluating predictions after the example
#' network has been updated using `train_data`.
#'
#' @format A data frame in which rows represent observations and columns
#' represent ecological variables corresponding to vertices in the example
#' PPMN.
#'
#' @examples
#' data(val_data)
#' head(val_data)
#'
#' @keywords datasets
"val_data"
