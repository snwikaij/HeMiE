#' Helper function for binary/prop log loss
#'
#' @param x_raw Raw observations.
#' @param x_hat Predictions on these observations.
#'
#' @keywords internal
#'
.ppmn_logistic_loss <- function(x_raw, x_hat){

  #constant
  eps <- .Machine$double.eps

  x_hat <- pmin(pmax(x_hat, eps), 1-eps)
  l <- -(x_raw*log(x_hat)+(1-x_raw)*log(1-x_hat))

  #if prop
  mid <- which(is.finite(x_raw) & x_raw > 0 & x_raw < 1)
  l[mid] <- l[mid]+x_raw[mid]*log(x_raw[mid])+(1-x_raw[mid])*log(1-x_raw[mid])

  return(l)}
