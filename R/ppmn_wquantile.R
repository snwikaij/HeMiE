#' Helper function for to calculated weighted quantiles.
#'
#' @param x Values (particles/simulations).
#' @param w GBU weights.
#' @param probs Probabilities.
#'
#' @keywords internal
#'
.ppmn_wquantile <- function(x, w, probs){

  ok <- is.finite(x) & is.finite(w)
  x  <- x[ok]
  w  <- w[ok]

  if(length(x)==0){return(rep(NA,length(probs)))}

  o  <- order(x)
  x  <- x[o]
  w  <- w[o]
  w  <- w/sum(w)
  cw <- cumsum(w)

  sapply(probs, function(p) x[which(cw>=p)[1]])}
