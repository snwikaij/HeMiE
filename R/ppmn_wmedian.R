#' Helper function for to calculated weighted median
#'
#' @param x values of bootstrap.
#' @param w weights appointed to x.
#'
#' @keywords internal
#'
.ppmn_wmedian <- function(x, w){

  ok <- is.finite(x) & is.finite(w) & w>0

  if(!any(ok)){return(NA_real_)}

  x <- x[ok]
  w <- w[ok]

  o <- order(x)

  x <- x[o]
  w <- w[o]

  x[which(cumsum(w)>=sum(w)/2)[1]]}
