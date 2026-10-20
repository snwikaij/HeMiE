#' Helper function to derive GBU weights
#'
#' @param loss The loss.
#' @param lambda Lambda.
#' @param prior Prior.
#'
#' @keywords internal
#'
.gbmn_stable_weights <- function(loss, lambda, prior=NULL){

  if(is.infinite(lambda) && lambda>0){

    w <- rep(0, length(loss))
    w[which.min(loss)] <- 1

    return(w)}

  #constant
  eps <- .Machine$double.eps

  if(is.null(prior)){prior <- rep(1/length(loss),length(loss))}
  prior <- prior/sum(prior)

  #weighting
  a <- log(pmax(prior, eps))-lambda*loss
  a <- a-max(a[is.finite(a)])
  w <- exp(a)
  w[!is.finite(w)] <- 0
  w/sum(w)}
