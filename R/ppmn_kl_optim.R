#' Helper function to Kullenback-Leibler divergence optimizer
#'
#' @param lambda Lambda
#' @param loss Loss
#' @param kl_target KL divergence target
#' @param prior Prior
#'
#' @keywords internal
#'
.ppmn_kl_optim <- function(lambda, loss, kl_target, prior=NULL){

  #constant
  eps <- .Machine$double.eps

  if(is.null(prior)){prior <- rep(1/length(loss),length(loss))}
  prior <- prior/sum(prior)

  w  <- .ppmn_stable_weights(loss, lambda, prior)
  KL <- sum(w*log(pmax(w,eps)/pmax(prior,eps)))

  KL-kl_target}
