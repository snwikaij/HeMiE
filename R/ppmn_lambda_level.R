#' Helper function finding fraction of loss
#'
#' @param loss A vector of losses.
#' @param level The level of loss.
#' @param max_lambda Maximum lambda.
#'
#' @importFrom stats uniroot
#'
#' @keywords internal
#'
.ppmn_lambda_level <- function(loss, best_frac, max_lambda){

  #extract loss
  loss <- as.numeric(loss)
  S    <- length(loss)

  #number of best particles is best_frac
  n_best_part <- max(1, ceiling(best_frac*S))
  best        <- order(loss)[seq_len(n_best_part)]

  #stable weight fun
  flamb <- function(lambda){
    w <- .ppmn_stable_weights(loss, lambda)
    sum(w[best])-(1-best_frac)}

  #set upper limit
  upper <- max_lambda

  while(flamb(upper)<0 &&upper<1e6){upper <- upper*2}

  if(flamb(upper)<0){return(upper)}

  uniroot(flamb, interval=c(0,upper))$root}
