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
.ppmn_lambda_level <- function(loss, ess_frac, max_lambda){

  #get loss and nsim
  loss <- as.numeric(loss)
  S <- length(loss)

  #t smal or full loss matrix is the whole index
  if(S<2){return(0)}
  if(all(loss == loss[1])){return(0)}

  #if ess_frac
  if(ess_frac == 0){return(Inf)}
  if(ess_frac == 1){return(0)}
  target_ess     <- ess_frac*S

  #if target_ess it too small
  if(target_ess<= 1){return(Inf)}

  #get my essy
  get_ess <- function(lambda){
    w <- .ppmn_stable_weights(loss, lambda)
    1/sum(w^2)}

  #if ess
  if(get_ess(max_lambda) > target_ess){return(max_lambda)}

  #find root
  uniroot(function(lambda)get_ess(lambda)-target_ess,
    interval = c(0, max_lambda))$root}
