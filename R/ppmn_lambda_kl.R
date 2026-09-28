#' Helper function for Kullenback-Leiber optimizer for lambda
#'
#' @param loss The loss.
#' @param kl_target KL-divergence target value.
#' @param prior Prior.
#' @param upper Upper boundary for lambda.
#'
#' @keywords internal
#'
#' @importFrom stats uniroot
#'
.ppmn_lambda_kl <- function(loss, kl_target, prior=NULL, upper=max_lambda){

  f_up <- .ppmn_kl_optim(upper, loss, kl_target, prior)
  if(!is.finite(f_up)||f_up<0){return(upper)}

  uniroot(.ppmn_kl_optim, interval=c(0,upper), loss=loss,
          kl_target=kl_target, prior=prior)$root}
