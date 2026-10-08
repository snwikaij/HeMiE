
#' Helper transformation function
#'
#' @param x_raw Raw observations.
#' @param x_hat Predictions on x_raw.
#' @param type Type of distribution for x_raw.
#' @param MAD Median absolute deviation.
#'
#' @keywords internal
#'
.ppmn_transform_loss <- function(x_raw, x_hat, type, MAD=NULL){

  #constant
  eps <- .Machine$double.eps

######################
#loss and resid funcs#
######################

  if(type %in% c("binary", "proportional")){
    p_hat    <- pmin(pmax(x_hat, eps), 1-eps)
    r        <- (x_raw-p_hat)/sqrt(p_hat*(1-p_hat)+eps)
    l        <- .ppmn_logistic_loss(x_raw, x_hat)
  }else if(type %in% "discrete"){
    r        <- 2*(sqrt(x_raw+3/8)-sqrt(pmax(x_hat, eps)+3/8))
    l        <- r^2
  }else if(type=="continuous"){
    if(is.null(MAD)){stop("Why is MAD not there?")}
    r        <- (x_raw-x_hat)/MAD
    l        <- r^2
  }else{stop("Unknown response type: ", type)}


  #check for inf and na and stuff
  bad <- (!is.finite(l) | !is.finite(r))

  if(any(bad)){
    l[bad] <- NA
    r[bad] <- NA}

  list(loss=l, residuals=r)}
