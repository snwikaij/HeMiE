#' Helper function that defines and predicts with different functions
#'
#' @param func Name of the function.
#' @param beta Parameters of the functions.
#' @param trans Transformation type applied over all independent variables.
#' @param x matrix of independent variables.
#'
#' @keywords internal
#'
#' @importFrom stats plogis dnorm
#'
#' @export
.ppmn_predict_edge <- function(func, beta, trans, x){

  #transformation functions under estimation
  if(!is.na(trans)){
    if(trans == "log"){
      x[x <= 0] <- NA
      x <- log(x)}
    if(trans == "sqrt"){
      x[x < 0] <- NA
      x <- sqrt(x)}
  if(trans == "log10"){
    x[x < 0] <- NA
    x <- log10(x)}}

  #link function of standard (g)lms
  if(func %in% c("identity", "log", "logit", "log10")){eta <- beta[1] + drop(x %*% beta[-1])}
  if(func == "identity"){return(eta)}
  if(func == "log"){eta <- pmax(pmin(eta, 30), -30); return(exp(eta))}
  if(func == "logit"){return(plogis(eta))}
  if(func == "log10"){eta <- pmax(pmin(eta, 30), -30); return(10^(eta))}

  #non standard models that are useful
  if(func == "asymptotic"){
    b0 <- beta[1]
    b1 <- beta[2]
    b2 <- beta[3]

    x  <- as.matrix(x)
    x1 <- x[, 1]

  if(ncol(x) > 1){gamma <- drop(x[, -1, drop = FALSE] %*% beta[-(1:3)])}else{gamma <- 0}
    eta <- b2 * x1 + gamma

    return(b0 - b1 * exp(-eta))}

  if(func == "mvlogistic"){

    #asymptote
    b0 <- beta[1]

    #intercept
    b1 <- beta[2]

    x <- as.matrix(x)

    #coefs
    params <- beta[-(1:2)]

    if(ncol(x) != length(params)){
      stop("Number of predictor parameters does not match number of predictors.")}

    eta <- b1+drop(x %*% params)

    return(b0*plogis(eta))}

  if(func == "gompertz"){
    b0 <- beta[1]
    b1 <- beta[2]
    b2 <- beta[3]

    x <- as.matrix(x)
    x1 <- x[, 1]

    # sum(bj * xj) for j >= 3 corresponds to x columns
  if (ncol(x) > 1) {gamma <- drop(x[, -1, drop = FALSE] %*% beta[-(1:3)])}else{gamma <- 0}

    eta <- b1 + gamma - b2 * x1
    return(b0 * exp(-exp(eta)))}

  if(func == "sigmoidal"){
    b0     <- beta[1]
    b1     <- beta[2]
    b2     <- beta[3]

    x <- as.matrix(x)
    x1 <- x[, 1]

  if (ncol(x) > 1){gamma <- drop(x[, -1, drop = FALSE] %*% beta[-(1:3)])}else{gamma <- 0}
    eta <- (x1 - b1 + gamma) / b2
    return(b0 / (1 + exp(-eta)))}

  if(func == "gaussian"){return(beta[1] * exp(-0.5 * ((x - beta[2]) / beta[3])^2))}

  if(func == "class"){p <- exp(dnorm(x, beta[1], beta[2], log = T) - dnorm(x, beta[3], beta[4], log = T))
  return(as.numeric((p * beta[5]) / ((p * beta[5]) + (1 - beta[5]))))}
  }
