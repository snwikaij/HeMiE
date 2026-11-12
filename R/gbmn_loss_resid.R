#' Helper function to calculate loss and residuals
#'
#' @param obs_mat Observed values matrix.
#' @param preds Predictions on this observed matrix.
#' @param obs_nodes Observed nodes of the matrix.
#' @param nsim Number of simulations.
#' @param mad_constant Constant.
#'
#' @keywords internal
#'
#' @importFrom stats mad median
#'
.gbmn_loss_resid <- function(obs_mat, preds, obs_nodes, nsim, mad_constant=1.4826, MAD_used=NULL){

  #constant
  eps      <- .Machine$double.eps

  #number of rows
  i        <- nrow(obs_mat)

  #number of nodes
  n_vert   <- length(obs_nodes)

  #store null mad
  null_mad <- is.null(MAD_used)

  #create named vector
  if(null_mad){MAD_used <- setNames(rep(NA, n_vert), obs_nodes)}

  #number of functions
  n_func  <- dim(preds$FunctionPredictions)[4]

  #child
  child <- dimnames(preds$FunctionPredictions)[[2]]

  #combined arrays [i,j,s]
  loss1  <- array(NA, dim=c(i, n_vert, nsim), dimnames=list(NULL, obs_nodes, NULL))
  resid1 <- array(NA, dim=c(i, n_vert, nsim), dimnames=list(NULL, obs_nodes, NULL))

  for(j in seq_along(obs_nodes)){

    node  <- obs_nodes[j]
    x_raw <- obs_mat[,node]
    type  <- .gbmn_detect_type(x_raw)

    #derive mad directly from prior predictions to scale the residuals
    x_exp <- preds$Expected[[node]]

    #for continuous variables.
    if(type=="continuous"){

      #use mad provided if null
      if(null_mad){

      r0  <- x_raw-x_exp
      MAD <- mad(r0, center=median(r0, na.rm=T), constant=mad_constant, na.rm=T)

    #if mini set 0
      if(!is.finite(MAD) || MAD==0){MAD <- eps}

    #re-use mad from previous (I hope)
    }else{

      MAD <- MAD_used[[node]]}

    }else{
      MAD <- 1
    }

    #store mad
    MAD_used[node] <- MAD

    #all predictions
    xhat_mat <- do.call(cbind,lapply(seq_len(nsim), function(s){preds$Variance[[node]][[s]]}))

    for(s in seq_len(nsim)){

      #data-model data variance
      a1 <- .gbmn_transform_loss(x_raw, xhat_mat[,s], type, MAD)

      #loss and resid excluding the function
      loss1[,j,s]  <- a1$loss
      resid1[,j,s] <- a1$residuals}
  }

  return(list(loss1=loss1, resid1=resid1, MAD_used=MAD_used))}
