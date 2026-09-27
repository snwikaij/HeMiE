#' Helper function to calculate loss and residuals
#'
#' @param obs_mat Observed values matrix.
#' @param preds Predictions on this observed matrix.
#' @param obs_nodes Observed nodes of the matrix.
#' @param nsim Number of simulations.
#' @param mad_constant Constant.
#'
#' @importFrom stats mad median
#'
.ppmn_loss_resid <- function(obs_mat, preds, obs_nodes, nsim, mad_constant=1){

  #constant
  eps     <- .Machine$double.eps

  #number of rows
  i       <- nrow(obs_mat)

  #number of nodes
  n_vert  <- length(obs_nodes)

  #number of functions
  n_func  <- dim(preds$FunctionPredictions)[4]

  #child
  child <- dimnames(preds$FunctionPredictions)[[2]]

  #combined arrays [i,j,s]
  loss1  <- array(NA,dim=c(i,n_vert,nsim),dimnames=list(NULL,obs_nodes,NULL))
  resid1 <- array(NA,dim=c(i,n_vert,nsim),dimnames=list(NULL,obs_nodes,NULL))
  resid0 <- array(NA,dim=c(i,n_vert,nsim),dimnames=list(NULL,obs_nodes,NULL))

  #function arrays [i,j,s,f]
  loss1_fun  <- array(NA,dim=c(i,n_vert,nsim,n_func),dimnames=list(NULL,obs_nodes,NULL,NULL))
  resid1_fun <- array(NA,dim=c(i,n_vert,nsim,n_func),dimnames=list(NULL,obs_nodes,NULL,NULL))
  resid0_fun <- array(NA,dim=c(i,n_vert,nsim,n_func),dimnames=list(NULL,obs_nodes,NULL,NULL))

  for(j in seq_along(obs_nodes)){

    node  <- obs_nodes[j]
    x_raw <- obs_mat[,node]
    type  <- .ppmn_detect_type(x_raw)

    #derive mad directly from prior predictions to scale the residuals
    x_exp <- preds$Expected[[node]]

    if(type=="continuous"){
      r0  <- x_raw-x_exp
      MAD <- mad(r0, center=median(r0, na.rm=T), constant=mad_constant, na.rm=T)
      if(!is.finite(MAD)||MAD==0){MAD <- eps}
    }else{
      MAD <- 1}

    #all predictions
    xhat_mat <- do.call(cbind,lapply(seq_len(nsim), function(s){preds$Variance[[node]][[s]]}))

    #median prior predictions
    x0       <- apply(xhat_mat, 1, median, na.rm=T)

    for(s in seq_len(nsim)){

      #data-model data variance
      a1 <- .ppmn_transform_loss(x_raw, xhat_mat[,s], type, MAD)

      #model-mean(model) model variance
      a0 <- .ppmn_transform_loss(x0, xhat_mat[,s], type, MAD)

      #loss and resid excluding the function
      loss1[,j,s]  <- a1$loss
      resid1[,j,s] <- a1$residuals
      resid0[,j,s] <- a0$residuals}

    #########################
    #separate edge functions#
    #########################

    #if node/vertex is a child select out the predictions of the edge functions
    loc_child   <- match(node, child)
    fun_rows    <- preds$EdgeFunctions[preds$EdgeFunctions$child==node,,drop=F]

    for(k in seq_len(nrow(fun_rows))){

      f         <- fun_rows$f[k]

      xfun      <- preds$FunctionPredictions[,loc_child,,f, drop=F]
      dim(xfun) <- c(i, nsim)
      x0_fun    <- apply(xfun, 1, median,na.rm=T)

      for(s in seq_len(nsim)){

        a1 <- .ppmn_transform_loss(x_raw, xfun[,s], type, MAD)
        a0 <- .ppmn_transform_loss(x0_fun, xfun[,s], type, MAD)

        loss1_fun[,j,s,f]  <- a1$loss
        resid1_fun[,j,s,f] <- a1$residuals
        resid0_fun[,j,s,f] <- a0$residuals}}
  }

  list(loss1=loss1,resid1=resid1,resid0=resid0,
       loss1_fun=loss1_fun,resid1_fun=resid1_fun,resid0_fun=resid0_fun)}
