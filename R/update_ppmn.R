#' Update PPMN function
#'
#' @param object A foundational or updated PPMN.
#' @param new_data New data set with all root vertices and least one child vertex present.
#' @param nsim Number of simulations used for updating.
#' @param level Generalized credibility level (default level = 0.9).
#' @param method The method being used to set the updating strength of the model. The default
#' method is 'kernel', The other method 'manual' uses the Kullenback-Leibler divergence applied
#' over the number of simulations to express the concentration strength relative to no change relative to the prior.
#' Here the default of the update strength is 95%.
#' @param up_strength The fraction of information that is allowed to concentrate in the generalized posterior
#'  when using the method 'manual'.
#' @param max_lambda Maximum lambda that can be set for Generalized Bayesian Updating (default max_lambda = 1000).
#' @param scale Whether each edge should be updated locally ("local") or together globally (default scale = "global").
#' @param fun_w Determines the method to weigh between multiple functions (multi-edge function) predicting a child.
#' The method to be chosen is "Exceedance", "Relative_loss" or "GBU" (default fun_w = "Relative_loss").
#' @param covar Covariance matrix that is calculated between model parameters (default covar = T).
#' @param seed Seed value 123.
#'
#' @description
#' The updating function uses loss-based Generalized Bayesian Updating (GBU) to update the prior to a generalized
#' posterior (often called GIBBS posterior). In simplicity it simulates parameters (particles) from prior model
#' These particles are used in a function to predict. Each prediction the mismatch with the data (loss) is calculated.
#' To updated the prior these losses each particles is then weighted inversely relative to the losses.
#'
#' @importFrom stats median cov2cor
#'
#' @export
update_ppmn <- function(object, new_data, nsim=3000, level=0.9,
                        method="kernel", up_strength=0.95, max_lambda=1000,
                        scale="global", fun_w = "Exceedance",
                        covar=T, seed=123){

  if(nsim<1){stop("number of simulations cannot be smaller than 1.")}
  if(level>0.99999){stop("level cannot be larger than .99999.")}
  if(level<0.00001){stop("level cannot be smaller than .00001.")}
  if(!is.null(up_strength)){
    if(up_strength>1){stop("update fraction cannot be larger than 1.")}
    if(up_strength<0){stop("update fraction cannot be smaller than 0.")}}
  if(!inherits(object,"ppmn")){class(object) <- c("ppmn","list")}

  #seed
  set.seed(seed)

  #constant
  eps <- .Machine$double.eps

  #check structure of the new data
  new_data <- new_data[,colnames(new_data) %in% object$Structure$exec_dag$order, drop=F]
  keep     <- apply(new_data,1,function(r) any(is.finite(as.numeric(r))))
  new_data <- new_data[keep,,drop=F]
  new_data <- new_data[,apply(new_data,2,function(x) any(is.finite(as.numeric(x)))), drop=F]

  #stop when columns are not provided
  if(nrow(new_data)==0){stop("After filtering new_data has 0 rows with finite values.")}

#########################
#predictions on new data#
#########################

  preds          <- predict_ppmn(object, new_data, nsim=nsim)

  pred_nodes_all <- setdiff(names(preds$Variance), object$Roots)
  obs_mat        <- as.matrix(new_data)
  obs_nodes      <- intersect(colnames(obs_mat), pred_nodes_all)

  if(length(obs_nodes)==0){stop("No observed child vertices in new_data match the model.")}

  lossandresid   <- .ppmn_loss_resid(obs_mat, preds, obs_nodes, nsim)
  scale          <- match.arg(scale,c("local","global"))

  edge_info      <- object$EdgeFunctions

#########################
#local summary per child#
#########################

  Fl   <- dim(lossandresid$loss1_fun)[4]
  n_J  <- length(obs_nodes)

  loss_fun <- array(NA,dim=c(nsim, n_J ,Fl), dimnames=list(NULL, obs_nodes,NULL))

  for(j in seq_along(obs_nodes)){

    for(f in seq_len(Fl)){

      L      <- lossandresid$loss1_fun[,j , ,f , drop=F]
      dim(L) <- c(nrow(obs_mat), nsim)

      loss_fun[ ,j ,f] <- apply(L,2,function(x){x <- x[is.finite(x)]
      if(length(x)==0){NA}else{median(x)}})}}

  function_particle_weights <- list()
  alpha                     <- list()
  weights                   <- list()
  lambda                    <- list()
  ESS                       <- c()

  #per edge summary
  if(scale=="local"){

    for(j in seq_along(obs_nodes)){

      node     <- obs_nodes[j]
      fun_rows <- edge_info[edge_info$child==node,,drop=F]
      ids      <- fun_rows$id
      fs       <- fun_rows$f

      pi_f <- object$FunctionWeights[[node]][ids]
      if(length(pi_f)==0 || any(!is.finite(pi_f)) || sum(pi_f)<=0){pi_f <- rep(1/length(ids),length(ids))}
      pi_f <- as.numeric(pi_f/sum(pi_f))

      L         <- loss_fun[,j,fs,drop=F]
      dim(L)    <- c(nsim, length(fs))

      loss_vec  <- as.numeric(L)
      prior_vec <- rep(pi_f,each=nsim)/nsim

      if(method=="kernel"){

        r1      <- lossandresid$resid1_fun[,j,,fs,drop=F]
        r0      <- lossandresid$resid0_fun[,j,,fs,drop=F]
        dim(r1) <- c(nrow(obs_mat),nsim*length(fs))
        dim(r0) <- c(nrow(obs_mat),nsim*length(fs))

        kl_target <- .ppmn_kl_calibrate(r1,r0,n_boot=2000,n_grid=3000)
        lambda_j  <- .ppmn_lambda_kl(loss_vec, kl_target, prior_vec, max_lambda)

      }else if(method=="manual"){

        kl_target <- log(1/(1-up_strength))
        lambda_j  <- .ppmn_lambda_kl(loss_vec, kl_target, prior_vec, max_lambda)

      }else{stop("Method should be 'kernel' or 'manual'.")}

      q <- .ppmn_stable_weights(loss_vec,lambda_j,prior_vec)
      q <- matrix(q,nrow=nsim,ncol=length(fs))

      if(fun_w=="GBU"){

        alpha[[node]]        <- colSums(q)
        names(alpha[[node]]) <- ids

      }else if(fun_w=="Relative_loss"){

        if(ncol(L)==1){
          alpha[[node]] <- 1
        }else{
        rs <- rowSums(L,na.rm=T)
        rs[!is.finite(rs) | rs==0] <- NA

        L_rel <- L/rs
        alpha[[node]]        <- 1-colMeans(L_rel,na.rm=T)}

        names(alpha[[node]]) <- ids

      }else if(fun_w=="Exceedance"){

        alpha[[node]] <- rep(0,ncol(L))

        for(s in seq_len(nrow(L))){

          z  <- L[s,]
          ok <- is.finite(z)

          if(any(ok)){

            best <- which(ok & z==min(z[ok]))
            alpha[[node]][best] <- alpha[[node]][best]+1/length(best)}
        }

        names(alpha[[node]]) <- ids

      }else{stop("fun_w should be 'GBU', 'Relative_loss' or 'Exceedance'.")}

      object$FunctionWeights[[node]] <- alpha[[node]]/sum(alpha[[node]])

      for(k in seq_along(ids)){

        w <- q[,k]
        if(sum(w)>0){w <- w/sum(w)}else{w <- rep(1/nsim,nsim)}
        function_particle_weights[[ids[k]]] <- w}

      weights[[node]] <- rowSums(q)
      lambda[[node]]  <- lambda_j
      ESS[node]       <- 1/sum(q^2)}

    ###################################
    #covar if scale is local bit wonky#
    ###################################

    logw_covar <- rep(0,nsim)

    for(node in obs_nodes){
      w <- weights[[node]]
      logw_covar <- logw_covar+log(pmax(w,eps))}

    logw_covar <- logw_covar-max(logw_covar)
    w_covar    <- exp(logw_covar)
    w_covar    <- w_covar/sum(w_covar)}

########################
#global combined update#
########################

  #is easy by summarizing the median per s for total loss and resid mat
  if(scale=="global"){

    object$EdgeFunctions$f

    loss_global <- apply(lossandresid$loss1, 3, function(z){
      z <- z[is.finite(z)]
      if(length(z)==0){NA}else{median(z)}})

    resid_global1 <- apply(lossandresid$resid1, 3, function(z){
      z <- z[is.finite(z)]
      if(length(z)==0){NA}else{median(z)}})

    resid_global0 <- apply(lossandresid$resid0, 3, function(z){
      z <- z[is.finite(z)]
      if(length(z)==0){NA}else{median(z)}})

    if(method=="kernel"){

      kl_target     <- .ppmn_kl_calibrate(resid_global1, resid_global0, n_boot=2000, n_grid=3000)
      global_lambda <- .ppmn_lambda_kl(loss_global, kl_target, prior = NULL, max_lambda)

    }else if(method=="manual"){

      kl_target     <- log(1/(1-up_strength))
      global_lambda <- .ppmn_lambda_kl(loss_global, kl_target, prior = NULL, max_lambda)

    }else{stop("Method should be 'kernel' or 'manual'.")}

    weights <- .ppmn_stable_weights(loss_global, global_lambda)
    lambda  <- global_lambda
    ESS     <- 1/sum(weights^2)
    w_covar <- weights}

#########################################
#funct weights alpha using global lambda#
#########################################

  for(j in seq_along(obs_nodes)){

    #select node id and func
    node     <- obs_nodes[j]
    fun_rows <- edge_info[edge_info$child==node,,drop=F]
    ids      <- fun_rows$id
    fs       <- fun_rows$f

    #previous function weights are the function prior
    pi_f <- object$FunctionWeights[[node]][ids]
    if(length(pi_f)==0||any(!is.finite(pi_f))||sum(pi_f)<=0){
      pi_f <- rep(1/length(ids),length(ids))}
      pi_f <- as.numeric(pi_f/sum(pi_f))

    #s,f loss matrix sxf
    L      <- loss_fun[,j,fs,drop=F]
    dim(L) <- c(nsim,length(fs))

    #gbu weighting usefull for similair functions horrible for
    #mixtures of functions as it uses the global lambda
    loss_vec  <- as.numeric(L)
    prior_vec <- rep(pi_f,each=nsim)/nsim

    q <- .ppmn_stable_weights(loss_vec,global_lambda,prior_vec)
    q <- matrix(q,nrow=nsim,ncol=length(fs))

    if(fun_w=="GBU"){

      alpha[[node]]        <- colSums(q)
      names(alpha[[node]]) <- ids

    }else if(fun_w=="Relative_loss"){

      if(ncol(L)==1){
        alpha[[node]] <- 1
      }else{
        rs <- rowSums(L,na.rm=T)
        rs[!is.finite(rs) | rs==0] <- NA

        L_rel <- L/rs
        alpha[[node]]        <- 1-colMeans(L_rel,na.rm=T)}

      names(alpha[[node]]) <- ids

    }else if(fun_w=="Exceedance"){

      alpha[[node]] <- rep(0,ncol(L))

      for(s in seq_len(nrow(L))){

        z  <- L[s,]
        ok <- is.finite(z)

        if(any(ok)){

          best                <- which(ok & z==min(z[ok]))
          alpha[[node]][best] <- alpha[[node]][best]+1/length(best)}}

      names(alpha[[node]]) <- ids

    }else{stop("fun_w should be 'GBU', 'Relative_loss' or 'Exceedance'.")}

    object$FunctionWeights[[node]] <- alpha[[node]]/sum(alpha[[node]])

    #Particle weights for each multi-edge function not within each
    for(k in seq_along(ids)){

      w <- q[,k]
      w <- w/sum(w)

      function_particle_weights[[ids[k]]] <- w}}

###################################
#sumarize params per edge function#
###################################

  function_summaries <- list()
  function_draws     <- list()

  for(id in names(object$Parameters)){

    node  <- edge_info$child[match(id,edge_info$id)]
    draws <- preds$ParametersPerFunction[[id]]

    if(is.null(draws)){
      stop("No parameter draws found for edge function: ",id)}

    draws <- do.call(rbind,draws)
    draws <- as.matrix(draws)
    colnames(draws) <- names(object$Parameters[[id]])

    if(scale=="local"){
      w <- function_particle_weights[[id]]
      if(is.null(w)){w <- rep(1/nsim,nsim)}
    }else{
      w <- weights}

    w <- w/sum(w)

    mu     <- colSums(draws*w)
    diffmu <- sweep(draws,2,mu,"-")
    se     <- sqrt(colSums((diffmu^2)*w))
    llul   <- apply(draws,2,function(x){.ppmn_wquantile(x, w, c((1-level)/2,1-(1-level)/2))})

    function_summaries[[id]] <- list(mu=mu,se=se,ll=llul[1,],ul=llul[2,])
    function_draws[[id]]     <- draws}

  ###################
  #derive covariance#
  ###################

  sigma_weighted <- NULL

  if(covar==T){

    if(covar==T && scale=="local"){warning("Covariance is derived from weighted edge lambda's, not one matrix lambda.")}

    params                     <- do.call(cbind,function_draws)
    params                     <- as.matrix(params)
    params[!is.finite(params)] <- 0

    mu_par      <- colSums(params*w_covar)
    diff_mu_par <- sweep(params,2,mu_par,"-")
    diff_mu_par[!is.finite(diff_mu_par)] <- 0

    sigma_weighted <- t(diff_mu_par)%*%sweep(diff_mu_par,1,w_covar,"*")
    sigma_weighted <- sigma_weighted+diag(eps,ncol(sigma_weighted))}

  ##########################
  #move summaries to object#
  ##########################

  old_params <- object$Parameters

  for(id in names(old_params)){

    summ    <- function_summaries[[id]]
    sub_old <- old_params[[id]]

    for(b in names(sub_old)){

      if(!b %in% names(summ$mu)){next}

      sub_old[[b]]["mu"] <- as.numeric(summ$mu[[b]])
      sub_old[[b]]["se"] <- as.numeric(summ$se[[b]])
      sub_old[[b]]["ll"] <- as.numeric(summ$ll[[b]])
      sub_old[[b]]["ul"] <- as.numeric(summ$ul[[b]])}

    old_params[[id]] <- sub_old}

  object$Update           <- object$Update+1
  object$`Old parameters` <- object$Parameters
  object$Parameters       <- old_params

  object$Info <- list(scale=scale, method=method, weights=weights,
                      function_particle_weights=function_particle_weights,
                      function_weights=object$FunctionWeights,
                      lambda=unlist(lambda), ESS=ESS)

  ##########################
  #rebuild covariance Sigma#
  ##########################

  theta_names <- unlist(Map(function(id,pars) paste0(id,"_",names(pars)),
                            names(object$Parameters),object$Parameters))

  theta_var <- unlist(lapply(object$Parameters,function(e){
    sapply(e,function(p) as.numeric(p["se"])^2)}))

  theta_var[!is.finite(theta_var) | theta_var<1e-8] <- 1e-8

  Sigma_new <- diag(theta_var)
  dimnames(Sigma_new) <- list(theta_names,theta_names)

  if(covar==T && !is.null(sigma_weighted)){

    rownames(sigma_weighted) <- theta_names
    colnames(sigma_weighted) <- theta_names

    sigma_weighted[!is.finite(sigma_weighted)] <- 1e-8
    diag(sigma_weighted) <- pmax(diag(sigma_weighted),eps)

    R                <- cov2cor(sigma_weighted)
    R[!is.finite(R)] <- 0
    diag(R)          <- 1

    sd_new                                 <- sqrt(pmax(theta_var,eps))

    object$Sigma                           <- diag(sd_new)%*%R%*%diag(sd_new)
    object$Sigma[!is.finite(object$Sigma)] <- 0
    diag(object$Sigma)                     <- pmax(diag(object$Sigma),eps)
    dimnames(object$Sigma)                 <- list(theta_names,theta_names)

  }else{
    object$Sigma <- Sigma_new}

  #############################
  #update residual diagnostics#
  #############################

  object <- .ppmn_resid_diagnostic(object,new_data)

  ########################
  #rebuild parameter table#
  ########################

  par_tab <- lapply(seq_len(nrow(edge_info)),function(i){

    id       <- edge_info$id[i]
    sub      <- object$Parameters[[id]]
    edge     <- edge_info$edge[i]
    edge_fun <- edge_info$fun[i]
    dep      <- edge_info$child[i]
    indep    <- object$Structure$predict_dag$parents[[id]]

    if(edge_fun %in% c("asymptotic","sigmoidal","mvlogistic","gompertz","gaussian")){
      indep <- c(NA,NA,indep)
    }else if(edge_fun=="class"){
      indep <- c(NA,indep,NA,indep,NA)
    }else{
      indep <- c(NA,indep)}

    cbind(
      edge=edge,
      `edge function`=edge_fun,
      dependent=dep,
      independent=indep,
      parameters=names(sub),
      do.call(rbind,sub))})

  object$`Parameter table` <- do.call(rbind.data.frame,par_tab)
  rownames(object$`Parameter table`) <- NULL

  return(object)}
