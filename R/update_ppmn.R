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

  edge_info      <- object$EdgeFunctions

######
#loss#
######

  loss <- apply(lossandresid$loss1, 3, function(z){
    z <- z[is.finite(z)]
    if(length(z)==0){NA}else{median(z)}})

  resid_1 <- apply(lossandresid$resid1, 3, function(z){
    z <- z[is.finite(z)]
    if(length(z)==0){NA}else{median(z)}})

  resid_0 <- apply(lossandresid$resid0, 3, function(z){
    z <- z[is.finite(z)]
    if(length(z)==0){NA}else{median(z)}})

  if(method=="kernel"){

    kl_target     <- .ppmn_kl_calibrate(resid_1, resid_0, n_boot=2000, n_grid=3000)
    lambda        <- .ppmn_lambda_kl(loss, kl_target, prior = NULL, max_lambda)

  }else if(method=="manual"){

    kl_target     <- log(1/(1-up_strength))
    lambda <- .ppmn_lambda_kl(loss, kl_target, prior = NULL, max_lambda)

  }else{stop("Method should be 'kernel' or 'manual'.")}

  weights <- .ppmn_stable_weights(loss, lambda)
  ESS     <- 1/sum(weights^2)
  function_particle_weights <- list()

############################################
#nodes that contribute to observed vertices#
############################################

  #select graph
  g <- object$Structure$exec_dag$graph

  relevant_nodes <- unique(unlist(lapply(obs_nodes, function(node){
        igraph::as_ids(igraph::subcomponent(g,node,mode = "in"))})))

    #retain childs that ahave edge functions
    relevant_nodes <- intersect(relevant_nodes, unique(edge_info$child))

    #and that have function selections from predict_ppmn()
    relevant_nodes <- intersect(relevant_nodes, colnames(preds$FunctionSelection))

###############################
#update alpha per child vertex#
###############################

  for(node in relevant_nodes){

    ids <- edge_info$id[edge_info$child == node]

    #function selected in each simulation
    selected <- preds$FunctionSelection[, node]

##################################
#alpha_m(new) = sum_s w(s) z_m(s)#
##################################

  alpha_new <- vapply(ids, function(id){
    active <- !is.na(selected) & selected == id
    sum(weights[active], na.rm = TRUE)}, numeric(1))

    names(alpha_new) <- ids

    #numerical normalization
    if(sum(alpha_new) > 0){

    alpha_new <- alpha_new/sum(alpha_new)

    }else{
    #if nothing usable was sampled, #retain previous function weights
    alpha_new <- object$FunctionWeights[[node]][ids]

    if(length(alpha_new) == 0||any(!is.finite(alpha_new))||sum(alpha_new) <= 0){

    alpha_new <- rep(1/length(ids), length(ids))}

    alpha_new <- alpha_new/sum(alpha_new)

    names(alpha_new) <- ids}

    object$FunctionWeights[[node]] <- alpha_new

############################
#for all selected functions#
############################

   for(id in ids){
    active <- !is.na(selected) & selected == id
    w_id   <- weights*as.numeric(active)

   if(sum(w_id) > 0){
    function_particle_weights[[id]] <-w_id/sum(w_id)

   }else{
    #function not selected in any particle
    function_particle_weights[[id]] <- NULL}}}

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

      #global update:
      #use only particles in which this function was active
      w <- function_particle_weights[[id]]

      #if this function was never selected, there is no
      #information with which to update its parameters
      if(is.null(w) || length(w)==0 || sum(w,na.rm=TRUE)<=0){

        old <- object$Parameters[[id]]

        function_summaries[[id]] <- list(
          mu = vapply(old,function(x) as.numeric(x["mu"]),numeric(1)),
          se = vapply(old,function(x) as.numeric(x["se"]),numeric(1)),
          ll = vapply(old,function(x) as.numeric(x["ll"]),numeric(1)),
          ul = vapply(old,function(x) as.numeric(x["ul"]),numeric(1))
        )

        function_draws[[id]] <- draws

        next
      }

    w[!is.finite(w)] <- 0
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

    params <- do.call(cbind,function_draws)
    params <- as.matrix(params)

    #which function and child belongs to each parameter
    npars       <- sapply(function_draws, ncol)
    fun_par     <- rep(names(function_draws), npars)
    child_par   <- edge_info$child[match(fun_par, edge_info$id)]

    #parameter active or inactive within each simulation
    active_par <- matrix(F,nrow=nsim,ncol=ncol(params))

    for(k in seq_len(ncol(params))){

      selected <- preds$FunctionSelection[,child_par[k]]

      active_par[,k] <- !is.na(selected) & selected==fun_par[k]}

      #use conditional updated parameter means
      mu_par <- unlist(lapply(names(function_summaries),function(id){function_summaries[[id]]$mu}))
      mu_par <- as.numeric(mu_par)

      #deviation only contributes when function was selected
      diff_mu_par                            <- sweep(params, 2, mu_par,"-")
      diff_mu_par[!active_par]               <- 0
      diff_mu_par[!is.finite(diff_mu_par)]   <- 0

      sigma_weighted <- t(diff_mu_par)%*%sweep(diff_mu_par,1,weights,"*")
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

  object <- .ppmn_resid_diagnostic(object, new_data)

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
