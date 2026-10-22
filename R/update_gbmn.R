#' Update GBMN function
#'
#' @param object A foundational or updated GBMN.
#' @param new_data New data set with all root vertices and least one child vertex present.
#' @param nsim Number of simulations used for updating.
#' @param level Generalized credibility level (default level = 0.9).
#' @param method The method being used to set the updating strength of the model. The default
#' method is 'holdout', The method 'manual' concentrates by selecting the best x% of the simulations. The method coverage
#' calculates the coverage of the expected value using bootstrapped subset. The method 'holdout' is a version
#' that looks similar to SafeBayes.
#' @param up_character A parameter controlling the updating method.
#' For method "holdout", it specifies the number of repetitions (default = 20).
#' For method "manual", it specifies the fraction of best-performing simulations.
#' For method "coverage", it specifies the number of bootstraps.
#' @param max_lambda Maximum lambda that can be set for Generalized Bayesian Updating (default max_lambda = 200).
#' @param diagnostics If residual diagnostics should be returned (default diagnostics = T).
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
update_gbmn <- function(object, new_data, nsim=3000, level=0.9,
                        method="holdout", up_character=NULL, max_lambda=200,
                        diagnostics=T, covar=T, seed=123){

  if(nsim<1){stop("number of simulations cannot be smaller than 1.")}
  if(level>0.99999){stop("level cannot be larger than 0.99999.")}
  if(level<0.00001){stop("level cannot be smaller than 0.00001.")}

  if(method=="holdout"){
    if(is.null(up_character)){up_character <- 20}
    if(up_character<1){
      stop("The update character for the method 'holdout' is the number of hold-out repetitions and cannot be smaller than 1.")}}
  if(method=="coverage"){
    if(is.null(up_character)){up_character <- 50}
    if(up_character<1){
      stop("The update character for the method 'coverage' is the number of bootstraps and cannot be smaller than 1.")}}
  if(method=="manual"){
    if(is.null(up_character) || up_character>1 || up_character<0){
      stop("The update character for method 'manual' is the best percentage of particles and cannot be smaller than 0 or larger than 1.")}}

  if(!inherits(object,"gbmn")){class(object) <- c("gbmn","list")}

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

  preds          <- predict_gbmn(object, new_data, nsim=nsim)

  pred_nodes_all <- setdiff(names(preds$Variance), object$Roots)
  obs_mat        <- as.matrix(new_data)
  obs_nodes      <- intersect(colnames(obs_mat), pred_nodes_all)

  if(length(obs_nodes)==0){stop("No observed child vertices in new_data match the model.")}

  lossandresid   <- .gbmn_loss_resid(obs_mat, preds, obs_nodes, nsim)

  edge_info      <- object$EdgeFunctions

######
#loss#
######

  loss <- apply(lossandresid$loss1, 3, function(z){
    z <- z[is.finite(z)]
    if(length(z)==0){NA}else{mean(z)}})

  #left from previous coding might be useful in the future?
  #resid <- apply(lossandresid$resid1, 3, function(z){
  #  z <- z[is.finite(z)]
  #  if(length(z)==0){NA}else{median(z)}})

  if(method=="manual"){

    lambda        <- .gbmn_lambda_level(loss, ess_frac=up_character, max_lambda=max_lambda)

  }else if(method=="holdout"){

    lambda        <- .gbmn_lambda_fast(lossandresid$loss1, n_rep=up_character, max_lambda=max_lambda)

  }else if(method=="coverage"){

    lambda        <- .gbmn_coverage_fast(lossandresid$loss1, preds=preds,
                                 obs_mat=obs_mat, obs_nodes=obs_nodes,
                                 max_lambda=max_lambda,
                                 level=level, nsim=nsim, nboot=up_character)

  }else{stop("Method should be 'holdout' 'coverage' or 'manual'.")}

  weights                   <- .gbmn_stable_weights(loss, lambda)
  ESS                       <- 1/sum(weights^2)
  function_particle_weights <- list()

  if(lambda>=max_lambda){warning("Lambda is max_lambda, consider increasing max_lambda.")}
  if(ESS<(nsim*0.05)){warning("The ESS is smaller than 5% of the simulations, indicating the posterior concentration may be too strong.")}

############################################
#nodes that contribute to observed vertices#
############################################

  #select graph
  g <- object$Structure$exec_dag$graph

  #select relevant nodes
  relevant_nodes <- unique(unlist(lapply(obs_nodes, function(node){
        igraph::as_ids(igraph::subcomponent(g,node,mode = "in"))})))

  #retain childs that ahave edge functions
  relevant_nodes <- intersect(relevant_nodes, unique(edge_info$child))

  #and that have function selections from predict_gbmn()
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

    #use only particles in which function was active
    w <- function_particle_weights[[id]]

    #if this function was never selected, there is no
    #information with which to update its parameters
    if(is.null(w) || length(w)==0 || sum(w,na.rm=TRUE)<=0){

      old <- object$Parameters[[id]]

      function_summaries[[id]] <- list(
        mu = vapply(old,function(x) as.numeric(x["mu"]),numeric(1)),
        se = vapply(old,function(x) as.numeric(x["se"]),numeric(1)),
        ll = vapply(old,function(x) as.numeric(x["ll"]),numeric(1)),
        ul = vapply(old,function(x) as.numeric(x["ul"]),numeric(1)))

        function_draws[[id]] <- draws

        next}

    w[!is.finite(w)] <- 0
    w      <- w/sum(w)

    mu     <- colSums(draws*w)
    diffmu <- sweep(draws,2,mu,"-")
    se     <- sqrt(colSums((diffmu^2)*w))
    llul   <- apply(draws,2,function(x){.gbmn_wquantile(x, w, c((1-level)/2,1-(1-level)/2))})

    function_summaries[[id]] <- list(mu=mu,se=se,ll=llul[1,],ul=llul[2,])
    function_draws[[id]]     <- draws}

###################
#derive covariance#
###################

 sigma_weighted <- NULL

 if(covar==T){

    #create param matrix [S,K]
    params <- do.call(cbind, function_draws)
    params <- as.matrix(params)

    #all weights
    mu_vec <- colSums(params*weights)

    #center parameters
    diff_mu_par <- sweep(params, 2, mu_vec, "-")

    # Weighted covariance
    sigma_weighted <- t(diff_mu_par) %*% diag(weights) %*% diff_mu_par

    #add small constant to prevent issues
    sigma_weighted <- sigma_weighted+diag(eps, ncol(sigma_weighted))}

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
  object$OldParameters    <- object$Parameters
  object$Parameters       <- old_params

  object$Info <- list(method=method, weights=weights,
                      function_particle_weights=function_particle_weights,
                      function_weights=object$FunctionWeights,
                      lambda=lambda, ESS=ESS)

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

  if(diagnostics){object <- .gbmn_resid_diagnostic(object, new_data)}

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

  object$ParameterTable <- do.call(rbind.data.frame, par_tab)
  rownames(object$ParameterTable) <- NULL

return(object)}
