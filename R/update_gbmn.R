#' Update GBMN function
#'
#' @param object A foundational or updated GBMN.
#' @param new_data New data set with all root vertices and least one child vertex present.
#' @param nsim Number of simulations used for updating.
#' @param level Generalized credibility level (default level = 0.9).
#' @param method The method being used to set the updating strength of the model (default =
#' = 'hold-amis').
#' The method 'manual' concentrates by selecting the best x% of the simulations.
#' The method coverage calculates the coverage of the expected value using bootstrapped subset.
#' The method 'holdout' is a version that looks similar to SafeBayes.
#' The method 'amis' uses Adaptive Multiple Importance Sampling.
#' The method 'hold-amis' uses a fixed 5 repetitions for the 'holdout' method
#' and uses the derived learning rate as a fixed value for amis.
#' @param up_character A parameter controlling the updating method.
#' For method "hold-amis", it specifies the number adaptive iterations (default = 5).
#' For method "amis", it specifies the number adaptive iterations (default = 5).
#' For method "holdout", it specifies the number of repetitions (default = 20).
#' For method "manual", it specifies the fraction of best-performing simulations (default = NULL).
#' For method "coverage", it specifies the number of bootstraps (default = 50).
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
                        method="hold-amis", up_character=NULL, max_lambda=200,
                        diagnostics=T, covar=T, seed=123){

  if(nsim<1){stop("number of simulations cannot be smaller than 1.")}
  if(level>0.99999){stop("level cannot be larger than 0.99999.")}
  if(level<0.00001){stop("level cannot be smaller than 0.00001.")}

  if(method=="hold-amis"){
    if(is.null(up_character)){up_character <- 5}
    if(up_character<1){
      stop("The update character for the method 'hold-amis' is the number of adaptive repetitions and cannot be smaller than 1.")}
    niter      <- as.integer(up_character)
    first_nsim <- 1000
    split_nsim <- floor(nsim/niter)
    if(split_nsim<50){
      stop("At least 50 particles per AMIS iteration are required.")}}
  if(method=="amis"){
    if(is.null(up_character)){up_character <- 5}
    if(max_lambda==200){message("AMIS uses a (heuristic) lambda = 1 by default. Set max_lambda to change the learning rate.")
      max_lambda <- 1}
      if(up_character<1){
        stop("The update character for the method 'amis' is the number of addaptive repetitions and cannot be smaller than 1.")}
         niter      <- as.integer(up_character)
         split_nsim <- floor(nsim/(niter+1))
         first_nsim <- nsim-niter*split_nsim
   if(split_nsim<50){
      stop("Heuristic boundary of at least 50 particles per iteration is required.")}}else if(method!="hold-amis"){
        first_nsim <- nsim}
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

  preds          <- predict_gbmn(object, new_data, nsim=first_nsim)

  pred_nodes_all <- setdiff(names(preds$Variance), object$Roots)
  obs_mat        <- as.matrix(new_data)
  obs_nodes      <- intersect(colnames(obs_mat), pred_nodes_all)

  if(length(obs_nodes)==0){stop("No observed child vertices in new_data match the model.")}

  lossandresid   <- .gbmn_loss_resid(obs_mat, preds, obs_nodes, first_nsim)

  edge_info      <- object$EdgeFunctions

######
#loss#
######

  loss <- apply(lossandresid$loss1, 3, function(z){
    z <- z[is.finite(z)]
    if(length(z)==0){NA}else{mean(z)}})

  if(method=="manual"){

    lambda        <- .gbmn_lambda_level(loss, ess_frac=up_character, max_lambda=max_lambda)

  }else if(method=="holdout"){

    lambda        <- .gbmn_lambda_fast(lossandresid$loss1, n_rep=up_character, max_lambda=max_lambda)

  }else if(method=="coverage"){

    lambda        <- .gbmn_coverage_fast(lossandresid$loss1, preds=preds,
                                 obs_mat=obs_mat, obs_nodes=obs_nodes,
                                 max_lambda=max_lambda,
                                 level=level, nsim=nsim, nboot=up_character)
  }else if(method=="amis"){

    lambda        <- max_lambda

  }else if(method=="hold-amis"){

    lambda        <- .gbmn_lambda_fast(lossandresid$loss1, n_rep=5, max_lambda=max_lambda)

  }else{stop("Method should be 'hold-amis', 'amis', 'holdout' 'coverage' or 'manual'.")}

  if(method %in% c("amis", "hold-amis")){

    amis_result   <- .gbmn_amis_easy(object=object, preds=preds, loss=loss,
                                 new_data=new_data, obs_nodes=obs_nodes,
                                 first_mad=lossandresid$MAD_used, obs_mat=obs_mat,
                                 niter=up_character, lambda=lambda, nsim=split_nsim)

    preds_post <- amis_result$preds
    weights    <- amis_result$weights

  }else{

    preds_post <- preds
    weights    <- .gbmn_stable_weights(loss, lambda)
  }

  ESS <- 1/sum(weights^2)

  if(method!="amis" && lambda>=max_lambda){warning("Lambda is max_lambda, consider increasing max_lambda.")}
  if(ESS < length(weights)*0.05){warning("The ESS is smaller than 5% of the simulations.")}

#######################
#summarize all results#
#######################

  particle_summ <- .gbmn_summarize_particles(object=object, preds=preds_post, weights=weights,
                                             obs_nodes=obs_nodes, level=level, covar=covar)

  #extract posterior res
  function_summaries        <- particle_summ $function_summaries
  function_particle_weights <- particle_summ $function_particle_weights
  sigma_weighted            <- particle_summ $sigma_weighted

  #update function alpha weights
  object$FunctionWeights    <- particle_summ$function_weights

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

  object$ParameterTable           <- do.call(rbind.data.frame, par_tab)
  rownames(object$ParameterTable) <- NULL

return(object)}
