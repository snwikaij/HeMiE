#' Helper function for Adapted Multiple Importance Sampling
#'
#' @param object An foundational or updated gbmn object.
#' @param preds Predictions of that object.
#' @param loss Loss based on those predictions.
#' @param new_data The new data used for that loss.
#' @param first_mad The Maximum Absolute Deviation of that prediction
#' @param obs_nodes The observed vertices (still inconsistent with terminology)
#' @param obs_mats Observed matrix after cleaning new data.
#' @param niter Number of adaptive iterations (T).
#' @param lambda Lambda which is by default 1.
#' @param nsim Number of simulations per iteration.
#'
#' @description
#' This is based on the article "Adaptive Multiple Importance Sampling" by Cornuet et al. (2012). I keep
#' lambda fixed at 1. and the method is still experimental but looks promising.
#' CORNUET, J. M., MARIN, J. M., Mira, A., & Robert, C. P. (2012). Adaptive multiple importance sampling.
#' Scandinavian Journal of Statistics, 39(4), 798-812.
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom tmvtnorm dtmvnorm
#'
#' @keywords internal
#'
.gbmn_amis_easy <- function(object, preds, loss, new_data,
                            first_mad, obs_nodes, obs_mat,
                            niter, lambda, nsim){

  #fix preds and loss as iter
  pred_iter <- preds
  loss_iter <- loss

  #set constant
  eps <- .Machine$double.eps

  #create list for objects to store
  amis      <- vector("list", (niter+1))

  #set prior as q0
  amis[[1]] <- list(proposal=object, preds=preds, loss=loss)

  for(i in seq_len(niter)){

    #obtain weights
    weights            <- .gbmn_amis_reweight(object=object,
                                              amis=amis[seq_len(i)],
                                              lambda=lambda)

    #produce next proposal
    proposal_particles <- .gbmn_proposal_particles(object,
                                                   preds=weights$preds,
                                                   weights=weights$weights,
                                                   obs_nodes=obs_nodes)

    #predict using the new proposal
    pred_iter          <- predict_gbmn(proposal_particles$proposal,
                                       new_data = new_data, nsim = nsim)

    #calculate losses
    loss_iter          <- .gbmn_loss_resid(obs_mat=obs_mat, preds=pred_iter,
                                           obs_nodes=obs_nodes, nsim=nsim,
                                           MAD_used=first_mad)$loss1

    #sumarized per s simulation
    loss_iter          <- apply(loss_iter, 3, function(z){
      z <- z[is.finite(z)]; if(length(z) == 0){NA}else mean(z)})

    #store props, preds and losses
    amis[[i+1]]        <- list(proposal=proposal_particles$proposal,
                               preds=pred_iter,
                               loss=loss_iter)}

  #select parameter names e.
  ids    <- names(object$Parameters)
  ef_par <- rownames(object$Sigma)

  #geth all thetas in matrix format
  thetas <- do.call(rbind, lapply(amis, function(x){

    par_list  <- lapply(ids, function(id){
      draws <- x$preds$ParametersPerFunction[[id]]

      if(is.list(draws)){
        draws           <- do.call(rbind, draws)}
      draws           <- as.matrix(draws)
      colnames(draws) <- paste0(id, "_", names(x$proposal$Parameters[[id]]))

      draws})

    params <- do.call(cbind, par_list)
    params[,ef_par, drop=F]}))

  #all functions
  all_funcs <- do.call(rbind, lapply(amis, function(x){x$preds$FunctionSelection}))

  #all losses
  loss_all <- unlist(lapply(amis, function(x){x$loss}), use.names = FALSE)

  #all particles per proposal
  N_t <- vapply(amis, function(x){length(x$loss)}, integer(1))

  #extract function to get proper edges
  get_field <- function(prop, field){

    x <- unlist(lapply(prop$Parameters, function(edge){
      vapply(edge, function(p){as.numeric(p[field])}, numeric(1))}), use.names = F)

    setNames(x, ef_par)}

  theta_prob <- lapply(amis, function(x){

    prop         <- x$proposal

    mu_global    <- get_field(prop, "mu")
    lower_global <- get_field(prop, "a")
    upper_global <- get_field(prop, "b")

    amis_sigma <- prop$Sigma[ef_par, ef_par, drop=F]

    # Continuous log-density
    if(all(is.infinite(lower_global)) && all(is.infinite(upper_global))){
      log_q <- mvtnorm::dmvnorm(thetas, mean=mu_global, sigma=amis_sigma, log=T)

    }else{
      log_q <- tmvtnorm::dtmvnorm(thetas, mean=mu_global, sigma=amis_sigma,
                                  lower=lower_global, upper=upper_global, log=T)}

    #select probabilities per function
    for(node in colnames(all_funcs)){

      #id name per node
      ids_node <- prop$EdgeFunctions$id[prop$EdgeFunctions$child == node]

      if(length(ids_node) == 0){next}

      alpha <- prop$FunctionWeights[[node]][ids_node]

      #of 0 then set equal probs
      if(length(alpha) == 0 ||any(!is.finite(alpha)) ||sum(alpha) <= 0){

        alpha <- rep(1/length(ids_node), length(ids_node))}

      alpha <- alpha/sum(alpha)

      index <- match(all_funcs[,node], ids_node)

      log_q <- log_q+log(alpha[index])}

    log_q})

  #probabilities per simualtion
  prob_t <- do.call(cbind, theta_prob)

  #proposal draws correction for n
  log_terms <- sweep(prob_t, 2, log(N_t/sum(N_t)),"+")

  #get rowmeans and weights
  row_max <- apply(log_terms, 1, max)
  log_m   <- row_max+log(rowSums(exp(sweep(log_terms, 1, row_max, "-"))))

  #get the log weights
  log_w <- prob_t[, 1] -lambda*loss_all -log_m

  #set to neg -inf because on log scale
  log_w[!is.finite(log_w)] <- -Inf

  #set max log weight as minimal boundary
  weights <- exp(log_w-max(log_w))

  #normalize weights
  weights <- weights/sum(weights)

  #weighted parameter means
  mu_vec <- colSums(thetas*weights)

  #weighted covariance
  diff_mu_par <- sweep(thetas, 2, mu_vec, "-")

  #weighted covariance
  sigma_weighted <- t(diff_mu_par) %*% diag(weights) %*% diff_mu_par

  #add small constant to prevent issues
  sigma_weighted <- sigma_weighted+diag(eps, ncol(sigma_weighted))

  ############################
  #return all draws and preds#
  ############################

  #all parameter draws
  all_params <- setNames(lapply(ids, function(id){
    par_names       <- names(object$Parameters[[id]])
    draws           <- thetas[, paste0(id, "_", par_names),drop = FALSE]
    colnames(draws) <- par_names
    draws}), ids)

  #combined prediction structure
  preds_all <- list(FunctionSelection = all_funcs,
                    ParametersPerFunction = all_params)

  return(list(
    preds   = preds_all,
    weights = weights,
    info    = list(ESS=1/sum(weights^2), lambda=lambda, nsim=length(weights), niter=niter)))}
