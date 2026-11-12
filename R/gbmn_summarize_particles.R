#' Helper function to summarize the output of different calibration methods
#'
#' @param object An foundational or updated GBMN
#' @param preds Predictions of the GBMN
#' @param weights Weights of the GBMN
#' @param obs_nodes Observed nodes (names)
#' @param level credibility level
#' @param covar return covariance.
#'
#' @keywords internal
#'
.gbmn_summarize_particles <- function(object, preds, weights, obs_nodes, level, covar){

  #covar offset constant
  eps        <- 1e-8

  #normalize weights
  weights   <- weights/sum(weights)

  #edge info
  edge_info <- object$EdgeFunctions

  #function weights
  function_weights          <- object$FunctionWeights
  function_particle_weights <- list()

  #select graph
  g <- object$Structure$exec_dag$graph

  #select relevant nodes
  relevant_nodes <- unique(unlist(lapply(obs_nodes, function(node){
    igraph::as_ids(igraph::subcomponent(g,node,mode = "in"))})))

  #retain childs that have edge functions
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

    function_weights[[node]] <- alpha_new

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

    if(is.list(draws)){
      draws <- do.call(rbind, draws)}
      draws <- as.matrix(draws)

    #use only particles in which function was active
    w               <- function_particle_weights[[id]]

    #if this function was never selected, there is no
    #information with which to update its parameters
    if(is.null(w) || length(w)==0 || sum(w, na.rm=TRUE)<=0){

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

    #weighted covariance
    sigma_weighted <- t(diff_mu_par) %*% diag(weights) %*% diff_mu_par

    #add small constant to prevent issues
    sigma_weighted <- sigma_weighted+diag(eps, ncol(sigma_weighted))}

  #give names to covar mat
  theta_names <- unlist(Map(function(id, pars)paste0(id,"_",names(pars)),
                            names(object$Parameters), object$Parameters))
  #set names
  dimnames(sigma_weighted) <- list(theta_names, theta_names)

return(list(function_weights  = function_weights,
    function_particle_weights = function_particle_weights,
    function_summaries        = function_summaries,
    sigma_weighted            = sigma_weighted))}
