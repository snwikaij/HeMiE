#' Helper function to generate proposals
#'
#' @param object A foundational or updated GBMN.
#' @param preds Predictions from the GBMN
#' @param weights Weights based on the loss
#' @param obs_nodes Observed nodes/vertices (I am so inconsistent with nodes/vertices).
#'
#' @keywords internal
#'
.gbmn_proposal_particles <- function(object, preds, weights, obs_nodes){

  #constant
  eps              <- .Machine$double.eps

  #extract edge information
  edge_info        <- object$EdgeFunctions

  #extract weights multi-edge
  function_weights <- object$FunctionWeights

###################
#observed vertices#
###################

  #extract graph
  g              <- object$Structure$exec_dag$graph

  #select nodes in the graph, I am so inconsistent with nodes and vertex
  relevant_nodes <- unique(unlist(lapply(obs_nodes, function(node){
    igraph::as_ids(igraph::subcomponent(g, node, mode = "in"))})))

  #select the children
  relevant_nodes <- intersect(relevant_nodes, unique(edge_info$child))

  #select of the observed vertices select whou are present
  relevant_nodes <- intersect(relevant_nodes, colnames(preds$FunctionSelection))

#############################################
#summarize multi-edge function probabilities#
#############################################

  for(node in relevant_nodes){

    ids      <- edge_info$id[edge_info$child == node]
    selected <- preds$FunctionSelection[, node]

    alpha_new <- vapply(ids, function(id){
      active <- !is.na(selected) & selected == id
      sum(weights[active])}, numeric(1))

    names(alpha_new) <- ids

    if(sum(alpha_new) > 0){

      alpha_new <- alpha_new/sum(alpha_new)

    }else{

      alpha_new <- function_weights[[node]][ids]

      if(length(alpha_new) != length(ids) || any(!is.finite(alpha_new)) || sum(alpha_new) <= 0){
        alpha_new <- rep(1/length(ids), length(ids))}

      alpha_new        <- alpha_new/sum(alpha_new)
      names(alpha_new) <- ids}

    function_weights[[node]] <- alpha_new}

#############################
#extract params per function#
#############################

  function_draws     <- list()

  #for all edge functions select the parameters
  for(id in names(object$Parameters)){

    #all parameter draws
    draws <- preds$ParametersPerFunction[[id]]

    #bind to matrix
    if(is.list(draws)){draws <- do.call(rbind, draws)}

    #set as matrix preventing ambiguity (ugh)
    draws <- as.matrix(draws)

    #set as b0, b1, b2, ..., bn
    colnames(draws)      <- names(object$Parameters[[id]])
    function_draws[[id]] <- draws}

#########
#results#
#########

  #select parameters
  params <- as.matrix(do.call(cbind, function_draws))

  theta_names <- unlist(Map(
    function(id, pars) paste0(id, "_", names(pars)),
    names(object$Parameters), object$Parameters), use.names = FALSE)

  colnames(params)         <- theta_names

  #based on the weights make new proposal
  mu_vec                   <- colSums(params*weights)

  #substract mean
  diff_mu_par              <- sweep(params, 2, mu_vec, "-")

  #standard error
  se_vec                   <- sqrt(colSums((diff_mu_par^2)*weights))

  #weighted covariance
  sigma_weighted           <- crossprod(diff_mu_par, diff_mu_par*weights)

  #add small constant to prevent issues
  sigma_weighted           <- sigma_weighted+diag(eps, ncol(sigma_weighted))

  #add names
  dimnames(sigma_weighted) <- list(theta_names, theta_names)

############################
#Generate a proposal object#
############################

  #proposal
  proposal_object <- object

  #update parameter means and se
  count_dracula   <- 0

  for(id in names(proposal_object$Parameters)){

    for(b in names(proposal_object$Parameters[[id]])){

      count_dracula <- count_dracula+1

      proposal_object$Parameters[[id]][[b]]["mu"] <- mu_vec[count_dracula]
      proposal_object$Parameters[[id]][[b]]["se"] <- se_vec[count_dracula]}}

  #update covariance
  proposal_object$Sigma <- sigma_weighted

  #fun weights probs
  proposal_object$FunctionWeights <- function_weights

  return(list(proposal=proposal_object,
              function_weights=function_weights,
              mu=mu_vec,
              sigma=sigma_weighted))}
