#' Helper function to reweight AMIS particles
#'
#' @param object A foundational or updated gbmn object.
#' @param amis List containing proposals, predictions and losses.
#' @param lambda Learning rate of the GBU.
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom tmvtnorm dtmvnorm
#'
#' @keywords internal
#'
.gbmn_amis_reweight <- function(object, amis, lambda){

  #select parameter names
  ids    <- names(object$Parameters)
  ef_par <- rownames(object$Sigma)

  #get thetas into matrix
  thetas <- do.call(rbind, lapply(amis, function(x){

    par_list <- lapply(ids, function(id){

      draws <- x$preds$ParametersPerFunction[[id]]

      if(is.list(draws)){
        draws <- do.call(rbind, draws)}

      draws           <- as.matrix(draws)
      colnames(draws) <- paste0(id, "_", names(x$proposal$Parameters[[id]]))

      draws})

    params <- do.call(cbind, par_list)
    params[,ef_par,drop=F]}))

  #all functions
  all_funcs <- do.call(rbind, lapply(amis, function(x){
    as.matrix(x$preds$FunctionSelection)}))

  #all losses
  loss_all <- unlist(lapply(amis, function(x){
    x$loss}), use.names=F)

  #all particles per proposal
  N_t <- vapply(amis, function(x){length(x$loss)}, integer(1))

#############################
#extract parameter densities#
#############################

  #extract function to get proper edges
  get_field <- function(prop, field){

    x <- unlist(lapply(names(prop$Parameters), function(id){

      edge <- prop$Parameters[[id]]

      values <- vapply(edge, function(p){as.numeric(p[field])}, numeric(1))

      names(values) <- paste0(id, "_", names(edge))

      values}), use.names=T)

    as.numeric(x[ef_par])}

  #calculate densities of all particles under all proposals
  theta_prob <- lapply(amis, function(x){

    #new proposal mu and se
    prop         <- x$proposal

    #mu and ll and ul truncs
    mu_global    <- get_field(prop, "mu")
    lower_global <- get_field(prop, "a")
    upper_global <- get_field(prop, "b")

    #get sigma
    amis_sigma <- prop$Sigma[ef_par,ef_par,drop=F]

    #continuous log-density
    if(all(is.infinite(lower_global) & lower_global<0) && all(is.infinite(upper_global) & upper_global>0)){

      #mvn normal
      log_q <- mvtnorm::dmvnorm(thetas, mean=mu_global, sigma=amis_sigma, log=T)

    }else{

      #use trunc mvn normal
      log_q <- tmvtnorm::dtmvnorm(thetas, mean=mu_global, sigma=amis_sigma, lower=lower_global, upper=upper_global, log=T)}

    #select probabilities per function
    for(node in colnames(all_funcs)){

      #id name per node
      ids_node <- prop$EdgeFunctions$id[prop$EdgeFunctions$child == node]

      if(length(ids_node)==0){next}

      alpha <- prop$FunctionWeights[[node]][ids_node]

      #if 0 then set equal probabilities
      if(length(alpha)==0 || any(!is.finite(alpha)) ||
         sum(alpha)<=0){

        alpha <- rep(1/length(ids_node),length(ids_node))}

      alpha <- alpha/sum(alpha)

      index <- match(all_funcs[,node],ids_node)

      #only add density when a function was selected
      if(any(is.na(index) & !is.na(all_funcs[,node]))){
        stop("Function selection does not match proposal.")}

      active <- !is.na(index)

      log_q[active] <- log_q[active]+log(alpha[index[active]])}

    log_q})

##############
#AMIS weights#
##############

  #probs per simulation
  prob_t    <- do.call(cbind, theta_prob)

  #proposal draws correction for inequality in integers
  log_terms <- sweep(prob_t, 2, log(N_t/sum(N_t)),"+")

  #get rowmeans and weights
  row_max <- apply(log_terms, 1, max)

  #calculate log mixture density
  log_m <- rep(-Inf, length(row_max))

  #select finite values
  valid <- is.finite(row_max)

  #calculate stable log mean
  if(any(valid)){

    log_m[valid] <- row_max[valid]+
      log(rowSums(exp(sweep(
        log_terms[valid,,drop=F],
        1,row_max[valid],"-"))))}

  #get the log weights re-apply gbu
  log_w <- prob_t[,1] -lambda*loss_all -log_m

  #set to neg -inf because on log scale
  log_w[!is.finite(log_w)] <- -Inf

  #check weights
  if(!any(is.finite(log_w))){
    stop("No finite AMIS weights.")}

  #set max log weight as minimal boundary
  weights <- exp(log_w-max(log_w))

  #normalize weights
  weights <- weights/sum(weights)

  ############################
  #return all draws and preds#
  ############################

  #all parameter draws
  all_params <- setNames(lapply(ids, function(id){

    par_names <- names(object$Parameters[[id]])

    draws <- thetas[,paste0(id,"_",par_names),drop=F]
    colnames(draws) <- par_names

    draws}),ids)

  #combined prediction structure
  preds_all <- list(
    FunctionSelection=all_funcs,
    ParametersPerFunction=all_params)

  return(list(
    preds   = preds_all,
    weights = weights,
    ESS     = 1/sum(weights^2)))}

