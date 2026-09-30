#' Predict function for a PPMN object
#'
#' @param object A foundational or updated PPMN.
#' @param new_data New data for prediction.
#' @param nsim Number of simulations.
#' @param seed Seed value 123.
#'
#' @export
predict_ppmn <- function(object, new_data, nsim = 1000, seed = 123){

  if(nsim<1){stop("Number of simulations cannot be smaller than 1.")}
  if(!inherits(object, "ppmn")){class(object) <- c("ppmn", "list")}

  #set seed
  set.seed(seed)

  #select nodes
  node_names  <- igraph::V(object$Structure$visual_dag$graph)$name
  given_nodes <- intersect(names(new_data), node_names)

  #store only roots as starting values
  input_values <- as.list(new_data)
  input_values <- input_values[intersect(object$Roots, names(input_values))]

  ###################################
  #edge function and child structure#
  ###################################

  edge_info <- object$EdgeFunctions
  if(is.null(edge_info)||nrow(edge_info)==0){stop("No edge function information present in object.")}

  child_nodes <- unique(edge_info$child)
  n_child     <- length(child_nodes)
  max_f       <- max(edge_info$f)

  #parents per edge function
  if(!is.null(object$Structure$predict_dag$parents)){
    parent_list <- object$Structure$predict_dag$parents
  }else if("parents" %in% names(edge_info)){
    parent_list <- strsplit(edge_info$parents, "\\+")
    names(parent_list) <- edge_info$id
  }else{stop("No parent information present in object.")}

  #############################
  #global parameter extraction#
  #############################

  mu_global    <- unlist(lapply(object$Parameters, function(edge){sapply(edge, function(p) p["mu"])}))
  lower_global <- unlist(lapply(object$Parameters, function(edge){sapply(edge, function(p) p["a"])}))
  upper_global <- unlist(lapply(object$Parameters, function(edge){sapply(edge, function(p) p["b"])}))

  #set names
  if(!is.null(rownames(object$Sigma))){
    names(mu_global)    <- rownames(object$Sigma)
    names(lower_global) <- rownames(object$Sigma)
    names(upper_global) <- rownames(object$Sigma)}

  ################
  #MC simulations#
  ################

  if(all(is.infinite(lower_global)) && all(is.infinite(upper_global))){
    theta_draws <- MASS::mvrnorm(n=nsim, mu=mu_global, Sigma=object$Sigma)
  }else{
    theta_draws <- tmvtnorm::rtmvnorm(n=nsim, mean=mu_global, sigma=object$Sigma,
                                      lower=lower_global, upper=upper_global)}

  theta_draws   <- matrix(theta_draws, nrow = nsim, dimnames = list(seq_len(nsim), names(mu_global)))

  ##################################
  #store predictions and parameters#
  ##################################

  predictions <- vector("list", nsim)
  for(s in seq_len(nsim)){predictions[[s]] <- input_values}

  Variance <- list()
  Expected <- input_values

  #store predictions per edge function [i,j,s,f] i=row (observatino) j=col (vertex), s=simulation, f=function
  pred_fun <- array(NA_real_, dim = c(nrow(new_data), n_child, nsim, max_f),
                    dimnames = list(NULL, child_nodes, NULL, seq_len(max_f)))

  #this stores the selected functions f per simulations s
  FunctionSelection <- matrix(NA_character_, nrow = nsim, ncol = n_child,
                              dimnames = list(seq_len(nsim), child_nodes))

  #store expected values per edge function [i,j,f]
  exp_fun  <- array(NA_real_, dim = c(nrow(new_data), n_child, max_f),
                    dimnames = list(NULL, child_nodes, seq_len(max_f)))

  #store parameter draws per edge function
  Parsbyfun                   <- vector("list", nrow(edge_info))
  names(Parsbyfun)            <- edge_info$id

  #for each edge function store the parameters in a list
  for(k in seq_len(nrow(edge_info))){

    #subset edge function
    id  <- edge_info$id[k]

    #select columns with edge fun id
    sel <- grepl(paste0("^",id,"_"), colnames(theta_draws))

    #for simulation s select the column with the particles
    #and store it in the parsbyfun (parameters by function)
    Parsbyfun[[id]] <- lapply(seq_len(nsim), function(s){
      as.numeric(theta_draws[s, sel, drop=F])})}

  ######################################
  #weighted average over edge functions#
  ######################################

  avg_fun <- function(mat, w){

    out <- rep(NA_real_, nrow(mat))
    for(i in seq_len(nrow(mat))){

      ok <- is.finite(mat[i,])
      if(any(ok)){
        out[i] <- sum(mat[i,ok]*w[ok])/sum(w[ok])}}

    out}

  ##########
  #DAG loop#
  ##########

  topo_order <- object$Structure$exec_dag$order

  #pain in the ... to code and figure out how to correctly
  #predict down the chains
  for(node in topo_order){

    #select edge functions belonging to child node
    fun_select <- edge_info$id[edge_info$child==node]
    if(length(fun_select)==0){next}

    #which parent node fit
    j <- match(node, child_nodes)

    #function weights for averaging
    w_fun <- object$FunctionWeights[[node]][fun_select]
    if(length(w_fun)==0||any(!is.finite(w_fun))||sum(w_fun)<=0){
      w_fun <- rep(1/length(fun_select), length(fun_select))
    }else{
      w_fun <- as.numeric(w_fun/sum(w_fun))}

    #create a vector of selecte functions based on the function weight
    selected_fun <- sample(fun_select, size = nsim, replace = TRUE, prob = w_fun)

    #stor the vector per node
    FunctionSelection[,node] <- selected_fun

    ############################
    #expected function response#
    ############################

    for(f in seq_along(fun_select)){

      id      <- fun_select[f]
      parents <- parent_list[[id]]

      x <- do.call(cbind, lapply(parents, function(p){Expected[[p]]}))
      x <- as.matrix(x)

      beta <- mu_global[grepl(paste0("^",id,"_"), names(mu_global))]
      beta <- as.numeric(beta)

      yhat <- .ppmn_predict_edge(edge_info$fun[edge_info$id==id], beta,
                           edge_info$trans[edge_info$id==id], x)

      exp_fun[,j,f]  <- yhat}

    mat_exp          <- exp_fun[,j,seq_along(fun_select),drop=F]
    dim(mat_exp)     <- c(nrow(new_data), length(fun_select))
    Expected[[node]] <- avg_fun(mat_exp, w_fun)

    ##########################
    #stochastic mc prediction#
    ##########################

    Variance[[node]] <- vector("list", nsim)

    for(s in seq_len(nsim)){

      for(f in seq_along(fun_select)){

        #select fun and id for parents per fun
        id      <- fun_select[f]
        parents <- parent_list[[id]]

        #input values (parents)
        x <- do.call(cbind, lapply(parents, function(p){predictions[[s]][[p]]}))
        x <- as.matrix(x)

        #parameters per edge and simulation
        beta <- Parsbyfun[[id]][[s]]

        #predict with pred edge
        yhat <- .ppmn_predict_edge(edge_info$fun[edge_info$id==id], beta,
                             edge_info$trans[edge_info$id==id], x)

        #not really fun though
        pred_fun[,j,s,f] <- yhat}

      mat_pred      <- pred_fun[,j,s,seq_along(fun_select),drop=F]
      dim(mat_pred) <- c(nrow(new_data), length(fun_select))

      m_s <- match(FunctionSelection[s, node], fun_select)

      predictions[[s]][[node]] <- mat_pred[, m_s]
      Variance[[node]][[s]] <- predictions[[s]][[node]]}
  }

  ###########################
  #old style combined output#
  ###########################

  Parameters           <- list()
  for(node in child_nodes){
    fun_select         <- edge_info$id[edge_info$child==node]
    Parameters[[node]] <- lapply(seq_len(nsim), function(s){

      mats <- lapply(fun_select, function(id){Parsbyfun[[id]][[s]]})

      names(mats) <- fun_select
      mats})}

  out <- list(
    Expected = Expected,
    Variance = Variance,
    Parameters = Parameters,
    ParametersPerFunction = Parsbyfun,
    FunctionPredictions = pred_fun,
    ExpectedFunctions = exp_fun,
    FunctionSelection = FunctionSelection,
    FunctionWeights = object$FunctionWeights,
    EdgeFunctions = edge_info,
    GivenNodes = given_nodes)

  class(out) <- c("ppmn_prediction", "list")

  return(out)}
