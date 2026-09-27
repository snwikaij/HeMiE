#' Plots the marginal expectations of each edge function with credibility intervals
#'
#' @param object A foundational or updated PPMN which will be used for predictions.
#' @param child The dependent (child) variable to display.
#' @param data Data that contains x and y coordinates (points) to plot.
#' @param level Credibility level (default level = 0.9).
#' @param nsim Number of simulations to construct the intervals.
#' @param n_grid Grid size of the added variable plot.
#' @param smooth_f A value for the smoother span only used to display specific non-linear models.
#' @param ylim Limit of the y-axis.
#' @param ylab Title of the y-axis.
#' @param lwd_exp Line width of the expected value.
#' @param lwd_col Colour of the expected value.
#' @param size_title Title size of the y-axis label.
#' @param size_text Text size of the title.
#' @param int_col Colour of the credibility intervals around the expected value.
#' @param int_alpha Transparency of the credibility intervals around the expected value.
#' @param pt_size Size of the points in the plot.
#' @param pt_alpha Transparency of the points in the plot.
#' @param seed Seed value 123.
#'
#' @description
#' This is essentially an added variable plot. It displays the marginal relation of the
#' dependent variable and all independent variables for each edge-function belonging to
#' the selected child. If multiple edge-functions predict the same child, a separate
#' marginal plot is returned for every edge-function and parent combination.
#'
#' @importFrom stats lowess approx median quantile
#' @importFrom ggplot2 ggplot geom_line geom_ribbon xlab element_text geom_point aes aes_string ggplot_build
#'
#' @export
marginal_ppmn <- function(object, child, data,
                          level=0.9, nsim=1000, n_grid = 100, smooth_f=0.5,
                          ylim=NULL, ylab=NULL, lwd_exp = 0.8, lwd_col="dodgerblue4", size_title=8, size_text=8,
                          int_col="dodgerblue3", int_alpha=0.3, pt_size=2, pt_alpha=0.25, seed=123){

  #sets seed
  set.seed(seed)

  #expected graphical structure
  g <- object$Structure$exec_dag$graph

#################
#draw parameters#
#################

  #select means
  mu_global        <- unlist(lapply(object$Parameters, function(e) sapply(e, `[[`, "mu")))
  names(mu_global) <- rownames(object$Sigma)

  #selectu upper lower truncation
  lower <- unlist(lapply(object$Parameters, function(e) sapply(e, `[[`, "a")))
  upper <- unlist(lapply(object$Parameters, function(e) sapply(e, `[[`, "b")))

  #if yes truncated or if not trucated
  if(all(is.infinite(lower)) && all(is.infinite(upper))){
    theta <- MASS::mvrnorm(nsim, mu_global, object$Sigma)
  }else{
    theta <- tmvtnorm::rtmvnorm(nsim, mu_global, object$Sigma, lower, upper)}

  #store in matrix
  theta <- matrix(theta, nrow = nsim,
                  ncol = length(mu_global),
                  dimnames = list(seq_len(nsim), names(mu_global)))

#####################################
#select edge functions for the child#
#####################################

  #select edge info and stor it
  edge_info <- object$EdgeFunctions
  edge_info <- edge_info[edge_info$child==child,,drop=F]

  if(nrow(edge_info)==0){stop("No edge function found for child: ",child)}

  #make predictions select expected outcome
  pred_exp <- predict_ppmn(object, data, nsim=1)$Expected

###################
#get parent values#
###################

#do I include this as helper function yes/no????
get_parent_values <- function(p, data, pred_exp, trans = NA){

  #if not a parent get data
  if(!p %in% colnames(data)||!is.numeric(data[[p]])){
      x <- pred_exp[[p]]
  }else{
      x <- data[[p]]}

  if(!is.na(trans) && trans=="log"){
      x[x<=0] <- NA
      x <- log(x)}

  if(!is.na(trans) && trans=="sqrt"){
      x[x<0] <- NA
      x <- sqrt(x)}

  if(!is.na(trans) && trans=="log10"){
      x[x<=0] <- NA
      x <- log10(x)}

      x[!is.finite(x)] <- NA
      x}

#########################
#loop over edge function#
#########################

  plots         <- list()
  Y_weighted    <- list()
  grid_weighted <- list()

  for(f in seq_len(nrow(edge_info))){

    #get all infor of edge
    id        <- edge_info$id[f]
    edge_name <- edge_info$edge[f]
    func      <- edge_info$fun[f]
    trans     <- edge_info$trans[f]

    #select parent of the edge fun
    parents <- object$Structure$predict_dag$parents[[id]]
    if(length(parents)==0){next}

    #edge fun should be parents of child
    graph_parents <- igraph::V(g)$name[igraph::neighbors(g, child, mode="in")]
    if(!all(parents %in% graph_parents)){stop("Parent mismatch for edge-function: ",id)}

    #select parameter for edge
    edge_pars <- object$Parameters[[id]]
    param_idx <- paste0(id,"_",names(edge_pars))
    param_idx <- param_idx[param_idx %in% colnames(theta)]

    if(length(param_idx)==0){stop("No parameters found for edge-function: ",id)}

#######################
#loop over each parent#
#######################

    for(i in seq_along(parents)){

      vary <- parents[i]

#############
#create grid#
#############

   #for log and log10 0 cannot be included in the grid
   if(!is.na(trans) && trans %in% c("log", "log10")){

      if(!vary %in% colnames(data)){
         xgrid <- pred_exp[[vary]]
      }else{
         xgrid <- data[[vary]]}
         xgrid <- xgrid[is.finite(xgrid) & xgrid>0]
      }else{
       if(!vary %in% colnames(data)){
         xgrid <- pred_exp[[vary]]
      }else{
         xgrid <- data[[vary]]}
         xgrid <- xgrid[is.finite(xgrid)]}

      if(length(xgrid)<2){next}

      grid <- seq(min(xgrid,na.rm=T), max(xgrid,na.rm=T), length.out=n_grid)

#######################################
#mvlogistic uses observed combinations#
#######################################

  #gives wonky curve when keeping one variable fixed at mean E(y)=f(x,mean(x),theta)
  #the median gives also a more stable output then mean
  if(func=="mvlogistic"){

      X_actual <- do.call(cbind, lapply(parents, function(p){

    if(p %in% colnames(data) && is.numeric(data[[p]])){
         x <- data[[p]]
    }else{
         x <- pred_exp[[p]]}

         x[!is.finite(x)] <- NA
         x}))

    if(vary %in% colnames(data) && is.numeric(data[[vary]])){
         x_vary <- data[[vary]]
    }else{
         x_vary <- pred_exp[[vary]]}

         Y <- matrix(NA, n_grid, nsim)

    for(s in seq_len(nsim)){
          beta <- theta[s,param_idx]
          y    <- .ppmn_predict_edge(func,beta,trans,X_actual)
          ok   <- is.finite(x_vary) & is.finite(y)

    if(sum(ok)>=3 && length(unique(x_vary[ok]))>=2){

          #this was chatGPTs suggestion I could not figure out the issue
          #with the curve it actually quite a nice solution.
          sm    <- lowess(x_vary[ok], y[ok], f=smooth_f)
          Y[,s] <- approx(x=sm$x, y=sm$y, xout=grid, rule=2, ties=median)$y}}}

##########################
#all other edge-functions#
##########################

  else{

    X <- do.call(cbind,lapply(parents,function(p){

    if(p==vary){
       grid
    }else{

    xc <- get_parent_values(p, data, pred_exp, trans)

    if(!is.na(trans) && trans=="log"){
       rep(exp(median(xc, na.rm=T)), n_grid)
    }else if(!is.na(trans) && trans=="log10"){
       rep(10^(median(xc, na.rm=T)), n_grid)
    }else if(!is.na(trans) && trans=="sqrt"){
       rep((median(xc, na.rm=T))^2, n_grid)
    }else{
       rep(median(xc, na.rm=T),n_grid)}}}))
        Y <- matrix(NA, n_grid, nsim)

     for(s in seq_len(nsim)){

          beta <- theta[s, param_idx]

          y <- .ppmn_predict_edge(func, beta, trans, X)
          y[!is.finite(y)] <- NA

          Y[,s] <- y}}

########################################
#store prediction of multiple functions#
########################################

if(is.null(Y_weighted[[vary]])){Y_weighted[[vary]] <- list()}

           Y_weighted[[vary]][[id]] <- Y
           grid_weighted[[vary]]    <- grid

#######################
#summarize simulations#
#######################

df <- data.frame(
      x  = grid,
      mu = apply(Y,1,median,na.rm=T),
      ll = apply(Y,1,quantile,probs=(1-level)/2,na.rm=T),
      ul = apply(Y,1,quantile,probs=1-(1-level)/2,na.rm=T))

    if(!is.null(ylim)){
       df$ll <- ifelse(ylim[1]>df$ll,ylim[1],df$ll)
       df$ul <- ifelse(ylim[2]<df$ul,ylim[2],df$ul)
       ylim_func <- ggplot2::ylim(ylim)
    }else{
       ylim_func <- NULL}

#######################
#observed marginal fit#
#######################

      if(!vary %in% colnames(data)){
        x_obs <- pred_exp[[vary]]
      }else{
        x_obs <- data[[vary]]}

      if(child %in% colnames(data)){
        y_obs <- data[[child]]
      }else{
        y_obs <- rep(NA,length(x_obs))}

      if(func=="mvlogistic"){

        X_obs <- do.call(cbind,lapply(parents,function(p){

          if(p %in% colnames(data) && is.numeric(data[[p]])){
            x <- data[[p]]
          }else{
            x <- pred_exp[[p]]}
            x[!is.finite(x)] <- NA
            x}))
      }else{
        X_obs <- do.call(cbind,lapply(parents,function(p){
          if(p==vary){
            x_obs
          }else{

            xc <- get_parent_values(p,data,pred_exp,trans)
            n  <- length(x_obs)

            if(!is.na(trans) && trans=="log"){
              rep(exp(median(xc, na.rm=T)), n)
            }else if(!is.na(trans) && trans=="log10"){
              rep(10^(median(xc, na.rm=T)), n)
            }else if(!is.na(trans) && trans=="sqrt"){
              rep((median(xc, na.rm=T))^2, n)
            }else{
              rep(median(xc, na.rm=T),n)}}}))}

    #posterior predictive expectation at observed x
    y_hat <- rowMeans(
      sapply(seq_len(nsim),function(s){

      beta <- theta[s,param_idx]

      yh <- .ppmn_predict_edge(func, beta, trans, X_obs)
      yh[!is.finite(yh)] <- NA

        yh}),
      na.rm=T)

###############
#create figure#
###############

    if(!is.null(ylab)){ylab <- ylab}else{ylab <- child}

    p <- ggplot(df,aes(x,mu))+
      geom_line(lwd=lwd_exp,col=lwd_col)+ylim_func+
      geom_ribbon(aes(ymin=ll,ymax=ul),fill=int_col,alpha=int_alpha)+
      xlab(vary)+ylab(ylab)+
      theme_classic()+
      theme(axis.title=element_text(size=size_title),
            axis.text=element_text(size=size_text))+
      if(!vary %in% colnames(data) || !child %in% colnames(data)){NULL}else{
        geom_point(data=data,
                   aes_string(x=vary,y=child),
                   inherit.aes=F,pch=19,alpha=pt_alpha,size=pt_size)}
      ggplot_build(p)

      plots[[paste0(id,"_",vary)]] <- p}}

#############################
#weighted edge function plot#
#############################

  if(nrow(edge_info)>1){

    for(vary in names(Y_weighted)){

     Y_list <- Y_weighted[[vary]]
     ids    <- names(Y_list)

     #function weights
     w_fun <- object$FunctionWeights[[child]][ids]

     #if not weights then likely 1 or 0.5/0.5 at least uniform
     if(length(w_fun)==0||any(!is.finite(w_fun))||sum(w_fun)<=0){
        w_fun <- rep(1/length(ids),length(ids))
     }else{
        w_fun <- w_fun/sum(w_fun)}

      #select the grids
      grid <- grid_weighted[[vary]]

      #weighted prediction per simulation
      Y <- matrix(NA, nrow=n_grid, ncol=nsim)

      for(s in seq_len(nsim)){

        for(i in seq_len(n_grid)){

          y <- sapply(ids,function(id){
          Y_list[[id]][i,s]})

          ok <- is.finite(y)

          if(any(ok)){
            Y[i,s] <- sum(y[ok]*w_fun[ok])/sum(w_fun[ok])}}}

######################
#collapse simulations#
######################

#summarize eahs simulation at each grid point
df <- data.frame(
      x  = grid,
      mu = apply(Y,1, median,na.rm=T),
      ll = apply(Y,1, quantile, probs=(1-level)/2, na.rm=T),
      ul = apply(Y,1, quantile, probs=1-(1-level)/2, na.rm=T))

    if(!is.null(ylim)){
      df$ll <- ifelse(ylim[1]>df$ll,ylim[1],df$ll)
      df$ul <- ifelse(ylim[2]<df$ul,ylim[2],df$ul)
      ylim_func <- ggplot2::ylim(ylim)
    }else{
      ylim_func <- NULL}

    if(!is.null(ylab)){
      ylab_use <- ylab
    }else{
      ylab_use <- child}

    p <- ggplot(df,aes(x,mu)) +
        geom_line(lwd=lwd_exp,col=lwd_col) +
        ylim_func +
        geom_ribbon(aes(ymin=ll,ymax=ul), fill=int_col, alpha=int_alpha) +
        xlab(vary) + ylab(ylab_use)+
        theme_classic()+
        theme(axis.title=element_text(size=size_title),
              axis.text=element_text(size=size_text))+
        if(!vary %in% colnames(data)||!child %in% colnames(data)){NULL}else{
          geom_point(data=data, aes_string(x=vary,y=child),
            inherit.aes=F, pch=19, alpha=pt_alpha,
            size=pt_size)}

      ggplot_build(p)

      plots[[paste0("weighted_",vary)]] <- p}}

  return(plots)}

