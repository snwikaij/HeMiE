#' Bootstrap coverage of generalized posterior expected predictions
#'
#' @param loss_full Full loss array (N x J x S).
#' @param preds Predictions
#' @param obs_nodes Observed children
#' @param obs_mat Observed new data
#' @param max_lambda Maximum value for lambda.
#' @param level Credibility level.
#' @param nboot Number of bootstraps
#' @param nsim Number of simulations
#' @param shrink number of contractions around lambda
#'
#' @description
#' Very experimental it concentrates the generalized posterior too fast. It checks whether the
#' average of all expectations is covered x% by the credibility under bootstrapping.
#'
#' @keywords internal
#'
.gbmn_coverage_fast <- function(loss_full, preds, obs_nodes, obs_mat, max_lambda, level, nboot, nsim, shrink=3){

  #which row to bootstrap
  selection_boot <- lapply(1:nboot, function(i) sample(seq_len(dim(loss_full)[1]), dim(loss_full)[1], T))

  #loss per bootstrap
  loss_b <- t(sapply(selection_boot, function(i){
    apply(loss_full[i,,,drop=F], 3, function(z){z <- z[is.finite(z)]
    if(length(z)==0){NA}else{mean(z)}})}))

  #set a loss per boot as matrix
  loss_b    <- as.matrix(loss_b)

  #set loss per simulation/particle
  loss_s <- apply(loss_full, 3, function(z){
    z <- z[is.finite(z)]; if(length(z)==0){NA}else{mean(z)}})

  #all predictions
  pred_s <- simplify2array(lapply(obs_nodes, function(node){
    do.call(cbind, preds$Variance[[node]])}))

  #what is observed
  observed <- is.finite(obs_mat[, obs_nodes, drop=FALSE])

  #coverage prob is
  notcover <- (1-level)/2

  #create sequence
  seq_lambda <- seq(0, max_lambda, length.out=20)

  for(i in 1:shrink){

   coverage <- vapply(seq_lambda, function(learn_rate){

    #calculate weights for lambda in sequence a, ..., A for
    #the full model and data
    w_total <- .gbmn_stable_weights(loss_s, learn_rate)

    #and bootstrapped data
    w_boot  <- t(vapply(seq_len(nboot), function(b){
                 .gbmn_stable_weights(loss_b[b,], learn_rate)}, numeric(nsim)))

    #calculate expected value for ij
    xhat_ij <- sapply(seq_along(obs_nodes), function(j){rowSums(sweep(pred_s[,,j], 2, w_total, "*"))})

    #predictions lower and higher than xhat ij
    lower_pred  <- sweep(pred_s, c(1,3), xhat_ij, "<=")
    higher_pred <- sweep(pred_s, c(1,3), xhat_ij, "<")

    cover <- sapply(seq_along(obs_nodes), function(j){

      ix <- which(observed[,j] & is.finite(xhat_ij[,j]))
      if(length(ix)==0){return(NA)}

      #true/false matrices
      le <- matrix(lower_pred[ix,,j,drop=FALSE], ncol=nsim)
      lt <- matrix(higher_pred[ix,,j,drop=FALSE], ncol=nsim)

      #cumulative bootstrap weights
      cdf_lower  <- w_boot %*% t(le*1)
      cdf_higher <- w_boot %*% t(lt*1)

      #coverage
      mean(cdf_lower>=notcover & cdf_higher<notcover+level, na.rm=TRUE)})

    #weight coverage according to finite predictions per vertex
    cover_weight <- colSums(observed & is.finite(xhat_ij))
    cover        <- weighted.mean(cover, cover_weight, na.rm=TRUE)

   }, numeric(1))

   #create lambda coverage df
   df_cover   <- data.frame(lambda=seq_lambda, coverage=coverage)

   #calculated difference
   diff_level <- df_cover$coverage-level

   #select index x
   ix         <- (which(head(diff_level, -1) >= 0 & tail(diff_level, -1) <= 0))

   if(length(ix)>0){

     #select lower and higher learning rates than level
     lambda_bound <- df_cover$lambda[c(ix[1], ix[1]+1)]

     #Search within those bounds
     seq_lambda <- seq(min(lambda_bound), max(lambda_bound), length.out=20)}

   }

  df_cover$lambda[which.min(abs(df_cover$coverage-level))]}
