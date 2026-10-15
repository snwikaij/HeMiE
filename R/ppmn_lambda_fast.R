#' Helper function oob fast
#'
#' @param loss_array An array fo losses.
#' @param loss Compressed loss per s.
#' @param max_lambda Maximum lambda.
#' @param n_split Number of repetitions.
#' @param n_lambda Number of steps for lambda.
#' @param ess_frac Effective sample size boundary.
#'
#' @keywords internal
#'
.ppmn_lambda_fast <- function(loss_array,  max_lambda,
                              n_rep, n_lambda=200){
  #number of samples
  n_samp <- dim(loss_array)[1]

  #total loss per [s]
  loss_s <- apply(loss_array, 3, function(x){
    x <- x[is.finite(x)]; if(length(x)==0){NA_real_}else{mean(x)}})

  #get the ESS to prevent collapse
  get_ess <- function(l){
    w <- .ppmn_stable_weights(loss_s, l)
    1/sum(w^2)}

  lower_ess <- 0.05*sum(is.finite(loss_s))

  if(get_ess(max_lambda)<lower_ess){
    max_lambda <- uniroot(function(l){get_ess(l)-lower_ess}, c(0, max_lambda))$root}

  #A sequence of learning rates
  lambda_seq <- seq(0, max_lambda, length.out=n_lambda)

  #for each bootstrap derive a lambda
  score <- matrix(NA, nrow=n_rep, ncol=n_lambda)

  for(b in seq_len(n_rep)){

    #subset
    inb_data   <- sample.int(n_samp, floor(0.8*n_samp))
    oob_data   <- setdiff(seq_len(n_samp), inb_data)

    L_inb <- apply(loss_array[inb_data,,,drop=F],3,function(x){
      x <- x[is.finite(x)]
      if(length(x)==0){NA_real_}else{mean(x)}})

    L_oob <- apply(loss_array[oob_data,,,drop=F],3,function(x){
      x <- x[is.finite(x)]
      if(length(x)==0){NA_real_}else{mean(x)}})

    ok <- (is.finite(L_inb) & is.finite(L_oob))

    if(!any(ok)){next}

    for(k in seq_along(lambda_seq)){

      w          <- .ppmn_stable_weights(L_inb[ok], lambda_seq[k])
      score[b,k] <- sum(w * L_oob[ok])}}

  avg <- colMeans(score, na.rm=T)
  if(all(!is.finite(avg))){return(0)}

  lambda_seq[which.min(avg)]}
