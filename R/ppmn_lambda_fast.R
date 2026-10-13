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
                              n_boot=20, n_lambda=100, ess_frac=0.1){

  #loss per per [i,s]
  loss_is <- apply(loss_array, c(1,3), function(x){
    x <- x[is.finite(x)]; if(length(x)==0){NA}else{mean(x)}})

  #total loss per [s]
  loss_s <- apply(loss_array, 3, function(x){
    x <- x[is.finite(x)]
    if(length(x)==0){NA_real_}else{mean(x)}})

  #minimum effective particle count
  ess_min <- ceiling(ess_frac*length(is.finite(loss_s)))

  #get the ESS to prevent collapse
  get_ess <- function(l){
    w <- .ppmn_stable_weights(loss_s, l)
    1/sum(w^2)}

  if(get_ess(max_lambda)<ess_min){
    max_lambda <- uniroot(function(l){get_ess(l)-ess_min}, c(0, max_lambda))$root}

  #A sequence of learning rates
  lambda_seq <- seq(0, max_lambda, length.out=n_lambda)

  #for each bootstrap derive a lambda
  score <- matrix(NA, nrow=n_boot, ncol=n_lambda)

  for(b in seq_len(n_boot)){

    #subset
    inb_data   <- sample.int(nrow(loss_is), floor(0.8*nrow(loss_is)))
    oob_data   <- setdiff(seq_len(nrow(loss_is)), inb_data)

    L_inb <- colMeans(loss_is[inb_data,,drop=F], na.rm=T)
    L_oob <- colMeans(loss_is[oob_data,,drop=F], na.rm=T)

    ok <- (is.finite(L_inb) & is.finite(L_oob))

    if(!any(ok)){next}

    for(k in seq_along(lambda_seq)){

      w          <- .ppmn_stable_weights(L_inb[ok], lambda_seq[k])
      score[b,k] <- sum(w * L_oob[ok])}}

  avg <- colMeans(score, na.rm=T)

  if(all(!is.finite(avg))){return(0)}

  lambda_seq[which.min(avg)]}
