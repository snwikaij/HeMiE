#' Helper OOB optimizer
#'
#' @param obs_mat Observed matrix.
#' @param loss Loss matrix.
#' @param max_lambda Maximum lambda.
#' @param seed seed.
#' @param seq_len Length of lambda sequence.
#' @param n_boot number of bootstraps.
#'
#' @keywords internal
#'
.ppmn_oob_lambda_pred <- function(obs_mat, loss, max_lambda, seed, seq_len, n_boot=1000){

  set.seed(seed)

  n_samp     <- nrow(obs_mat)
  lambda_seq <- seq(0, max_lambda, length.out=seq_len)

  loss_is <- apply(loss, c(1,3), function(z){
    z <- z[is.finite(z)]
    if(length(z)==0){NA}else{mean(z)}})

  n_samp <- nrow(loss_is)

#################
#bootstrap count#
#################

  boot_count <- t(replicate(
    n_boot, tabulate(sample(seq_len(n_samp),
                            n_samp, replace=T), nbins=n_samp)))

  oob_count <- 1*(boot_count==0)

#######################
#within and oob losses#
#######################

  loss_in <- boot_count%*%loss_is
  loss_in <- loss_in/rowSums(boot_count)

  loss_oob <- oob_count%*%loss_is
  loss_oob <- loss_oob/rowSums(oob_count)

  loss_lambda <- sapply(lambda_seq, function(lambda){

    loss_b <- sapply(seq_len(n_boot), function(b){

      w <- .ppmn_stable_weights(loss_in[b,], lambda)

      ok <- is.finite(w)&is.finite(loss_oob[b,])

      w <- w[ok]
      w <- w/sum(w)

      sum(w*loss_oob[b,ok])})

    mean(loss_b, na.rm=T)})

  lambda_seq[which.min(loss_lambda)]}
