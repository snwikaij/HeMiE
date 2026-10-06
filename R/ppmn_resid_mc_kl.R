#' Helper MC approx. for mv KL-divergence
#'
#' @param resid1 Residuals.
#' @param resid0 Model prior residuals
#' @param loss Loss.
#' @param n_boot Number of bayesian bootstraps.
#' @param n_mc Number of mc.
#'
#' @keywords internal
#'
.ppmn_resid_mc_kl <- function(resid1, resid0, loss, n_boot=1000, n_mc=5000){

  #best fit particle
  best     <- which.min(loss)
  best_res <- resid1[,,best]

  #select j,s resids j=vertices, s=3000=nsim
  prior_res <- apply(resid0, c(2,3),
    function(z){
      z <- z[is.finite(z)]
      if(length(z)==0){NA}else{mean(z)}})

  #particles x nodes
  prior_res <- t(prior_res)

  #bayes boot
  i_samp <- nrow(best_res)
  j_vert <- ncol(best_res)

  #weights for bboot
  weights <- rgamma(i_samp*n_boot, shape=1, rate=1)
  w_boot  <- matrix(weights, nrow=i_samp, ncol=n_boot)
  w_boot  <- sweep(w_boot, 2, colSums(w_boot), "/")

  #one vec per bootstrap
  data_res <- t(sapply(seq_len(n_boot), function(b){colSums(best_res*w_boot[,b])}))

  #set bw
  bw_prior <- apply(prior_res, 2, bw.nrd0)
  bw_data  <- apply(data_res, 2, bw.nrd0)

  #boot the resids
  id <- sample(seq_len(nrow(prior_res)), n_mc, replace=T)
  mc_res <- prior_res[id,,drop=F]

  log_fdata <- sapply(seq_len(n_mc), function(m){

      z <- sweep(data_res, 2, mc_res[m,], "-")
      z <- sweep(z, 2, bw_data, "/")

      a  <- -0.5*rowSums(z^2)-sum(log(bw_data))-j_vert*log(2*pi)/2
      mx <- max(a)
      mx+log(mean(exp(a-mx)))})

  mx    <- max(log_fdata)
  log_Z <- mx+log(mean(exp(log_fdata-mx)))

  a <- log_fdata-max(log_fdata)
  w <- exp(a)
  w <- w/sum(w)

  KL_target <- sum(w*(log_fdata-log_Z))

  #no uppy if kl is inf
  if(!is.finite(KL_target)||KL_target<0){
  KL_target <- 0}

  KL_target}
