#' Helper function for Bayesian bootstrap
#'
#' @param resid1 Residuals for particle with best fit.
#' @param resid0 Residuals for the prior.
#' @param n_grid Size of the grid.
#' @param n_boot Number of bootstraps.
#'
#' @keywords internal
#'
.ppmn_resid_bboot <- function(resid1, resid0, n_grid=3000, n_boot=1000){

  #constant
  eps <- .Machine$double.eps

  #set resids as numeric
  resid1 <- as.numeric(resid1)
  resid0 <- as.numeric(resid0)

  #select finite
  resid1 <- resid1[is.finite(resid1)]
  resid0 <- resid0[is.finite(resid0)]

  if(length(resid1)<2||length(resid0)<2){return(0)}

  #nsamp
  n_samp <- length(resid1)

  #############
  #Bayes bootp#
  #############

  #weights
  weights <- rgamma(n_samp*n_boot, shape=1, rate=1)
  w_boot  <- matrix(weights, nrow=n_samp, ncol=n_boot)
  w_boot  <- sweep(w_boot, 2,colSums(w_boot), "/")

  #resid boots
  r_boot <- matrix(resid1, nrow=n_samp, ncol=n_boot)

  #mean resids
  mean_resid_boot <- colSums(r_boot*w_boot)

  #trim edges
  q_prior      <- quantile(resid0, c(.005,.995), na.rm=T)
  prior_res    <- resid0[resid0>=q_prior[1] & resid0<=q_prior[2]]

  bw_data_mean <- density(mean_resid_boot, kernel="gaussian", bw="nrd0")$bw
  bw_prior     <- density(prior_res, kernel="gaussian", bw="nrd0")$bw

  #create grid
  q_grid <- quantile(c(mean_resid_boot, prior_res), c(.01,.99),na.rm=T)

  kd_prior <- density(prior_res, kernel="gaussian", bw=bw_prior, n=n_grid, from=q_grid[1], to=q_grid[2])
  kd_data  <- density(mean_resid_boot, kernel="gaussian", bw=bw_data_mean, n=n_grid, from=q_grid[1], to=q_grid[2])

  p_prior <- kd_prior$y
  p_data  <- kd_data$y
  grid    <- kd_prior$x
  dx      <- grid[2]-grid[1]

  p_prior <- p_prior/sum(p_prior*dx)
  p_data  <- p_data/sum(p_data*dx)

  p_post  <- pmax(p_prior,eps)*pmax(p_data,eps)
  p_post  <- p_post/sum(p_post*dx)

  KL_target <- sum(p_post*(log(pmax(p_post,eps))-log(pmax(p_prior,eps)))*dx)
  if(!is.finite(KL_target)||KL_target<0){KL_target <- 0}

  KL_target}
