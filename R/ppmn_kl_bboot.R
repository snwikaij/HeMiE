#' Helper function for multivariate Bayesian bootstrap KL calibration
#'
#' @param resid1 residual matrix for best particle.
#' @param resid0 residual array from prior.
#' @param n_grid size of the grid in grid approximation.
#' @param n_boot number of Bayesian bootstrap replicates.
#' @param shrink covariance shrinkage toward a diagonal matrix.
#'
#' @keywords internal
#'
.ppmn_kl_bboot <- function(resid1, resid0, method="multikernel", n_grid=5000, n_boot=1000, shrink=0.1){

  #constant
  eps <- .Machine$double.eps

  #dimensions
  d <- dim(resid0)
  if(length(d)!=3){stop("resid0 must be an [i,j,s] array.")}

  #select each dim
  n_samp <- d[1]
  n_vert <- d[2]
  n_sim  <- d[3]

###################################################
#prior and observed distribution of residual means#
###################################################

  #use medians to increase outlier robustness
  resid_0 <- t(apply(resid0, c(2, 3), function(x) median(x, na.rm = T)))

  resid_1 <- apply(resid1, c(1, 2), function(x) median(x, na.rm = T))

########################
#multivariate Bayes BB#
########################

  #derive weight per ith observation
  w_boot <- matrix(rgamma(n_samp*n_boot,shape=1,rate=1),
                   nrow=n_samp,ncol=n_boot)
  w_boot <- sweep(w_boot,2,colSums(w_boot),"/")

  #handle missing residuals
  observed           <- is.finite(resid_1)
  resid_1[!observed] <- 0

  #multiply w and resied
  num <- crossprod(w_boot,resid_1)

  #if not observed 0
  den <- crossprod(w_boot, observed*1)

  #results of bb
  boot <- num/den
  boot[den==0] <- NA

################
#create kernels#
################

  if(method == "multikernel"){
  #select bandwidths
  bw_pr_kernel <- apply(resid_0, 2, function(r) density(r, kernel="gaussian", bw="nrd0")$bw)
  bw_bb_kernel <- apply(boot, 2, function(b) density(b, kernel="gaussian", bw="nrd0")$bw)

  #create grid
  q_grid <- quantile(c(resid_0, boot), c(.01,.99), na.rm=T)

  #create kernels
  kd_prior <- lapply(seq_len(n_vert), function(r){
    density(resid_0[,r], kernel="gaussian", bw=bw_pr_kernel[r],
            n=n_grid, from=q_grid[1], to=q_grid[2])})
  kd_psuedolike <- lapply(seq_len(n_vert), function(r){
    density(boot[,r], kernel="gaussian", bw=bw_bb_kernel[r],
            n=n_grid, from=q_grid[1], to=q_grid[2])})

  #derive kl via grid
  KLt <- sapply(seq_len(n_vert), function(p){

    p_prior <- pmax(kd_prior[[p]]$y,eps)
    p_data  <- pmax(kd_psuedolike[[p]]$y,eps)

    grid    <- kd_prior[[p]]$x
    dx      <- grid[2]-grid[1]

    #normalize densities
    p_prior <- p_prior/sum(p_prior*dx)
    p_data  <- p_data/sum(p_data*dx)

    #calculate posterior
    p_post  <- p_prior*p_data
    p_post  <- p_post/max(p_post)
    p_post  <- p_post/sum(p_post*dx)

    #KL div
    KL <- sum(p_post*(log(p_post)-log(p_prior))*dx)

    KL})

  KLt <- sum(KLt)/n_vert

  }else if(method == "multivar"){

    #prior mean and covar
    mu_prior   <- colMeans(resid_0)
    Sig_prior  <- cov(resid_0)

    #bb mean and covar
    mu_bb      <- colMeans(boot)
    Sig_bb     <- cov(boot)

    #########################
    #analytical MVN posterior
    #########################

    #inverse matrix to get precision
    prec_mat_pr <- solve(Sig_prior)
    prec_mat_bb <- solve(Sig_bb)

    #posterior covariance
    post_covar  <- solve(prec_mat_pr+prec_mat_bb)

    #posterior mean
    mu_post <- drop(post_covar %*% (prec_mat_pr %*% mu_prior + prec_mat_bb%*%mu_bb))

    #########################
    #analytical KL divergence
    #########################

    dif <- mu_post-mu_prior

    KL <- 0.5*(
      sum(diag(prec_mat_pr%*%post_covar))+
        drop(crossprod(dif,prec_mat_pr%*%dif))-
        n_vert+
        as.numeric(determinant(Sig_prior,logarithm=T)$modulus)-
        as.numeric(determinant(post_covar,logarithm=T)$modulus))

    KLt <- KL/n_vert}else{stop("Either multikernel or multivar.")}

    KLt}
