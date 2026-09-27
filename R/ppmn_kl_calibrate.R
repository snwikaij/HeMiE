#' Grid approximation for Kullenback-Leibler divergence
#'
#' @param r_data Residuals of the data.
#' @param r_prior Residuals of the prior.
#' @param n_boot Number of bootstraps.
#' @param n_grid Size of the grid.
#'
#' @importFrom stats density rnorm quantile
#'
.ppmn_kl_calibrate <- function(r_data, r_prior, n_boot=2000, n_grid){

  eps <- .Machine$double.eps

  r_data  <- as.matrix(r_data)
  r_prior <- as.matrix(r_prior)

  keep_i  <- apply(r_data,1,function(x) any(is.finite(x)))
  r_data  <- r_data[keep_i,,drop=F]
  r_prior <- r_prior[keep_i,,drop=F]

  n  <- nrow(r_data)
  rd <- as.numeric(r_data)
  rp <- as.numeric(r_prior)
  rd <- rd[is.finite(rd)]
  rp <- rp[is.finite(rp)]

  if(length(rd)<2||length(rp)<2){return(0)}

  bw_data <- density(rd, kernel="gaussian", bw="nrd0")$bw

  r_boot <- sample(rd, size=n*n_boot, replace=T)
  r_boot <- r_boot+rnorm(n*n_boot, mean=0, sd=bw_data)
  r_boot <- matrix(r_boot, nrow=n, ncol=n_boot)

  mean_resid_boot <- colMeans(r_boot, na.rm=T)

  bw_data_mean <- density(mean_resid_boot, kernel="gaussian", bw="nrd0")$bw
  bw_prior     <- density(rp, kernel="gaussian", bw="nrd0")$bw

  q_grid <- quantile(c(mean_resid_boot, rp), c(.01,.99), na.rm=T)

  kd_prior <- density(rp, kernel="gaussian", bw=bw_prior,n=n_grid, from=q_grid[1], to=q_grid[2])
  kd_data  <- density(mean_resid_boot, kernel="gaussian", bw=bw_data_mean, n=n_grid, from=q_grid[1],to=q_grid[2])

  p_prior <- kd_prior$y
  p_data  <- kd_data$y
  grid    <- kd_prior$x
  dx      <- grid[2]-grid[1]

  p_prior <- p_prior/sum(p_prior*dx)
  p_data  <- p_data/sum(p_data*dx)

  p_post <- pmax(p_prior,eps)*pmax(p_data,eps)
  p_post <- p_post/max(p_post)
  p_post <- p_post/sum(p_post*dx)

  KL_target <- sum(p_post*(log(p_post)-log(p_prior))*dx)
  if(!is.finite(KL_target) || KL_target<0){KL_target <- 0}

  KL_target}
