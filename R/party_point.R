#' Assemblage partitioning point estimation
#'
#' @param x A single numeric vector of the environmental gradient.
#' @param y A community matrix with presence/absence or abundance.
#' @param nboot Number of bootstraps (default = 1000).
#' @param level Credibility level (default = 0.9).
#' @param method The method deriving the intervals for the partitioning point.
#'
#' @keywords internal
party_point <- function(x, y, nboot = 1000, level = 0.9, method = "weighted_x"){

if(!is.numeric(x)){stop("x has to be numeric.")}
if(!any(is.matrix(y), is.data.frame(y))){stop("y has to be a matrix or data frame.")}
if(nboot<1){stop("nboot has to be >= 1.")}
if(nboot<500){warning("Make sure the number of bootstraps is high enough.")}

#order x
xord        <- order(x)
x           <- x[xord]

#extract names
taxon_names <- colnames(y)

#order y by x
y           <- as.matrix(y[xord, , drop = FALSE])

#set dimensions
n           <- nrow(y)
t           <- ncol(y)
b           <- nboot

#generate weights
w_nb        <- matrix(rgamma(nboot*n, shape = 1, rate = 1), ncol=b, nrow=n)
w_nb        <- sweep(w_nb, 2, colSums(w_nb), "/")

if(method=="weighted_x"){
#create cumulative sum for x
dx          <- c(x[1], diff(x))
dxw         <- dx*w_nb
dxw         <- sweep(dxw, 2, colSums(dxw), "/")*sum(dx)
xhat_nb     <- apply(dxw, 2, cumsum)}

#create cumulative sume for y and Y
y_ntb       <- sweep(array(y, dim = c(n, t, b)), c(1, 3), w_nb, "*")
yhat_ntb    <- apply(y_ntb, c(2, 3), function(y) cumsum(y)/sum(y))
Yhat_nb     <- apply(yhat_ntb, c(1, 3), sum)
Yhat_nb     <- sweep(Yhat_nb , 2, apply(Yhat_nb, 2, max), "/")

#select max diff
P_ntb       <- sweep(-yhat_ntb, c(1, 3), Yhat_nb, "+")
nstar_tb    <- apply(P_ntb, c(2, 3), function(p) which.max(abs(p)))

#create party point matrix
if(method=="weighted_x"){
  bidx        <- matrix(rep(seq_len(b), each = t), nrow = t)
  theta_tb    <- matrix(xhat_nb[cbind(c(nstar_tb), c(bidx))], nrow = t, ncol = b)
}else if(method=="fixed_x"){
  theta_tb <- matrix(x[nstar_tb], nrow = t, ncol = b)
}else if(method=="random_xint"){
  theta_tb <- matrix(x[nstar_tb]+runif(length(nstar_tb))*(x[pmin(nstar_tb + 1L, n)] - x[nstar_tb]), nrow = t, ncol = b)
}else{stop("Not the correct method choosen.")}

#extract mean to display as point
score_tb   <- apply(P_ntb, c(2, 3), function(p){p[which.max(abs(p))]})
cles_t      <- apply(P_ntb, 2, function(p){mean(p>0, na.rm=T)+0.5*mean(p==0, na.rm=T)})

bounds      <- c((1-level)/2, level+(1-level)/2)

summary <- data.frame(
  taxon        = taxon_names,
  cles         = cles_t,
  part_mean    = rowMeans(theta_tb),
  part_median  = apply(theta_tb, 1, median),
  part_std.err = apply(theta_tb, 1, sd),
  part_ll      = apply(theta_tb, 1, quantile, probs = bounds[1]),
  part_ul      = apply(theta_tb, 1, quantile, probs = bounds[2]),
  score_mean   = rowMeans(score_tb),
  score_median = apply(score_tb, 1, median),
  score_std.err= apply(score_tb, 1, sd),
  score_ll     = apply(score_tb, 1, quantile, probs = bounds[1]),
  score_ul     = apply(score_tb, 1, quantile, probs = bounds[2]),
  n            = as.numeric(colSums(y>0, na.rm = T)))

x_return <- if(method=="weighted_x"){xhat_nb}else if(method=="fixed_x"){
  matrix(x, nrow=n, ncol = b)}else if(method=="random_xint"){matrix(x, nrow=n, ncol = b)}

out <- list(summary=summary, theta=theta_tb,
            Yhat=Yhat_nb, yhat=yhat_ntb, x=x_return, P=P_ntb,
            taxa=taxon_names, level=level)

class(out) <- c("party_point", "list")

return(out)}
