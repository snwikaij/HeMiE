#' Helper function to residual diagnostics
#'
#' @param object An updated PPMN.
#' @param data The data used to update the PPMN.
#' @param mad_constant Constant.
#'
#' @keywords internal
#'
#' @importFrom stats mad median na.omit quantile aggregate
#'
.ppmn_resid_diagnostic <- function(object, data, mad_constant=1.4826){

  eps <- .Machine$double.eps

  preds_resid <- predict_ppmn(object,data,nsim=1)

  obs_nodes <- intersect(
    setdiff(names(preds_resid$Expected),object$Roots),
    colnames(data))

  if(length(obs_nodes)==0){
    stop("No observed child vertices in data match the model.")}

  obs_mat <- as.matrix(data)

  x_hat           <- do.call(cbind,preds_resid$Expected[obs_nodes])
  x_hat           <- as.matrix(x_hat)
  colnames(x_hat) <- obs_nodes

  mat <- matrix(NA,nrow=nrow(obs_mat),ncol=length(obs_nodes))
  colnames(mat) <- obs_nodes

  for(j in seq_along(obs_nodes)){

    x_raw <- obs_mat[,obs_nodes[j]]
    y_hat <- x_hat[,j]
    type  <- .ppmn_detect_type(x_raw)

    if(type %in% c("binary","proportional")){

      p_hat    <- pmin(pmax(y_hat,eps),1-eps)
      mat[,j] <- (x_raw-p_hat)/sqrt(p_hat*(1-p_hat)+eps)

    }else if(type=="discrete"){

      mat[,j] <- 2*(sqrt(x_raw+3/8)-sqrt(pmax(y_hat,eps)+3/8))

    }else{

      mat[,j] <- x_raw-y_hat

      MAD <- mad(
        mat[,j],
        center=median(mat[,j],na.rm=T),
        constant=mad_constant,
        na.rm=T)

      if(!is.finite(MAD) || MAD==0){MAD <- eps}

      mat[,j] <- mat[,j]/MAD}
  }

  res_medsd <- apply(mat,2,function(x){
    c("median"=median(x,na.rm=T),
      "mu"=mean(x,na.rm=T))})

  forlabel <- round(res_medsd,2)
  labels   <- apply(forlabel,2,function(x){
    paste0("\nmedian=",x[1],"\nmu=",x[2])})

  labels <- paste0(obs_nodes,labels)

  list_long <- lapply(seq_len(ncol(mat)),function(x){
    data.frame(labels[x],na.omit(mat[,x]))})

  df_long <- do.call(rbind,list_long)
  colnames(df_long) <- c("Variable","Value")

  qy <- quantile(df_long$Value,c(.025,.975),na.rm=T)

  resid_boxplot <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x=Variable,y=Value))+
    ggplot2::ylim(qy)+
    ggplot2::geom_hline(yintercept=0, col="tomato3",lty=2)+
    ggplot2::ylab("Standardized residuals")+
    ggplot2::geom_boxplot(outlier.shape=NA)+
    ggplot2::geom_jitter(alpha=.2,pch=19,width=.2)+
    ggplot2::theme_classic()

  x_hat_df   <- as.data.frame(x_hat)
  x_hat_df$i <- seq_len(nrow(x_hat_df))

  x_hat_long <- tidyr::pivot_longer(
    x_hat_df,
    cols=-i,
    names_to="j",
    values_to="x_hat")

  obs_df   <- as.data.frame(obs_mat[,obs_nodes,drop=F])
  obs_df$i <- seq_len(nrow(obs_df))

  obs_rank <- obs_df
  obs_rank[,obs_nodes] <- lapply(obs_rank[,obs_nodes,drop=F],rank)

  x_hat_rank <- x_hat_df
  x_hat_rank[,obs_nodes] <- lapply(x_hat_rank[,obs_nodes,drop=F],rank)

  obs_long <- tidyr::pivot_longer(
    obs_df,
    cols=-i,
    names_to="j",
    values_to="x")

  obs_rank_long <- tidyr::pivot_longer(
    obs_rank,
    cols=-i,
    names_to="j",
    values_to="obs_rank")

  x_hat_rank_long <- tidyr::pivot_longer(
    x_hat_rank,
    cols=-i,
    names_to="j",
    values_to="x_hat_rank")

  x_hat_long$obs        <- obs_long$x
  x_hat_long$rank_obs   <- obs_rank_long$obs_rank
  x_hat_long$rank_x_hat <- x_hat_rank_long$x_hat_rank

  resid_pred_vs_fit <- ggplot2::ggplot(
    x_hat_long,
    ggplot2::aes(x=x_hat,y=obs))+
    ggplot2::geom_point()+
    ggplot2::xlab("Predicted")+
    ggplot2::ylab("Observed")+
    ggplot2::theme_classic()+
    ggplot2::geom_abline(intercept=0, slope=1, colour="tomato3", lwd=.6)+
    ggplot2::facet_wrap(.~j,scales="free")

  object$Residuals$list                     <- df_long
  object$Residuals$resid_boxplot            <- resid_boxplot
  object$Residuals$pred_observed            <- resid_pred_vs_fit

  return(object)}
