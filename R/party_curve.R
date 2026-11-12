#' Plot the emperical cumulative density curves
#'
#' @param object An partitioning object created by party_point.
#' @param taxon Select a taxon name from the assemblage (default is the first).
#' @param grid The size of the grid being interpolated.
#'
#'@keywords internal
party_curve <- function(object, taxon=NULL, grid=50){

  if(is.null(taxon)){
    warning("No taxon selected first is displayed.")
    taxon        <- unique(object$taxa)[1]}
  if(grid<=10){stop("A grid with size <= 10 is not allowed.")}

  #extract data y
  Ycurves        <- object$Yhat

  #extract data y
  array_ycurves  <- object$yhat

  #extract data x
  xcurves        <- object$x

  #extract taxon names
  tax_names      <- object$taxa

  #select level
  level          <- object$level

  #xtract taxon curves
  t              <- which(tax_names %in% taxon)
  ycurves        <- array_ycurves[,t,]

  #same x grid
  xgrid <- seq(min(xcurves), max(xcurves), length.out = grid + 1)

  #interpolate curves
  interpolate <- function(z){sapply(seq_len(dim(z)[2]), function(k){
    approx(xcurves[,k], z[,k], xout=xgrid, rule=2, ties=mean)$y})}

  #same grid size
  tax_grid <- interpolate(ycurves)
  com_grid <- interpolate(Ycurves)

  #get differences
  diff_grid <- com_grid-tax_grid

  #get bootstrap quantiles functions
  summarise_curve <- function(z){
    q <- t(apply(z, 1, quantile, probs=c((1-level)/2, 0.5, (1+level)/2), na.rm=T))
    data.frame(x=xgrid, y=q[, 2], ll=q[, 1], ul=q[, 3])}

  #get all quantiles
  df_taxa       <- summarise_curve(tax_grid)
  df_assemblage <- summarise_curve(com_grid)
  df_disc       <- summarise_curve(diff_grid)

  #plot the curves for Yhat and yhat
  pl1 <- ggplot(df_taxa, aes(x, y))+
    geom_line(col="dodgerblue4")+ylab("F(x)")+
    geom_ribbon(aes(ymin = ll, ymax = ul), fill = "dodgerblue", alpha = 0.2)+
    geom_line(data = df_assemblage, aes(x=x, y=y), inherit.aes = F, col="tomato3") +
    geom_ribbon(data=df_assemblage, aes(x=x, ymin=ll, ymax=ul), inherit.aes=F, fill="tomato", alpha=0.2)+
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())

  #plot the differences between both curves
  pl2 <- ggplot(df_disc, aes(x, y))+
    geom_line(col="dodgerblue4")+ylab(expression(P == F[T](x) - F[t](x)))+xlab("Environmental gradient")+
    geom_ribbon(aes(ymin = ll, ymax = ul), fill = "dodgerblue", alpha = 0.2)+
    geom_hline(yintercept = 0, lty=2, col = "tomato3")+
    theme_classic()

  combined   <- cowplot::plot_grid(pl1, pl2, ncol=1, align = "hv")
  out        <- list(plot=combined, plot1=pl1, plot=pl2)
  class(out) <- c("party_curve", "list")
  combined}
