#' Extract parameters from object
#'
#' @param object An updated or foundational PPMN.
#'
#' @description
#' This is a simple function that creates a data frame with the point, lower and upper interval. It also produce two figures
#' of this data frame so that the change in the estimated parameters can be checked. This function is build to make your life easier.
#' You do not have to use the figures, but it makes extraction a bit easier.
#'
#' @export
extract_param_ppmn <- function(object){

   if(object$Update>=1){

    OldParameters     <- do.call(rbind, lapply(names(object$OldParameters), function(ef){

    efun <- object$OldParameters[[ef]]

    do.call(rbind, lapply(names(efun), function(p){
      data.frame(ef=ef, parameter=p,
        mu=efun[[p]][["mu"]], ll=efun[[p]][["ll"]], ul=efun[[p]][["ul"]])}))}))

    UpdatedParameters <- do.call(rbind, lapply(names(object$Parameters), function(ef){

    efun <- object$Parameters[[ef]]

    do.call(rbind, lapply(names(efun), function(p){
      data.frame(ef=ef, parameter=p,
                 mu=efun[[p]][["mu"]], ll=efun[[p]][["ll"]], ul=efun[[p]][["ul"]])}))}))

    Old_params      <- OldParameters
    New_params      <- UpdatedParameters
    Old_params$type <- "Prior"
    New_params$type <- "Generalized posterior"

    df_params       <- rbind(Old_params, New_params)
    df_params$type  <- factor(df_params$type, levels=c("Prior", "Generalized posterior"))
    df_params$ef    <- factor(df_params$ef, levels=unique(df_params$ef))
  }else{

    OldParameters   <- do.call(rbind, lapply(names(object$Parameters), function(ef){

      efun <- object$Parameters[[ef]]

      do.call(rbind, lapply(names(efun), function(p){
        data.frame(ef=ef, parameter=p,
                   mu=efun[[p]][["mu"]], ll=efun[[p]][["ll"]], ul=efun[[p]][["ul"]])}))}))
    df_params       <- OldParameters
    df_params$type  <- "Prior"

    df_params$type  <- factor(df_params$type, levels=c("Prior"))
    df_params$ef    <- factor(df_params$ef, levels=unique(df_params$ef))
  }

simple_plot_fun <- function(df_params){

dodge <- position_dodge(width = 0.5)

pl1 <- ggplot(df_params, aes(x=parameter, y=mu, col=type, group=type))+
  geom_point(position=dodge)+ylab("Estimate")+
  scale_color_manual(breaks = c("Prior", "Generalized posterior"), values=c("orange2", "dodgerblue3"))+
  geom_errorbar(aes(x=parameter, ymin=ll, ymax=ul),
                width = 0, position=dodge)+
  facet_wrap(.~ef+parameter, scale="free", ncol=3)+
  theme_classic()+
  theme(strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1)),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "bottom")

pl2 <- ggplot(df_params, aes(x=parameter, y=mu, col=type, group=type))+
  geom_point(position=dodge)+ylab("Estimate")+
  scale_color_manual(breaks = c("Prior", "Generalized posterior"), values=c("orange2", "dodgerblue3"))+
  geom_errorbar(aes(x=parameter, ymin=ll, ymax=ul),
                width = 0, position=dodge)+
  facet_wrap(.~ef, scale="free", ncol=3)+
  theme_classic()+
  theme(strip.text = element_text(size = 6, margin = margin(1, 1, 1, 1)),
        strip.background = element_rect(linewidth = 0.2),
        legend.position = "bottom")

list(pl1=pl1, pl2=pl2)}

res <- simple_plot_fun(df_params)

out <- list(plot=list(plot1 = res$pl1, plot2 = res$pl2),
  OldParameters=OldParameters,
  UpdatedParameters=if(object$Update >= 1){
    UpdatedParameters}else{
    NULL})

class(out) <- c("ppmn.param", "list")

return(out)}
