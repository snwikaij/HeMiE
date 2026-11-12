#' Plot the partitioning points
#'
#' @param object An party_point object
#' @param cles_cut The selection cut-off for cles.
#'
#' @importFrom cowplot plot_grid
#'
#' @keywords internal
party_plot <- function(object, cles_cut=.95){

  if(cles_cut >= 1){stop("score_cut needs to be < 1.")}
  if(cles_cut <= 0){stop("score_cut needs to be > 0.")}
  level <- object$level

  left               <- object$summary[abs(object$summary$cles-0.5)+0.5 > abs(cles_cut-0.5)+0.5,]
  if(is.null(nrow(left))){stop("Nothing selected")}
  left$sign          <- as.factor(sign(left$score_mean))

  theta_df           <- as.data.frame(object$theta)
  rownames(theta_df) <- object$taxa
  theta_df           <- theta_df[left$taxon,]

  pos_neg <- split(theta_df, left$sign)
  bounds      <- c((1-level)/2, level+(1-level)/2)

  summary_pos_neg <- lapply(pos_neg, function(x){x <- as.matrix(x);
  data.frame(
    part_mean    = mean(x),
    part_median  = median(x),
    part_std.err = sd(x),
    part_ll      = quantile(x, probs = bounds[1]),
    part_ul      = quantile(x, probs = bounds[2]))})

  summary_pos_neg      <- do.call(rbind, summary_pos_neg)
  summary_pos_neg$sign <- as.factor(rownames(summary_pos_neg))

  pl1 <-  ggplot(left, aes(x=part_mean, y=reorder(taxon, part_mean), fill=sign))+
    xlim(min(left$part_ll), max(left$part_ul))+
    geom_errorbar(xmin=left$part_ll, xmax=left$part_ul, width=0)+
    geom_point(pch=21, size=abs(left$score_mean)*6)+
    scale_fill_manual(values = c("-1"="black", "1"="white"))+
    theme_classic()+
    theme(axis.title.y = element_blank(),
          legend.position = "none",
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          plot.margin = margin(5.5, 5.5, 0, 5.5))

  pl2 <- ggplot(summary_pos_neg, aes(x=part_mean, y=sign, fill=sign))+
    xlab("Environmental gradient")+
    xlim(min(summary_pos_neg$part_ll), max(summary_pos_neg$part_ul))+
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,fill = "grey90")+
    geom_errorbar(xmin=summary_pos_neg$part_ll, xmax=summary_pos_neg$part_ul, width=0)+
    geom_point(pch=21, size=3)+
    labs(fill="Sign")+
    scale_fill_manual(values = c("-1"="black", "1"="white"),
                      labels = c("-1" = "Negative", "1" = "Positive"))+
    scale_y_discrete(labels = c("-1" = "Negative", "1" = "Positive"))+
    theme_classic()+
    theme(legend.position = "bottom",
          plot.margin = margin(0, 5.5, 5.5, 5.5),
          axis.title.y = element_blank())

  plot_obj <- cowplot::plot_grid(pl1, pl2, ncol=1, rel_heights = c(0.8,0.2),align = "v", axis = "lr")
  out      <- list(plot=plot_obj, selected=left, assemblage=summary_pos_neg)

  class(out) <- c("party_plot", "list")

  return(out)}
