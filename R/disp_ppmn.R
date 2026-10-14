#' Displacement test for PPMN
#'
#' @param object A foundational or updated PPMN.
#' @param new_data New dataset with all root nodes and least one child node present.
#' @param nsim Number of simulations used to construct the intervals.
#' @param level Credibility level (default level = 0.9).
#' @param seed Seed value 123.
#'
#' @description
#' The displacement (D) quantifies the extent to which the parameters of current PPMN must change in order
#' to maximally align with the new data. The learning rate  is set sufficiently high so the centrality of the generalized
#' posterior is shifted fully dominated by the data. In comparison, the model predictions prior to the update 'behave' as not having
#' learned from the new data. The displacement measures how much information the data provides to shift the expected value relative
#' to the prediction and spread of the current PPMN. A value close to 0 indicates the data provides relatively limited information
#' that would result in the displacement of the expected value.
#'
#' @export
disp_ppmn <- function(object, new_data, nsim=3000, level=0.9, seed=123){

  pred1    <- predict_ppmn(object, new_data, nsim = nsim)
  ppmn_100 <- suppressWarnings(update_ppmn(object, up_character = 0, new_data=new_data, method="manual", nsim=nsim))
  pred_100 <- predict_ppmn(ppmn_100, new_data, nsim = nsim)

  #set seed
  set.seed(seed)

  #select predicted nodes
  pred_nodes_all <- setdiff(names(pred1$Variance), pred1$Roots)

  #select observed nodes/variables
  obs_nodes <- intersect(colnames(new_data), pred_nodes_all)

  #stack all prediction for standard mod
  pred1_stack <- do.call(rbind,
                         lapply(seq_len(nsim), function(i) {
                           do.call(cbind, lapply(pred1$Variance[obs_nodes], `[[`, i))}))

  #stack all predictions for updated post
  pred100_stack <- do.call(rbind,
                           lapply(seq_len(nsim), function(i) {
                             do.call(cbind, lapply(pred_100$Variance[obs_nodes], `[[`, i))}))

  #substract difference
  diffpred <- pred1_stack-pred100_stack
  diffpred <- na.omit(diffpred)

  #standardize per standard mod
  sdpred1  <- apply(pred1_stack, 2, function(x) sd(x, na.rm = T))

  #divide diff so by sd standard model
  stdzdiff <- sweep(diffpred, 2, sdpred1, "/")

  #map
  maxpost  <- function(x){d <- density(x); d$x[which.max(d$y)]}

  #density curves
  forestplot <- do.call(rbind, lapply(seq_len(ncol(stdzdiff)), function(i){
    d <- stdzdiff[,i]
    q <- quantile(d, c((1-level)/2, 1-(1-level)/2))
    data.frame(Variable=obs_nodes[i], median=median(d), map=maxpost(d), mu=mean(d), ll=q[1], ul=q[2])}))

  #rownames
  rownames(forestplot) <- NULL

  #create final results
  plot <- ggplot(forestplot, aes(Variable, median)) +
    geom_hline(yintercept = 0, linetype = 2,
               linewidth = .6, colour = "tomato3") +
    geom_point()+
    geom_errorbar(aes(ymin = ll, ymax = ul),
                  width = 0, linewidth = 0.5)+
    labs(x = NULL, y = "Displacement") +
    theme_classic()

  out <- list(plot=plot, results=forestplot)

  class(out) <- c("ppmn_disp","list")

  invisible(out)}
