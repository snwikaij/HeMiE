#' GBMN correlation matrix extraction
#'
#' @param object An updated or foundational GBMN.
#'
#' @description
#' This function extract the covariance matrix (Sigma) from an GBMN object. This is to make your life easier
#' It also creates an correlation matrix and biplot of a principal component analysis on this matrix.
#'
#' @export
extract_covar_gbmn <- function(object){

  cor_mat <- cov2cor(object$Sigma)
  pca_display <- eigen(cor_mat)

  pca_df <- sweep(
    pca_display$vectors[, 1:2, drop = FALSE],
    2,
    sqrt(pmax(pca_display$values[1:2], 0)),
    "*")

  rownames(pca_df) <- rownames(cor_mat)
  pca_df           <- data.frame(pc1=pca_df[,1], pc2=pca_df[,2], name=sub("^[^e]*e", "", rownames(pca_df)))

  pca_biplot <- ggplot(pca_df, aes(pc1, pc2)) +
    geom_segment(
      aes(x = 0, y = 0, xend = pc1, yend = pc2),
      arrow = arrow(length = grid::unit(0.2, "cm")),
      alpha=0.25) +
    geom_vline(xintercept = 0, lty=2, alpha=0.2)+
    geom_hline(yintercept = 0, lty=2, alpha=0.2)+
    geom_text(aes(label = name), hjust = -0.1, vjust = -0.3) +
    theme_classic()

  df           <- as.data.frame(as.table(cor_mat))
  colnames(df) <- c("Parameter1", "Parameter2", "Correlation")
  df           <- df[upper.tri(cor_mat), ]
  minmax       <- range(cor_mat[upper.tri(cor_mat)])

  up_diagcor <- ggplot(df, aes(Parameter1, Parameter2, fill = Correlation)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "magenta", mid = "yellow", high = "cyan",
      midpoint = 0, limits = c(minmax[1], minmax[2])) +
    coord_equal() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    labs(x = NULL, y = NULL, fill = "Correlation")+
    geom_text(aes(label = round(Correlation, 2)), size = 3)

  out <- list(cor_mat=cor_mat, corplot=up_diagcor, biplot=pca_biplot)

  class(out) <- c("gbmn.covar", "list")

  return(out)}
