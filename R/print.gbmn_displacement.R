#' Helper function to print GBMN displacement
#'
#' @param x An object returned by disp_gbmn.
#'
#' @keywords internal
#'
#' @export
print.gbmn_disp <- function(x, ...){

  cat("\nDisplacement\n")
  cat("-----------------\n\n")

  tab <- x$results

  #round numeric values
  num_cols <- c("median","map","mu","ll","ul")
  num_cols <- intersect(num_cols, names(tab))

  tab[num_cols] <- lapply(tab[num_cols], function(z){round(as.numeric(z),3)})

  print(tab, row.names=F)

invisible(x)}
