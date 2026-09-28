#' Helper function for root estimation
#'
#' @param x An object to summarize.
#'
#' @keywords internal
#'
#' @export
print.ppmn_root_estimation <- function(x, ...){

  cat("\nPPMN root estimation\n")
  cat("--------------------\n\n")

  tab <- x$summary

  #move row names into an explicit column if present
  if(!is.null(rownames(tab)) &&
     !all(rownames(tab) == seq_len(nrow(tab)))){

    tab <- data.frame(Root = rownames(tab), tab, row.names = NULL, check.names = F)}

  print(tab, row.names = F)

  invisible(x)
}
