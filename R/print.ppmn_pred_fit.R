#' Helper function for predictive fit
#'
#' @param x An object to summarize.
#'
#' @export
print.ppmn_pred_fit <- function(x, ...){

  print_tab <- function(tab){

    tab <- data.frame(
      Variable = rownames(tab),
      tab,
      row.names = NULL,
      check.names = FALSE
    )

    print(
      tab,
      row.names = FALSE
    )
  }

  cat("\nPPMN predictive fit\n")
  cat("-------------------\n")

  cat("\nBalance\n")
  cat("-------\n")
  print_tab(x$Balance)

  cat("\nCorrelation\n")
  cat("-----------\n")
  print_tab(x$Correlation)

  cat("\nStandardized residual error\n")
  cat("---------------------------\n")
  print_tab(x$Stdz)

  invisible(x)
}
