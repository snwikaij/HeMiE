#' Helper function for predictions
#'
#' @param x An object to summarize.
#'
#' @importFrom utils head
#'
#' @export
print.ppmn_prediction <- function(x, ...){

  cat("\nPPMN prediction\n")
  cat("---------------\n")

  cat("Simulations     :", x$nsim, "\n")
  cat("Roots           :", paste(x$Roots, collapse = ", "), "\n")

  pred_nodes <- setdiff(
    names(x$Expected),
    x$Roots
  )

  cat(
    "Predicted nodes :",
    paste(pred_nodes, collapse = ", "),
    "\n"
  )

  #####################
  #expected predictions
  #####################

  expected <- as.data.frame(
    x$Expected,
    check.names = FALSE
  )

  cat("\nExpected predictions\n")
  cat("--------------------\n")

  if(nrow(expected) <= 10){

    print(
      expected,
      row.names = FALSE
    )

  }else{

    print(
      head(expected, 6),
      row.names = FALSE
    )

    cat(
      "...",
      nrow(expected) - 6,
      "additional rows\n"
    )
  }

  invisible(x)
}
