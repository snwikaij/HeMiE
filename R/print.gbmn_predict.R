#' Helper function for predictions
#'
#' @param x An object to summarize.
#'
#' @importFrom utils head
#'
#' @keywords internal
#'
#' @export
print.gbmn_prediction <- function(x, ...){

  cat("\nGBMN prediction\n")
  cat("---------------\n")

  pred_nodes <- setdiff(names(x$Expected))

  cat("Predicted nodes :", paste(pred_nodes, collapse = ", "), "\n")

######################
#expected predictions#
######################

  expected <- as.data.frame(x$Expected, check.names = FALSE)

  cat("\nExpected predictions\n")
  cat("--------------------\n")

  if(nrow(expected) <= 10){
    print(expected, row.names = FALSE)
  }else{
    print(head(expected, 5), row.names = FALSE)
    cat("...", nrow(expected)-5, "additional rows\n")}

  invisible(x)}
