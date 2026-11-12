#' Print party_plot results
#'
#' @param x A party_plot object.
#' @param ... Additional arguments.
#'
#' @keywords internal
#'
#' @export
print.party_plot <- function(x, ...) {

  cat("\nPartitioning summary\n")
  cat("--------------------\n")

  cat("\nPositive and Negative\n")
  print(x$assemblage, row.names = FALSE)

  cat("\nSelected taxa\n")
  cat("--------------------\n")

  selected <- x$selected

  if (nrow(selected) <= 10) {
    print(selected, row.names = FALSE)
  } else {
    print(head(selected, 5), row.names = FALSE)
    cat("...", nrow(selected) - 5, "additional rows\n")
  }

  print(x$plot)

  invisible(x)}
