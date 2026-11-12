#' Print party_point results
#'
#' @param x A party_point object.
#' @param ... Additional arguments.
#'
#' @keywords internal
#'
#' @export
print.party_point <- function(x, ...) {

  cat("\nPartitioning summary\n")
  cat("--------------------\n")

  cat("\nTotal taxa:", nrow(x$summary), "\n")
  cat("Bootstrap replicates:", ncol(x$theta), "\n")
  cat("Credibility level:", x$level, "\n")

  cat("\nTaxa summary\n")
  cat("------------\n")

  summary <- x$summary

  if (nrow(summary) <= 10) {
    print(summary, row.names = FALSE)
  } else {
    print(head(summary, 5), row.names = FALSE)
    cat("...", nrow(summary) - 5, "additional rows\n")
  }

  invisible(x)
}
