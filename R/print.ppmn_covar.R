#' Print PPMN correlation matrix summary
#'
#' @param x An object returned by extract_covar_ppmn.
#'
#' @keywords internal
#'
#' @export
print.ppmn.covar <- function(x, ...) {

  cat(sprintf("Correlation matrix with %d x %d parameters.\n",
    nrow(x$cor_mat),
    ncol(x$cor_mat)))

  invisible(x)}
