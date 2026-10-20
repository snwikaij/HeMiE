#' Print GBMN parameter estimates
#'
#' @param x An object from gbmn_extract_param.
#' @param digits Number of decimal places.
#'
#' @keywords internal
#'
#' @export
print.gbmn.param <- function(x, digits = 3, ...) {

  print_table <- function(tab) {

    num_cols <- c("mu", "ll", "ul")
    num_cols <- intersect(num_cols, names(tab))

    tab[num_cols] <- lapply(tab[num_cols],function(z) round(as.numeric(z), digits))

    print(tab, row.names = FALSE)}

  cat("\nGBMN parameter estimates\n")
  cat("------------------------\n")

  cat("\nPrior parameters\n")
  cat("----------------\n")
  print_table(x$OldParameters)

  if(!is.null(x$UpdatedParameters)) {

    cat("\nGeneralized posterior parameters\n")
    cat("--------------------------------\n")
    print_table(x$UpdatedParameters)}

  invisible(x)}
