#' Helper function to detect dependent variable types
#'
#' @param y_raw Raw observations.
#'
#' @importFrom stats na.omit
#'
#' @keywords internal
#'
.ppmn_detect_type <- function(y_raw){

  #constant
  eps <- .Machine$double.eps

  #binary variables
  detect_binary <- function(x){ux <- unique(na.omit(x)); length(ux) <= 2 && all(ux %in% c(0, 1))}

  #count variables
  detect_count <- function(x){x <- na.omit(x); length(x) > 0 && all(x >= 0) && all(abs(x - round(x)) < eps)}

  #proportional data
  detect_beta <- function(x){x <- na.omit(x); if (length(x) == 0){return(F)}; in_range <- all(x >= 0 & x <= 1); many_unique <- length(unique(x)) > 2; in_range && many_unique}

  if(detect_binary(y_raw)){
    "binary"
  }else if(detect_beta(y_raw)){
    "proportional"
  }else if(detect_count(y_raw)){
    "discrete"
  }else("continuous")

}
