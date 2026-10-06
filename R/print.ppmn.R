#' Helper function to print summary table
#'
#' @param x An object to summarize.
#'
#' @keywords internal
#'
#' @export
print.ppmn <- function(x, ...){

  #number of vertices
  n_vertices <- igraph::vcount(
    x$Structure$visual_dag$graph)

  #number of unique synthesized parameters
  n_parameters <- sum(
    vapply(
      x$Parameters,
      length,
      integer(1)))

  #total number of parameter estimates synthesized
  n_estimates <- sum(
    unlist(
      lapply(
        x$Parameters,
        function(edge){
          vapply(
            edge,
            function(par) par["n"],
            numeric(1)
          )
        }
      )
    ),
    na.rm = TRUE
  )

  cat("\nPosterior Predictive Meta-Analytic Network\n")
  cat("-------------------------------------------\n\n")

  cat("Vertices               :", n_vertices, "\n")
  cat("Synthesized parameters :", n_parameters, "\n")
  cat("Parameter estimates    :", n_estimates, "\n")
  cat("Updates                :", x$Update, "\n")


  #################################
  # update information if present #
  #################################

  if(
    !is.null(x$Info) &&
    !is.null(x$Info$lambda) &&
    !is.null(x$Info$ESS)
  ){

    cat("\nGeneralized Bayesian update\n")
    cat("---------------------------\n")
    cat("Method   :", x$Info$method, "\n")

    update_table <- data.frame(
      Lambda = round(
        as.numeric(x$Info$lambda),
        3),
        ESS = round(
        as.numeric(x$Info$ESS),
        1
        )
      )

      print(
        update_table,
        row.names = FALSE
      )
  }


  ###################
  # parameter table #
  ###################

  cat("\n")

  par_table <- x$ParameterTable

  num_cols <- c("mu","se","ll","ul","Q","I2","a","b")
  num_cols <- intersect(num_cols, names(par_table))

  par_table[num_cols] <- lapply(
    par_table[num_cols],
    function(z){round(as.numeric(z),3)})

  print(
    par_table,
    row.names=FALSE
  )

  invisible(x)
}
