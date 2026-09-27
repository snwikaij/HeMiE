#' Helper function to print summary table
#'
#' @param x An object to summarize.
#'
#' @export
print.ppmn <- function(x, ...){

  #number of vertices
  n_vertices <- igraph::vcount(
    x$Structure$visual_dag$graph
  )

  #number of unique synthesized parameters
  n_parameters <- sum(
    vapply(
      x$Parameters,
      length,
      integer(1)
    )
  )

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

    #method
    if(!is.null(x$Info$method)){
      cat("Method :", x$Info$method, "\n")
    }

    #scale
    if(!is.null(x$Info$scale)){
      cat("Scale  :", x$Info$scale, "\n\n")
    }


    ################
    # LOCAL UPDATE #
    ################

    if(
      !is.null(x$Info$scale) &&
      x$Info$scale == "local"
    ){

      lambda_local <- unlist(x$Info$lambda)
      ESS_local    <- unlist(x$Info$ESS)

      #make sure order is identical
      node_names <- intersect(
        names(lambda_local),
        names(ESS_local)
      )

      update_table <- data.frame(
        Vertex = node_names,
        Lambda = as.numeric(lambda_local[node_names]),
        ESS    = as.numeric(ESS_local[node_names]),
        row.names = NULL
      )

      update_table$Lambda <- round(
        update_table$Lambda,
        3
      )

      update_table$ESS <- round(
        update_table$ESS,
        1
      )

      print(
        update_table,
        row.names = FALSE
      )
    }


    #################
    # GLOBAL UPDATE #
    #################

    if(
      !is.null(x$Info$scale) &&
      x$Info$scale == "global"
    ){

      update_table <- data.frame(
        Lambda = round(
          as.numeric(x$Info$lambda)[1],
          3
        ),
        ESS = round(
          as.numeric(x$Info$ESS)[1],
          1
        )
      )

      print(
        update_table,
        row.names = FALSE
      )
    }
  }


  ###################
  # parameter table #
  ###################

  cat("\n")

  print(
    x$`Parameter table`,
    row.names = FALSE
  )

  invisible(x)
}
