#' Helper function to create parameter table
#'
#' @param param_est Parameter estimate
#' @param e_fun Edge function
#'
.ppmn_tab_fun <- function(param_est, e_fun){

  #dependent variables from edge names
  dependent_vars <- sub("~.*", "", names(param_est))

  #build table per edge
  par_tab <- lapply(seq_along(param_est), function(i){

    sub_par <- param_est[[i]]
    edge    <- names(param_est)[i]

    #independent variables from edge
    indep <- strsplit(
      strsplit(edge, "~", fixed = TRUE)[[1]][2],
      "+",
      fixed = TRUE
    )[[1]]

    #match independent variables to parameters
    if(e_fun[[i]] %in%
       c("asymptotic", "sigmoidal", "mvlogistic",
         "gompertz", "gaussian")){

      indep <- c(NA_character_, NA_character_, indep)

    }else if(e_fun[[i]] == "class"){

      indep <- c(
        NA_character_,
        indep,
        NA_character_,
        indep,
        NA_character_
      )

    }else{

      indep <- c(NA_character_, indep)
    }

    #parameter estimates
    estimates <- as.data.frame(
      do.call(rbind, sub_par)
    )

    #check dimensions
    if(length(indep) != nrow(estimates)){
      stop(
        "Number of parameters and independent-variable labels ",
        "do not match for edge: ", edge
      )
    }

    #combine
    data.frame(
      edge = edge,
      `edge function` = e_fun[[i]],
      dependent = dependent_vars[i],
      independent = indep,
      parameters = names(sub_par),
      estimates,
      check.names = FALSE,
      row.names = NULL
    )
  })

  #combine edges
  full_table <- do.call(
    rbind,
    par_tab
  )

  rownames(full_table) <- NULL

  return(full_table)}
