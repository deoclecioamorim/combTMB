#' @title Simulate from a combTMB fitted model
#' @rdname simulate
#' @name simulate_combTMB
#' @param object combTMB fitted model
#' @param nsim number of response lists to simulate. Defaults to 1.
#' @param seed random number seed. Guarantees the reproducibility of results
#' @param \dots currently not used
#' @details In this function, random effects are simulated from their estimated
#' distribution. In the current version, it is not possible to condition the
#' estimated random effects.
#' @return A tibble. The list has length \code{nsim}.
#' Each simulated vector of observations is the same size as the vector of response
#' variables in the original data set.In the binomial family case each simulation
#' is a two-column matrix with success/failure.
#' @importFrom stats simulate
#' @export
simulate.combTMB <-function(object, nsim=1, seed = NULL, ...){
  #List of family allowed in simulate_combTMB function
  fnames <- .family_list_imple
  #Defensive programming
  if(!object$family$family %in% fnames) {
    cli_abort("Can't possible to generate simulated values for this family
              at the moment (sorry)")
  }

  #copied from stats:::simulate.lm
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  #Adapted from glmtmb and stats:::simulate.lm
  family <- object$family$family
  ret_sim <- replicate(nsim,
                   object$tmb_obj$simulate(par = object$fit$parfull)$Y,
                   simplify=FALSE)
  if (binomialType(family)) {
    size <- object$tmb_obj$env$data$size
    ret_sim <- lapply(ret_sim, function(x) cbind(x, size - x, deparse.level=0L))
    class(ret_sim) <- "data.frame"
    rownames(ret_sim) <- as.character(seq_len(nrow(ret_sim[[1L]])))
  } else {
    ret_sim <- as.data.frame(ret_sim)
  }

  names(ret_sim) <- paste0("sim_", seq_len(nsim))
  attr(ret_sim, "seed") <- RNGstate
  class(ret_sim) <- c("tbl_df", "tbl", "data.frame")
  ret_sim
}



