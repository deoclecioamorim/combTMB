#' Extract MCMC samples from an combTMB model fit with tmbstan.
#'
#' @param object Output from [tmbstan::tmbstan()] run on the `tmb_obj`
#'   element of an [combTMB::combTMB()] model. E.g.,
#'   `tmbstan::tmbstan(your_model$tmb_obj)`.
#'
#' @return
#' Returns a matrix of parameter samples. Rows correspond to the order
#' of `your_model$tmb_obj$env$last.par.best`. Columns correspond to
#' posterior samples.
#'
#' @importFrom stats predict
#'
#' @export
extract_mcmc <- function(object) {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    cli::cli_abort("rstan must be installed to use `extract_mcmc()`.")
  }
  post <- rstan::extract(object)
  p_names <- names(post)[-length(names(post))] # exclude "lp__"
  p <- lapply(seq_len(length(post[["lp__"]])), function(i) {
    post_pars <- list()
    for (j in seq_along(p_names)) {
      par_j <- p_names[j]
      if (is.matrix(post[[par_j]])) {
        post_pars[[j]] <- post[[par_j]][i, , drop = TRUE]
      } else {
        post_pars[[j]] <- post[[par_j]][i]
      }
    }
    post_pars
  })
  simplify2array(lapply(p, unlist))
}
