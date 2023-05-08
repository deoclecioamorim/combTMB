#' @title Extract residual standard deviation or dispersion parameter
#' @param object an object of class \dQuote{combTMB}
#' @param \dots currently not used.
#' @return .
#'
#' @aliases dispersion dispersion.combTMB
#' @export dispersion
dispersion <- function(object, ...) UseMethod("dispersion")
#' @export
#'
dispersion.combTMB <- function(object, ...) {
  pl <- .getParList(object)
  ff <- object$family$family
  q <- object$tmb_params$betaq #parameter bbn used with map q=0
  if (!usesDispersion(ff)) return(1.)
  if (length(pl$betad)>1) return(NA)
  switch(family(object)$family, gaussian = exp(0.5 * pl$betad),
         betabinomial = 1/(1+exp(pl$betad)),
         gnb=if(q==0){exp(pl$betad)}else{1/exp(pl$betad)},
         exp(pl$betad))
}

#Colocar a descrição das extrações do sigma

