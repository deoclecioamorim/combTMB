#' @title Extract or obtain generalized components from a glm or glmm model or combined model
#' @name extractors
#' @param object an object of class \dQuote{combTMB}
#' @param name of the component to be retrieved
#' @param \dots currently not used.
#' @return .
#'
NULL

#' @aliases getME
#' @rdname extractors
#' @seealso \code{\link[lme4]{getME}}
#' @importFrom lme4 getME
#' @export getME
#' @method getME combTMB
#' @export
#'
getME.combTMB <- function(object,
                              name = c("X","Xd","Z", "beta","betad" ,"b"),
                              ...) {
  if(missing(name)) stop("'name' must not be missing")
  ## Deal with multiple names -- "FIXME" is inefficiently redoing things
  if (length(name <- as.character(name)) > 1) {
    names(name) <- name
    return(lapply(name, getME, object = object))
  }
  if(name == "ALL") {
    return(sapply(eval(formals()$name),
                  getME.combTMB, object=object, simplify=FALSE))
  }

  stopifnot(inherits(object, "combTMB"))
  name <- match.arg(name)

  oo.env <- object$tmb_obj$env
  ### Start of the switch
  allpars <- oo.env$parList(object$fit$par, object$fit$parfull)
  switch(name,
         "X"     =  oo.env$data$X,
         "Xd"     =  oo.env$data$Xd,
         "Z"     = oo.env$data$Z,
         "beta"  = unlist(allpars[c("beta")]),
         "betad"  = unlist(allpars[c("betad")]),
         "b"  = unlist(allpars[c("b")]),
         "..foo.." = # placeholder!
           stop(gettextf("'%s' is not implemented yet",
                         sprintf("getME(*, \"%s\")", name))),
         ## otherwise
         stop(sprintf("Mixed-Effects extraction of '%s' is not available for class \"%s\"",
                      name, class(object))))

}






