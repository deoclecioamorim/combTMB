#' @title Extract Random Effects
#' @description Extract random effects from a fitted \code{combTMB} model
#' @param object a \code{combTMB} model
#' @param \dots Currently not used
#' @return
#' For \code{ranef}, an object of class \code{ranef.combTMB} with two components:
#' \describe{
#'   \item{cond}{a list of data frames, containing random effects}
#'     \item{fullranef}{a data frame with:}
#'    \describe{
#'       \item{group}{level of the grouping variable}
#'       \item{estimate}{value of the random-effects}
#'       \item{std.error}{standard deviation}
#'        }
#'}
#'
#' @note To show the complete structure use:
#'  \code{print(ranef(model),simplify=FALSE)}
#'
#' @examples
#' \donttest{
#' library(combTMB)
#' #Poisson-Normal model
#' m1 <- combTMB(OT ~ Period+(1|Donor), family=poisson(), data=embryos)
#' summary(m1)
#' r <- ranef(m1)
#' print(r, simplify=FALSE)
#' #Extract Donor
#' r$cond$Donor
#' }
#' @docType methods
#' @aliases ranef ranef.combTMB
#' @importFrom nlme ranef
#' @export ranef
#' @export
ranef.combTMB <- function(object, ...) {
  if (!missing(...)) {
    cli_alert_warning("Extra arguments discarded")
  }


  if(.has_random(object$formula)) {
    n_radom <- length(object$split_formula[[1]]$reTrmFormulas) #number of random
    out_ranef <- list()
    re_est <- as.list(object$tmb_report, "Estimate")$b
    re_ses <- as.list(object$tmb_report, "Std. Error")$b
    for(i in 1: n_radom) {
      level_names <- levels(object$frame[[object$split_formula[[1]]$reTrmFormulas[[i]][[3]]]])
      n_levels <- length(level_names)
      re_name <- object$split_formula[[1]]$reTrmFormulas[[i]][[3]]

      if(i==1) {
        start_pos <- 1
        end_pos <- n_levels
      } else {
        start_pos <- end_pos + 1
        end_pos <- start_pos + n_levels - 1
      }
      out_ranef[[i]] <- data.frame(
        group = paste0(re_name,"_", level_names),
        estimate = re_est[start_pos:end_pos],
        std.error = re_ses[start_pos:end_pos],
        stringsAsFactors = FALSE
      )

    }
    out_ranef <- do.call("rbind", out_ranef)
    row.names(out_ranef) <- NULL
  }

  .ranef <- out_ranef

  group <- unlist(lapply(strsplit(.ranef$group,"_"), getElement, 1))
  est <- .ranef$estimate

  cond <- list()
  for(j in 1:length(unique(group))) {
    cond[[unique(group)[j]]] = data.frame("Intercept" = est[which(group == unique(group)[j])])
  }

  structure(list(cond = cond, fullranef= out_ranef),
            class = "ranef.combTMB")


}


#' @method print ranef.combTMB
#' @export
print.ranef.combTMB <- function(x,  simplify=TRUE, ...) {
  print(if (simplify)
    unclass(x$cond) else unclass(x),
    ...)
  invisible(x)
}




