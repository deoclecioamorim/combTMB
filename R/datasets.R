#-----------------------------------------------------------------------
#' @title A clinical trial in toenail infection named onychomycosis
#'
#' @description This dataset comes from a randomized, multicenter, parallel-group
#' study to compare two oral treatments (coded as 0 and 1) for toenail dermatophytic
#' onychomycosis. For more details see Molenberghs et al. (2010).
#'
#' @format A data frame with 1908 rows and four variables:
#' \describe{
#' \item{idnum}{number of patients}
#' \item{time}{number of measurements per subject (in time)}
#' \item{treatn}{two oral treatments (0 and 1)}
#' \item{onyresp}{response variable: severity of infection, coded as 0
#' (not severe) or 1 (severe)}
#' }
#' @usage data(toenail, package = "combTMB")
#'
#' @references Molenberghs, G., Verbeke, G., Demetrio, C. G., and Vieira, A. M. (2010).
#' A family of generalized linear models for repeated measures with normal
#' and conjugate random effects. Statistical science, 25(3), 325-347.
#'
"toenail"


#-----------------------------------------------------------------------
#' @title Data on embryo production in Holstein cows
#'
#' @description In vitro embryo production (IVP) is the reproductive biotechnology
#'  with the greatest impact on the genetic improvement of cattle, according to
#'  the International Embryo Technology Society - IETS (Viana, 2020).
#'  This dataset describes the main variables used to measure quality/quantity
#'  of oocytes and embryos in Holstein cows. The composition of the database was
#'  obtained from follicular aspiration sessions, performed between January 6
#'  and December 28, 2020.
#' @format A data frame with 1148 rows and 15 variables:
#' \describe{
#' \item{Period}{periods of the year in which the follicular aspiration sessions were performed (P1 and P2). P1 (period of highest temperature comprising the months of June, July, August, September and October) and in period P2 (period of lowest temperature comprising the remaining seven months)}
#' \item{Interval}{interval between ovum pick-up performed in the same donor (2wks - two, 3wks - three and na - more than three weeks)}
#' \item{Injections}{number of injections of follicle stimulating hormone - FSH (0, 1 and 5),}
#' \item{Pregnant}{0 = empty and 1 = pregnant}
#' \item{Status}{H-heifers, M-lactating cows and D-dry cows}
#' \item{Donor}{variable referring to donor identification}
#' \item{O12}{number of grade one and two oocytes}
#' \item{OM}{number of oocytes from grade one to grade four}
#' \item{OT}{total number of oocytes}
#' \item{Sire}{identification of the semen donor bull}
#' \item{sexed}{sexed semen: 1 = yes and 0 = no}
#' \item{IVC}{number of oocytes in vitro culture}
#' \item{D3}{number of cleaved zygotes}
#' \item{Embryos}{number of embryos}
#' \item{id}{number of subjects}
#' }
#' @usage data(embryos, package = "combTMB")
#'
#' @references Viana, J. (2020). 2019 statistics of embryo production and transfer in domestic
#'farm animals, Embryo Technology Newsletter 38(4): 7–26

"embryos"

