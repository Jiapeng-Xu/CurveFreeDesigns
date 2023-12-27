#'
#' Terminate the phase Ia trial whenever the evidence for the current conclusion is considered sufficient
#'
#' @description  Stop if (1) out of patients n = n_{max} or; (2) the current recommended dose is very likely to be the MTD or; (3) all doses are evidently too toxic. Use working data (see function workdata) instead of real data for small samples
#'
#' @usage stopmtd(a.vec, b.vec, mtd, target, T.max, p1, p2)
#'
#' @param a.vec parameter alpha for current prior distribution
#' @param b.vec parameter beta for current prior distribution
#' @param mtd the MTD candidate
#' @param target the target toxicity level
#' @param T.max the absolute overly toxic level, usually .05 above target
#' @param p1 the error rate for concluding inadmissible; .10 is suggested
#' @param p2 the error rate for concluding admissible; .10 is suggested
#' 
#' @details \code{stopmtd()} performs the following stopping rule:
#'          (1) stop, if the current recommended dose is very likely to be the MTD:
#'              (a) if the current recommended dose is the highest dose: if the posterior probability of the current MTD being the true MTD exceeds a threshold, we conclude that the current MTD is the true MTD
#'              (b) if the current recommended dose is not the highest dose: if the posterior probability of the next higher dose being too toxic exceeds a threshold, we conclude that the current MTD is the true MTD
#'          (2) stop, if all doses are evidently too toxic:
#'              (a) if the posterior probability of the first dose being overly toxic exceeds a threshold, we conclude that the MTD is below the first dose
#'
#' @return \code{stopmtd()} returns (1) flag_stop: Boolean value indicating whether we should terminate the phase Ia trial and; (2) flag_mtd: Boolean value indicating whether we have sufficient evidence that the current MTD is the true MTD
#'
#' @seealso Fan, S., Lee, B. L., & Lu, Y. (2020). A curve-free bayesian decision-theoretic design for phase Ia/Ib trials considering both safety and efficacy outcomes. \emph{Statistics in Biosciences}, 12(2), 146–166. \url{https://doi.org/10.1007/s12561-020-09272-5}
#'          Fan SK, Lu Y, Wang YG. A simple Bayesian decision-theoretic design for dose-finding trials. \emph{Stat Med}, 2012 Dec 10;31(28):3719-30. doi: 10.1002/sim.5438. Epub 2012 Jul 5. PMID: 22763943. \url{https://doi.org/10.1002/sim.5438}
#'
#' @export

stopmtd <- function(a.vec, b.vec, mtd, target, T.max, p1, p2) {

  p.lowest <- pbeta(T.max, a.vec[1], b.vec[1])
  if (p.lowest < p1) return(2) # stop the trial and all doses likely to be overly toxic
  
  if (mtd == length(a.vec)) { # mtd is the highest dose
    if (pbeta(target, a.vec[mtd], b.vec[mtd]) > 1-p2) 
      return(1) # stop the trial and highest is highly likely to be safe
  }
  else{ # mtd is not the highest dose
    p.mtd.plus <- pbeta(T.max, a.vec[mtd + 1], b.vec[mtd + 1])
    if (p.mtd.plus < p2) return(1) # next higher dose likely to be overly toxic
  }
  return(0) # continue the trial 
}

