#'
#' Find intervals of biological efficacious doses (BEDs) for CFHD
#'
#' @description  Find the intervals of biological efficacious doses (BEDs) based on prior information, current MTD, and current data. 
#'
#' @usage findbeds_CFHD(mtd, n.eff, n.assign, a.pEff, b.pEff, E.min, gain.A, gain.AC, phi, lo)
#'
#' @param mtd the current maximum tolerated dose (MTD)
#' @param n.eff a list of no. of efficacies at all dose levels
#' @param n.assign a list of no. of patients at all dose levels
#' @param a.pEff parameter alpha for current prior distribution
#' @param b.pEff parameter beta for current prior distribution
#' @param E.min the minimum acceptable efficacy rate
#' @param gain.A the weight of reward that a BED is selected correctly; 1 is suggested
#' @param gain.AC the weight of reward that a non-BED is not selected correctly; 1 is suggested
#' @param phi the weight for how much more efficacious P_{E,i} is than E.min; 1 is suggested
#' @param lo the weight for how much less efficacious P_{E,i} is than E.min; 1 is suggested
#' 
#' @details \code{findbeds()} implement the following procedure:
#'          (1) find the interval with maximum utility A among all the possible intervals
#'          (2) find the subinterval B that maximizes the posterior probability
#'          In this function, only doses at or below MTD are investigated and we use real data rather than working data to address the BED interval.
#'
#' @return \code{findbeds()} returns (1) iL.B: the lower boundary of the interval of BEDs and; (2) iU.B: the upper boundary of the interval of BEDs and; (3) p.best: the posterior probability that interval B is acceptable
#'
#' @seealso Fan, S., Lee, B. L., & Lu, Y. (2020). A curve-free bayesian decision-theoretic design for phase Ia/Ib trials considering both safety and efficacy outcomes. \emph{Statistics in Biosciences}, 12(2), 146–166. \url{https://doi.org/10.1007/s12561-020-09272-5}
#'
#' @export

findbeds_CFHD <- function(mtd, n.eff, n.assign, a.pEff, b.pEff, E.min, gain.A, gain.AC, phi, lo) {
  max.Lik = - .Machine$double.xmax
  iL = 0
  iU = 0
  for (i in 1:mtd) {
    for (j in i:mtd) {
      Lik = 0 # reset likelihood for each potential interval
      for (k in 1:mtd) { # traverse all dose level
        if (n.assign[k] > 0) {
          if (i <= k & k <= j) { # if dose i belongs to interval A
            Lik = Lik + n.assign[k] * (ifelse(n.eff[k]/n.assign[k] >= E.min, 1, 0)*(gain.A + phi*(n.eff[k]/n.assign[k] - E.min)))
          }
          else { # if dose i does not belong to interval A
            Lik = Lik + n.assign[k] * (ifelse(n.eff[k]/n.assign[k] <  E.min, 1, 0)*(gain.AC + lo*(n.eff[k]/n.assign[k] - E.min)))
          }
        }
      }
      
      # if (Lik > max.Lik) {
      #   max.Lik = Lik
      #   iL = i
      #   iU = j
      # }
      if (Lik > max.Lik) {

        # Check whether the current BEDs interval is acceptable
        if (i == 1) {
          if (j == mtd) {
            if (sum(n.assign[i:j]) > 0) {
              if (sum(n.eff[i:j])/sum(n.assign[i:j]) >= E.min) {
                max.Lik = Lik
                iL = i
                iU = j}
            }
          }
          else {
            if (sum(n.assign[i:j]) > 0 & sum(n.assign[(j+1):mtd]) > 0) {
              if (sum(n.eff[(j+1):mtd])/sum(n.assign[(j+1):mtd]) < E.min &
                  sum(n.eff[i:j])/sum(n.assign[i:j]) >= E.min) {
                max.Lik = Lik
                iL = i
                iU = j}
            }
          }
        }
        else {
          if (j == mtd) {
            if (sum(n.assign[i:j]) > 0 & sum(n.assign[1:(i-1)]) > 0) {
              if (sum(n.eff[1:(i-1)])/sum(n.assign[1:(i-1)]) < E.min &
                  sum(n.eff[i:j])/sum(n.assign[i:j]) >= E.min) {
                max.Lik = Lik
                iL = i
                iU = j}
            }
          }
          else {
            if (sum(n.assign[i:j]) > 0 & sum(n.assign[1:(i-1)]) > 0 & sum(n.assign[(j+1):mtd]) > 0) {
              if (sum(n.eff[1:(i-1)])/sum(n.assign[1:(i-1)]) < E.min &
                  sum(n.eff[(j+1):mtd])/sum(n.assign[(j+1):mtd]) < E.min &
                  sum(n.eff[i:j])/sum(n.assign[i:j]) >= E.min) {
                max.Lik = Lik
                iL = i
                iU = j}
            }
          }
        }
      }
    }
  }
  return(c(iL, iU))
}

