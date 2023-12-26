#'
#' Find the maximum tolerated dose (MTD)
#'
#' @description Find the maximum tolerated dose (MTD) based on prior information and current data. Use working data instead of real data for small samples
#'
#' @usage findmtd(target, a.pTox, b.pTox, alpha, eta)
#' 
#' @param target the target toxicity level
#' @param a.pTox parameter alpha for current prior distribution
#' @param b.pTox parameter beta for current prior distribution
#' @param alpha the weight of penalty of toxicity rate below target level; 1 is suggested
#' @param eta the weight of penalty of toxicity rate above target level; 1 is suggested
#'
#' @return \code{findmtd()} returns: (1) gain and (2) mtd
#'
#' @seealso Fan, S., Lee, B. L., & Lu, Y. (2020). A curve-free bayesian decision-theoretic design for phase Ia/Ib trials considering both safety and efficacy outcomes. \emph{Statistics in Biosciences}, 12(2), 146–166. \url{https://doi.org/10.1007/s12561-020-09272-5}
#'          Fan SK, Lu Y, Wang YG. A simple Bayesian decision-theoretic design for dose-finding trials. \emph{Stat Med}, 2012 Dec 10;31(28):3719-30. doi: 10.1002/sim.5438. Epub 2012 Jul 5. PMID: 22763943. \url{https://doi.org/10.1002/sim.5438}
#'
#' @export

findmtd <- function(target, a.pTox, b.pTox, alpha, eta){

  # penalty function for mtd: E(alpha*|p-pi.T|_ +eta*|p-pi.T|+)

  mu <- a.pTox / (a.pTox + b.pTox) # a list of prior mean toxicity rate
  gain <- target * pbeta(target, a.pTox, b.pTox) - mu * pbeta(target, a.pTox + 1, b.pTox)
  gain <- -(alpha + eta) * gain - eta * (mu - target)
  mtd <- which.max(gain)
  
  return(list(gain = gain, mtd = which.max(gain)))
}
