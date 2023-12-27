#'
#' Decision about dose escalation based on the toxicities (Bayesian Decision-Theoretic Design)
#'
#' @description  Find the intervals of biological efficacious doses (BEDs) based on prior information, current MTD, and current data. 
#'
#' @usage dose_escalation(d, a.pTox, b.pTox, target, alpha, eta)
#'
#' @param d the current dose combination
#' @param a.pTox parameter alpha for current prior distribution
#' @param b.pTox parameter beta for current prior distribution
#' @param target the target toxicity level
#' @param alpha penalty term if the current DLT probability is below the target DLT probability
#' @param eta penalty term if the current DLT probability is above the target DLT probability
#' 
#' @details \code{dose_escalation()} calculate the values of utility function, the dose combination with greatest utility is the next dose allocation
#'          In this function, if the current dose combination is (r,s), next dose allocation is limited to dose combination (i,j) such that i <= r and j<= s, as well as the combinations (r+1,s) and (r,s+1).
#'
#' @return \code{dose_escalation()} returns the indexes of next dose allocation
#'
#' @seealso Fan, S. K., Lee, B. L. and Lu, Y.: A Curve-free Bayesian Decision-theoretic Design for Phase Ia/Ib Trials Considering both Safety and Efficacy Outcomes. Statistics in Biosciences, 12(2), 146-166 (2020). \url{https://doi.org/10.1007/s12561-020-09272-5}
#'
#' @export

dose_escalation = function (d, a.pTox, b.pTox, target, alpha, eta) {
  
  candidates = c()
  
  for (i in 1:d[1]) {
    for (j in 1:d[2]) {
      candidates = rbind (candidates, c(i,j))
    }
  }
  
  if(d[1] + 1 <= nrow(a.pTox)) {candidates = rbind (candidates, c(d[1] + 1, d[2]))}
  if(d[2] + 1 <= ncol(a.pTox)) {candidates = rbind (candidates, c(d[1], d[2] + 1))}
  
  utilities = sapply(1:nrow(candidates), function(index) {-(alpha + eta)*(target*pbeta(target, a.pTox[candidates[index,1], candidates[index,2]], b.pTox[candidates[index,1], candidates[index,2]]) - a.pTox[candidates[index,1], candidates[index,2]]/(a.pTox[candidates[index,1], candidates[index,2]] + b.pTox[candidates[index,1], candidates[index,2]]) * pbeta(target, a.pTox[candidates[index,1], candidates[index,2]] + 1, b.pTox[candidates[index,1], candidates[index,2]])) - eta*(a.pTox[candidates[index,1], candidates[index,2]]/(a.pTox[candidates[index,1], candidates[index,2]] + b.pTox[candidates[index,1], candidates[index,2]]) - target)})
  return (candidates[which.max(utilities),])
}