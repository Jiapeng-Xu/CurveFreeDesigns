#'
#' A recursion function that perform two dimensional test of the Rule-Based Design for Agents with Same Dose-Limiting Toxicities at dose combination (i,j)
#'
#' @description A recursion function that perform two dimensional test of the Rule-Based Design for Agents with Same Dose-Limiting Toxicities at dose combination (i,j), update the two-dimension isotonic estimation, and return the results
#'
#' @usage two.dim.search(p_hat, n1.dose, n2.dose, n, p, t, r, i, j)
#'
#' @param p_hat a matrix containing the isotonic estimation of DLT probability of each dose level
#' @param n1.dose the number of dose level of agent 1
#' @param n2.dose the number of dose level of agent 2
#' @param n a matrix containing the total number of patients of each dose level
#' @param p a matrix containing the true DLT probability of each dose level
#' @param t a matrix containing the number of DLT outcomes of each dose level
#' @param r a matrix containing the rejection decision of each dose level
#' @param i the current dose level of agent 1
#' @param j the current dose level of agent 2
#'
#' @return \code{two.dim.search()} returns:
#'         (1) \code{$t}: a matrix containing the number of DLT outcomes of each dose level
#'         (2) \code{$n}: a matrix containing the total number of patients of each dose level
#'         (3) \code{$r}: a matrix containing the rejection decision of each dose level
#'
#' @seealso Fan SK, Venook AP, Lu Y. Design issues in dose-finding Phase I trials for combinations of two agents. J Biopharm Stat. 2009;19(3):509-23. doi: 10.1080/10543400902802433. PMID: 19384692.
#' 
#' @export

two.dim.search = function(p_hat, n1.dose, n2.dose, n, p, t, r, i, j) {
  
  if (r[i+1,j]==0) {
    res = test(n1.dose, n2.dose,n, p, t, r, i+1, j)
    t = res$t
    n = res$n
    r = res$r
  } 
  if (r[i,j+1]==0) {
    res = test(n1.dose, n2.dose,n, p, t, r, i, j+1)
    t = res$t
    n = res$n
    r = res$r
  }
  # if both dose combination are accepted
  if (r[i+1,j]==2 & r[i,j+1]==2) {
    res = one.dim.search(p_hat, n1.dose, n2.dose, n, p, t, r, i+1,j+1)
    t = res$t
    n = res$n
    r = res$r
  }
  
  # if dose (i+1, j) is rejected
  if (r[i+1,j]==1) {
    if (j-1>=1) {
      res = one.dim.search(p_hat, n1.dose, n2.dose, n, p, t, r, i+1,j-1)
      t = res$t
      n = res$n
      r = res$r
    }
  }
  else {
    if (i+2<= n1.dose) {
      res = one.dim.search(p_hat, n1.dose, n2.dose, n, p, t, r, i+2,j)
      t = res$t
      n = res$n
      r = res$r
    }
  }
  # if dose (i, j+1) is rejected
  if (r[i,j+1]==1) {
    if (i-1>=1) {
      res = one.dim.search(p_hat, n1.dose, n2.dose, n, p, t, r, i-1,j+1)
      t = res$t
      n = res$n
      r = res$r
    }
  }
  else {
    if (j+2<=n2.dose) {
      res = one.dim.search(p_hat, n1.dose, n2.dose, n, p, t, r, i,j+2)
      t = res$t
      n = res$n
      r = res$r
    }
  }
  return(list(t = t,
              n = n,
              r = r))
}