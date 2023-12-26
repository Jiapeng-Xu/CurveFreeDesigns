#'
#' Determine whether (i,j) is the largest acceptable one
#'
#' @description Determine whether (i,j) is the largest acceptable one
#'
#' @usage low(p_hat, n1.dose, n2.dose, mtdp, i, j)
#'
#' @param p_hat a matrix containing the isotonic estimation of DLT probability of each dose level
#' @param n1.dose the number of dose level of agent 1
#' @param n2.dose the number of dose level of agent 2
#' @param mtdp the toxicity rate of MTD
#' @param i the current dose level of agent 1
#' @param j the current dose level of agent 2
#'
#' @return \code{low()} returns: 1 = no, 0 = yes
#'
#' @seealso Fan SK, Venook AP, Lu Y. Design issues in dose-finding Phase I trials for combinations of two agents. J Biopharm Stat. 2009;19(3):509-23. doi: 10.1080/10543400902802433. PMID: 19384692.
#' 
#' @export

low = function(p_hat, n1.dose, n2.dose, mtdp, i, j){
    answer=0
    for (k in i:n1.dose)
      for (l in j:n2.dose){
        if (p_hat[k,l]==mtdp & k+l>i+j)
          answer=1
      }
    return (answer)
  }