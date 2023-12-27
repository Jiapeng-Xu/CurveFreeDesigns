#'
#' Perform 2+1+3 test for the Rule-Based Design for Agents with Same Dose-Limiting Toxicities at dose combination (k,l), and return the results
#'
#' @description Perform 2+1+3 test for the Rule-Based Design for Agents with Same Dose-Limiting Toxicities at dose combination (k,l), and return the results
#'
#' @usage test(n1.dose, n2.dose, n, p, t, r, k, l)
#'
#' @param n1.dose the number of dose level of agent 1
#' @param n2.dose the number of dose level of agent 2
#' @param n a matrix containing the total number of patients of each dose level
#' @param p a matrix containing the true DLT probability of each dose level
#' @param t a matrix containing the number of DLT outcomes of each dose level
#' @param r a matrix containing the rejection decision of each dose level
#' @param k the current dose level of agent 1
#' @param l the current dose level of agent 2
#'
#' @return \code{test()} returns:
#'         (1) \code{$t}: a matrix containing the number of DLT outcomes of each dose level
#'         (2) \code{$n}: a matrix containing the total number of patients of each dose level
#'         (3) \code{$r}: a matrix containing the rejection decision of each dose level
#'
#' @seealso Fan SK, Venook AP, Lu Y. Design issues in dose-finding Phase I trials for combinations of two agents. J Biopharm Stat. 2009;19(3):509-23. doi: 10.1080/10543400902802433. PMID: 19384692.
#' 
#' @export

test = function(n1.dose, n2.dose, n, p, t, r, k, l) {

  # enroll 2 patients
  n[k,l]=2
  
  obs1 = runif(1)
  obs2 = runif(1)
  if (obs1<p[k,l]) t[k,l]=t[k,l]+1
  if (obs2<p[k,l]) t[k,l]=t[k,l]+1
  
  if (t[k,l]==2){ 
    rst=1 # (k,l) is rejected
  }
  else{
    if (t[k,l]==1){
      n[k,l] = n[k,l] + 1 # add 1 person
      obs3 = runif(1)
      if (obs3<p[k,l]){
        t[k,l]=t[k,l]+1
        rst=1
        } # (k,l) is rejected
      else rst=2
    }
    else{
      rst=2
    }
  }
  
  # update status for all 9 cells
  if (rst==1) {
    # If the current dose combination is rejected, we reject all the doses above the current dose
    r[k:n1.dose, l:n2.dose] = 1
  }
  else {
    r[1:k, 1:l] = 2
    }

  return (list(t = t,
               n = n,
               r = r))
}