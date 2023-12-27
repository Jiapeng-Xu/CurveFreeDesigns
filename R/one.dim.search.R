#'
#' A recursion function that perform one dimensional test of the Rule-Based Design for Agents with Same Dose-Limiting Toxicities at dose combination (i,j)
#'
#' @description A recursion function that perform one dimensional test of the Rule-Based Design for Agents with Same Dose-Limiting Toxicities at dose combination (i,j), update the two-dimension isotonic estimation, and return the results
#'
#' @usage one.dim.search(p_hat, n1.dose, n2.dose, n, p, t, r, i, j)
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
#' @return \code{one.dim.search()} returns:
#'         (1) \code{$t}: a matrix containing the number of DLT outcomes of each dose level
#'         (2) \code{$n}: a matrix containing the total number of patients of each dose level
#'         (3) \code{$r}: a matrix containing the rejection decision of each dose level
#'
#' @seealso Fan SK, Venook AP, Lu Y. Design issues in dose-finding Phase I trials for combinations of two agents. J Biopharm Stat. 2009;19(3):509-23. doi: 10.1080/10543400902802433. PMID: 19384692.
#' 
#' @export

one.dim.search = function(p_hat, n1.dose, n2.dose, n, p, t, r, i, j) {
  rowpt=0; colpt=0
  if (i+1<= n1.dose && r[i+1,j]!=1) rowpt=1;    
  if (j+1<= n2.dose && r[i,j+1]!=1) colpt=1;
  
  # if dose combination (i,j) hasn't been tested yet
  if(r[i,j] == 0) {
    res = test(n1.dose, n2.dose,n, p, t, r, i, j)
    t = res$t
    n = res$n
    r = res$r
  }
  
  # if dose combination (i,j) is rejected
  if (r[i,j] == 1) {
    if (j-1 >= 1 && r[i, j-1] != 1) {
      res = one.dim.search (p_hat, n1.dose, n2.dose, n, p, t, r, i, j-1)
      t = res$t
      n = res$n
      r = res$r
    }
    if (i-1 >= 1 && r[i-1, j] != 1) {
      res = one.dim.search (p_hat, n1.dose, n2.dose, n, p, t, r, i-1, j)
      t = res$t
      n = res$n
      r = res$r
    }
  }
  # if dose combination (i,j) is accepted
  else {
    # if can move up in both directions
    if (rowpt==1 & colpt==1) {
      # with lower dose in both directions
      if (j-1 >=1 & i-1 >=1) {
        p_hat = twodim.iso(n1.dose, n2.dose, p_hat, t, n)
          #twodim.iso(n1.dose, n2.dose, p_hat, t, n)
          #twodim.iso(p_hat, t, n, 4)
        prow = p_hat[i,j]-p_hat[i-1,j]
        pcol = p_hat[i,j]-p_hat[i,j-1]
        # having data to compare
        if (n[i-1,j]>0 & n[i,j-1]>0){
          if (prow > pcol) {
            res = one.dim.search(p_hat, n1.dose, n2.dose, n, p, t, r, i+1, j)
            t = res$t
            n = res$n
            r = res$r
          }
          if (prow < pcol) {
            res = one.dim.search(p_hat, n1.dose, n2.dose, n, p, t, r, i, j+1)
            t = res$t
            n = res$n
            r = res$r
          }
          if (prow == pcol) {
            res = two.dim.search(p_hat, n1.dose, n2.dose, n, p, t, r, i, j)
            t = res$t
            n = res$n
            r = res$r
          }
        }
      }
      else {
        res = two.dim.search(p_hat, n1.dose, n2.dose, n, p, t, r, i, j)
        t = res$t
        n = res$n
        r = res$r
      }
    }
    # only row moving up is possible
    if (rowpt==1 & colpt==0){
      res = one.dim.search(p_hat, n1.dose, n2.dose, n, p, t, r, i+1, j)
      t = res$t
      n = res$n
      r = res$r
    }
    # only row moving up is possible
    if (rowpt==0 & colpt==1){
      res = one.dim.search(p_hat, n1.dose, n2.dose, n, p, t, r, i, j+1)
      t = res$t
      n = res$n
      r = res$r
    }     
  }
  return(list(t = t,
              n = n,
              r = r))
}