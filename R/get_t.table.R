#'
#' Return a table containing the number of DLT1, DLT2, and DLT3 based on the conditional probability specified by users
#'
#' @description Return a table containing the number of DLT1, DLT2, and DLT3 based on the conditional probability specified by users
#'
#' @usage get_t.table(d, n_DLTs, p.conditional1, p.conditional2)
#'
#' @param d the current dose combination
#' @param n_DLTs the number of DLT outcomes
#' @param p.conditional1 The conditional probability of observing only DLT1, given that at least one DLT occurs and that DLT3 is absent
#' @param p.conditional2 The conditional probability of observing only DLT2, given that at least one DLT occurs and that DLT3 is absent
#'
#' @return \code{get_decision()} returns the result as a list,
#'         including:
#'         (1) \code{$t.table}: A matrix containing the number of DLT1, DLT2, and DLT3, respectively.
#'
#' @seealso Lee BL, Fan SK. A two-dimensional search algorithm for dose-finding trials of two agents. J Biopharm Stat. 2012;22(4):802-18. doi: 10.1080/10543406.2012.676587. PMID: 22651116.
#'
#' @export

get_t.table = function(d, n_DLTs, p.conditional1, p.conditional2) {
  
  if (n_DLTs == 0) {
    return (NULL)
  }
  
  # within the n DLTs, the times that at least one DLT3 is observed
  DLT3 = sum(runif(n_DLTs) < 0.2)
  
  if (DLT3 > 0) {
    
    # the four rows represent "DLT1", "DLT2", "DLT1 & DLT2", "DLT3 only" respectively
    res.DLT3 = rmultinom(DLT3, 1, c(0.25, 0.25, 0.25, 0.25))
    
    t.table.fun1 = function (col) {
      if (col[1] == 1) {
        res = c(1,0,1)
      }
      else if (col[2] == 1) {
        res = c(0,1,1)
      }
      else if (col[3] == 1) {
        res = c(1,1,1)
      }
      else {
        res = c(0,0,1)
      }
    }
    t.table1 = apply(res.DLT3,2,t.table.fun1)
    
  }
  else t.table1 = NULL
  
  
  # within the n DLTs, the times that no DLT3 is observed
  no.DLT3 = n_DLTs - DLT3
  
  if (no.DLT3 > 0) {
    
    # the three rows represent "DLT1 only", "DLT2 only", "DLT1 & DLT2" respectively
    res.no.DLT3 = rmultinom(no.DLT3, 1, c(p.conditional1[d[1], d[2]],
                                          p.conditional2[d[1], d[2]], 
                                          1-p.conditional1[d[1], d[2]]-p.conditional2[d[1], d[2]]))
    
    t.table.fun2 = function (col) {
      if (col[1] == 1) {
        res = c(1,0,0)
      }
      else if (col[2] == 1) {
        res = c(0,1,0)
      }
      else {
        res = c(1,1,0)
      }
    }
    t.table2 = apply(res.no.DLT3,2,t.table.fun2)
    
  }
  else t.table2 = NULL
  
  
  t.table = cbind(t.table1, t.table2)
  rownames(t.table) = c("T1", "T2", "T3")
  colnames(t.table) = 1:n_DLTs
  
  return (t.table)
}
