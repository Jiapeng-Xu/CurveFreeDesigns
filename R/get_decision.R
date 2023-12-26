#'
#' Perform A+B+C test for The Rule-Based Design for Agents with Non-Overlapping Dose-Limiting Toxicities
#'
#' @description Perform A+B+C test for The Rule-Based Design for Agents with Non-Overlapping Dose-Limiting Toxicities
#'
#' @usage get_decision(d, p.true, p.conditional1, p.conditional2, A, B, C, a.E, a.D, b.E, b.D, c.E)
#'
#' @param d the current dose combination
#' @param p.true a list of true DLT rates at all dose combination levels
#' @param p.conditional1 The conditional probability of observing only DLT1, given that at least one DLT occurs and that DLT3 is absent
#' @param p.conditional2 The conditional probability of observing only DLT2, given that at least one DLT occurs and that DLT3 is absent
#' @param A The sample size of first cohort of A+B+C design
#' @param B The sample size of first cohort of A+B+C design
#' @param C The sample size of first cohort of A+B+C design
#' @param a.E is the maximum number of DLTs needed out of A patients to to conclude the current dose is accepted;
#' @param a.D is the minimum number of DLTs needed out of A patients to to conclude the current dose is rejected;
#' @param b.E is the maximum number of DLTs needed out of B patients to to conclude the current dose is accepted;
#' @param b.D is the minimum number of DLTs needed out of B patients to to conclude the current dose is rejected;
#' @param c.E is the minimum number of DLTs needed out of C patients to to conclude the current dose is rejected, otherwise, the current dose will be accepted;
#'
#' @return \code{get_decision()} returns the result as a list,
#'         including:
#'         (1) \code{$flag.accepted}: A Boolean value indicating whether the current dose combination should be accepted
#'         (2) \code{$t.table}: A matrix containing the number of DLT1, DLT2, and DLT3, respectively.
#'         (3) \code{$n.DLTs}: The number of DLTs happens during the test
#'         (4) \code{$n.patients}: The number of patients enrolled for the test
#'
#' @seealso Lee BL, Fan SK. A two-dimensional search algorithm for dose-finding trials of two agents. J Biopharm Stat. 2012;22(4):802-18. doi: 10.1080/10543406.2012.676587. PMID: 22651116.
#'
#' @export

get_decision <- function (d, p.true, p.conditional1, p.conditional2, A, B, C, a.E, a.D, b.E, b.D, c.E) {
  
  # Enroll A patients, t.A DLTs are observed
  t.A = sum(runif(A) < p.true[d[1], d[2]])
  
  t.table = get_t.table(d, t.A, p.conditional1, p.conditional2)
  n.DLTs = t.A
  n.patients = A
  
  if (t.A <= a.E) {
    flag.accepted = TRUE
  }
  else if (a.E < t.A & t.A < a.D) {
    
    t.B = sum(runif(B) < p.true[d[1], d[2]])
    
    t.table = cbind(t.table, get_t.table(d, t.B, p.conditional1, p.conditional2))
    n.DLTs = n.DLTs + t.B
    n.patients = n.patients + B
    
    if (t.B <= b.E) {
      flag.accepted = TRUE
    }
    else if (b.E < t.B & t.B < b.D) {
      
      t.C = sum(runif(C) < p.true[d[1], d[2]])
      
      t.table = cbind(t.table, get_t.table(d, t.C, p.conditional1, p.conditional2))
      n.DLTs = n.DLTs + t.C
      n.patients = n.patients + C
      
      if (t.C <= c.E) {
        flag.accepted = TRUE
      }
      else {
        flag.accepted = FALSE
      }
    }
    else if (t.B >= b.D) {
      flag.accepted = FALSE
    }
  }
  else if (t.A >= a.D) {
    flag.accepted = FALSE
  }
  
  if (!is.null(t.table)) colnames(t.table) = 1:n.DLTs
  
  return (list(flag.accepted = flag.accepted,
               t.table = t.table,
               n.DLTs = n.DLTs,
               n.patients = n.patients))
}
