#'
#' Generate operating characteristics of The Rule-Based Design for Agents with Non-Overlapping Dose-Limiting Toxicities
#'
#' @description Obtain the operating characteristics of The Rule-Based Design for Agents with Non-Overlapping Dose-Limiting Toxicities by simulating trials
#'
#' @usage get_oc_2DLT(p.true, p.conditional1, p.conditional2, 
#' A, B, C, a.E, a.D, b.E, b.D, c.E, ntrial = 50000, seed = NULL)
#'
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
#' @param ntrial the total number of trials to be simulated
#' @param seed the random seed for simulation
#'
#' @return \code{get_oc_2DLT()} returns the operating characteristics of the proposted design as a list,
#'         including:
#'         (1) \code{$selection_prop}: the selection percentages of each toxicity level
#'         (2) \code{$recommended_prop}: the selection percentages of each dose combination
#'         (3) \code{$nDLT_list}: a list of the number of DLT outcomes of each trial
#'         (4) \code{$nPatient_list}: a list of sample size of each trial
#'
#' @seealso Lee BL, Fan SK. A two-dimensional search algorithm for dose-finding trials of two agents. J Biopharm Stat. 2012;22(4):802-18. doi: 10.1080/10543406.2012.676587. PMID: 22651116.
#' 
#' @export

get_oc_2DLT <- function (p.true, p.conditional1, p.conditional2, A, B, C, a.E, a.D, b.E, b.D, c.E, ntrial = 50000, seed = NULL) {
  
  set.seed(seed)
  n1_dose = nrow(p.true)
  n2_dose = ncol(p.true)
  MTD_list = c()
  nDLT_list = rep(0, ntrial)
  nPatient_list = rep(0, ntrial)
  
  for (trial in 1:ntrial) {
    
    d = c(1,1)
    state_table = matrix(-1, nrow = n1_dose, ncol = n2_dose)
    n.DLTs = 0
    n.patients = 0
    
    while (dplyr::between(d[1], 1, n1_dose) &
           dplyr::between(d[2], 1, n2_dose)) {
      res = get_decision(d, p.true, p.conditional1, p.conditional2, A, B, C, a.E, a.D, b.E, b.D, c.E)
      n.DLTs = n.DLTs + res$n.DLTs
      n.patients = n.patients + res$n.patients
      
      # If the current dose combination is rejected
      if (res$flag.accepted == FALSE) {
        # rejected all dose combinations with partial orders beyond the current one
        state_table[d[1]:n1_dose, d[2]:n2_dose] = FALSE
        # if the number of patients who experienced at least one Type 1 DLT exceeds the number of patients who experienced at least one Type 2 DLT, we move from (i,j) to (i-1,j)
        if (sum(res$t.table[1, ]) > sum(res$t.table[2, ])) {
          # if the current dose level for agent 1 is not the lowest level
          if (d[1] > 1) {
            # if (i-1,j) has not been accepted yet
            if (state_table[d[1] - 1, d[2]] != 1) d = d - c(1,0) # deescalate from (i,j) to (i-1,j) 
            else break; # else, break the while loop
          }
          else break; # else, break the while loop
        }
        # if the number of patients who experienced at least one Type 2 DLT exceeds the number of patients who experienced at least one Type 1 DLT, we move from (i,j) to (i,j-1)
        else if (sum(res$t.table[1, ]) < sum(res$t.table[2, ])) {
          # if the current dose level for agent 2 is not the lowest level
          if (d[2] > 1) {
            # if (i,j-1) has not been accepted yet
            if (state_table[d[1], d[2] - 1] != 1) d = d - c(0,1) # deescalate from (i,j) to (i,j-1)
            else break; # else, break the while loop
          }
          else break; # else, break the while loop
        }
        # if tie
        else if (sum(res$t.table[1, ]) == sum(res$t.table[2, ])) {
          if (d[1] > 1 & d[2] > 1) {
            if (runif(1) < 0.5) {
              if (state_table[d[1] - 1, d[2]] != 1) d = d - c(1,0)
              else break;
            }
            else {
              if (state_table[d[1], d[2] - 1] != 1) d = d - c(0,1)
              else break;
            }
          }
          else if (d[1] > 1) {
            if (state_table[d[1] - 1, d[2]] != 1) d = d - c(1,0)
            else break;
          }
          else if (d[2] > 1) {
            if (state_table[d[1], d[2] - 1] != 1) d = d - c(0,1)
            else break;
          }
          else break;
        }
      }
      
      # If the current dose combination is accepted
      else {
        
        # # accepted all dose combinations with partial orders below the current one
        # state_table[1:d[1], 1:d[2]] = TRUE
        # if (d[1] < n1_dose & d[2] < n2_dose) {
        #       # if (i,j+1) has not been rejected yet
        #       if (state_table[d[1] + 1, d[2] + 1] != 0) d = d + c(1,1) # escalate from (i,j) to (i,j+1)
        #       else break; # else, break the while loop
        #     }
        
        # accepted all dose combinations with partial orders below the current one
        state_table[1:d[1], 1:d[2]] = TRUE
        # if the number of patients who experienced at least one Type 1 DLT exceeds the number of patients who experienced at least one Type 2 DLT, we move from (i,j) to (i-1,j)
        if (sum(res$t.table[1, ]) > sum(res$t.table[2, ])) {
          # if the current dose level for agent 2 is not the highest level
          if (d[2] < n2_dose) {
            # if (i,j+1) has not been rejected yet
            if (state_table[d[1], d[2] + 1] != 0) d = d + c(0,1) # escalate from (i,j) to (i,j+1)
            else break; # else, break the while loop
          }
          else break; # else, break the while loop
        }
        # if the number of patients who experienced at least one Type 2 DLT exceeds the number of patients who experienced at least one Type 1 DLT, we move from (i,j) to (i,j-1)
        else if (sum(res$t.table[1, ]) < sum(res$t.table[2, ])) {
          # if the current dose level for agent 1 is not the highest level
          if (d[1] < n1_dose) {
            # if (i+1,j) has not been rejected yet
            if (state_table[d[1]+1, d[2]] != 0) d = d + c(1,0) # escalate from (i,j) to (i+1,j)
            else break; # else, break the while loop
          }
          else break; # else, break the while loop
        }
        # if tie
        else if (sum(res$t.table[1, ]) == sum(res$t.table[2, ])) {
          if (d[1] < n1_dose & d[2] < n2_dose) {
            if (runif(1) < 0.5) {
              if (state_table[d[1] + 1, d[2]] != 0) d = d + c(1,0)
              else break;
            }
            else {
              if (state_table[d[1], d[2] + 1] != 0) d = d + c(0,1)
              else break;
            }
          }
          else if (d[1] < n1_dose) {
            if (state_table[d[1] + 1, d[2]] != 0) d = d + c(1,0)
            else break;
          }
          else if (d[2] < n2_dose) {
            if (state_table[d[1], d[2] + 1] != 0) d = d + c(0,1)
            else break;
          }
          else break;
        }

        
      }
      
    }
    
    MTD_candidates = which(state_table == 1, TRUE)
    if (nrow(MTD_candidates) == 0) {
      MTD_list = rbind(MTD_list, matrix(c(NA,NA), nrow = 1, dimnames = list(trial)))
      nDLT_list[trial] = n.DLTs
      nPatient_list[trial] = n.patients
    }
    else {
      MTD = which(rowSums(MTD_candidates) == max(rowSums(MTD_candidates)), TRUE)
      MTD_list = rbind(MTD_list, matrix(MTD_candidates[MTD, ],
                                        nrow = length(MTD),
                                        dimnames = list(rep(trial, length(MTD)))))
      nDLT_list[trial] = n.DLTs
      nPatient_list[trial] = n.patients
    }
    
  }
  
  colnames(MTD_list) = c("drug1", "drug2")
  MTD_list = MTD_list[!rowSums(is.na(MTD_list)),]
  recommended_prop = table(MTD_list[,1], MTD_list[,2])
  recommended_prop = round(recommended_prop/length(unique(rownames(MTD_list[which(is.na(MTD_list[,1]) == FALSE),]))),2)
  
  DLT_list = c()
  for (i in 1:nrow(MTD_list)) {
    if (is.na(MTD_list[i,1])) {
      DLT_list = append(DLT_list, NA)
    }
    else {
      DLT_list = append(DLT_list, p.true[as.numeric(MTD_list[i,1]), as.numeric(MTD_list[i,2])])
    }
  }
  selection_prop = table(DLT_list)/(ntrial - sum(is.na(DLT_list)))
  
  return (list(selection_prop = selection_prop,
               recommended_prop = recommended_prop,
               nDLT_list = nDLT_list,
               nPatient_list = nPatient_list))
}
