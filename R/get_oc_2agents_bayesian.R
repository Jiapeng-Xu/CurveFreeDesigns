#'
#' Generate operating characteristics for the Bayesian Decision-Theoretic Design
#'
#' @description Obtain the operating characteristics of the Bayesian Decision-Theoretic Design
#'
#' @usage get_oc_2agents_bayesian(pTox, target, T.max, n.min, n.max, 
#' n.sim, seed, var.ratio = 4, alpha = 1.2, eta = 1, r1 = 0.5, r2 = 0.95, type = 1)
#'
#' @param pTox a list of true toxicity probabilities at all dose levels
#' @param target the target toxicity level
#' @param T.max the absolute overly toxic level, usually .05 above target
#' @param n.min the minimal sample size
#' @param n.max the maximal sample size
#' @param n.sim the total number of trials to be simulated
#' @param seed the random seed for simulation
#' @param var.ratio an equivalent no. of patients contained in the prior information, chosen such that prior.n.mtd * target is larger than but close to 1
#' @param alpha the weight of penalty of toxicity rate below target level; 1.2 is suggested
#' @param eta the weight of penalty of toxicity rate above target level; 1 is suggested
#' @param r1 the error rate for concluding inadmissible; .10 is suggested
#' @param r2 the error rate for concluding admissible; .10 is suggested
#' @param type the type of simulation, 1. Simulation study; 2. Sensitivity analysis based on random error (underestimate); 3. Sensitivity analysis based on random error (overestimate); 4. Sensitivity analysis based on fixed error
#'
#' @return \code{get_oc()} returns the operating characteristics of the CFBD design as a list,
#'         including:
#'         (1) \code{$n}: the average sample size over all simulated trials
#'         (2) \code{$selection_prop}: the selection percentages of each toxicity level
#'         (3) \code{$percentFound}: the percentage of trials recommending BEDs
#'         (5) \code{$percentTox}: the percentage of in-trial toxicities
#'         (6) \code{$percentMTD}: the selection percentages of each dose combination
#'         (7) \code{$percentPatients}: patient allocation for all doses under CFBD
#'
#' @seealso Fan, S., Lee, B. L., & Lu, Y. (2020). A curve-free bayesian decision-theoretic design for phase Ia/Ib trials considering both safety and efficacy outcomes. \emph{Statistics in Biosciences}, 12(2), 146–166. \url{https://doi.org/10.1007/s12561-020-09272-5}
#'
#' @export


get_oc_2agents_bayesian <- function (pTox, target, T.max, n.min, n.max, n.sim, seed, var.ratio = 4, alpha = 1.2, eta = 1, r1 = 0.5, r2 = 0.95, type = 1) {
  if (target < 0.05) {
    stop("the target is too low")
  }
  if (target > 0.6) {
    stop("the target is too high")
  }
  if (min(pTox) < 0 | max(pTox) > 1) {
    stop("the true toxicity probabilities are not valid, need to be larger than 0 and less than 1")
  }
  
  set.seed(seed)
  n1.dose <- nrow(pTox)
  n2.dose <- ncol(pTox)
  
  #--------------------
  # initialize storage
  #--------------------
  sample.size <- numeric(n.sim)      # sample size for each simulation
  is.MTD      <- numeric(n.sim)      # finding a MTD or not per simulation
  is.good     <- numeric(n.sim)      # BED.eff is above E.min or not per simulation
  DLT.list    <- numeric(n.sim)
  MTD.list    <- matrix(rep(0, 2*n.sim), nrow = n.sim)
  tot.tox     <- matrix(rep(0, n1.dose*n2.dose), nrow = n1.dose, ncol = n2.dose)     # total no. of toxicities at each dose in all simulations
  tot.assign  <- matrix(rep(0, n1.dose*n2.dose), nrow = n1.dose, ncol = n2.dose)     # total no. of patients assigned to each dose in all simulations
  tot.start   <- matrix(rep(0, n1.dose*n2.dose), nrow = n1.dose, ncol = n2.dose)     # total no. of selections as starting dose in all simulation
  tot.select  <- matrix(rep(0, n1.dose*n2.dose), nrow = n1.dose, ncol = n2.dose)     # total no. of selections as mtd in all simulations

  #-------------------
  # simulation begins
  #-------------------
  for (s in 1:n.sim) {
    if (type == 1) {
      a.pTox <- var.ratio * pTox
      b.pTox <- var.ratio * (1 - pTox)
      }
    else if (type == 2) {
      underestimated_pTox = (pTox - matrix(runif(n1.dose*n2.dose, 0, 1)*0.05, nrow = n1.dose, ncol = n2.dose))
      underestimated_pTox[underestimated_pTox < 0] = 0
      a.pTox = var.ratio * underestimated_pTox
      b.pTox <- var.ratio * (1 - underestimated_pTox)
    }
    else if (type == 3) {
      overestimated_pTox = (pTox + matrix(runif(n1.dose*n2.dose, 0, 1)*0.05, nrow = n1.dose, ncol = n2.dose))
      overestimated_pTox[overestimated_pTox > 1] = 1
      a.pTox = var.ratio * overestimated_pTox
      b.pTox <- var.ratio * (1 - overestimated_pTox)
    }
    else {
      misestimated_pTox = (pTox + matrix(runif(n1.dose*n2.dose, -1, 1)*0.05, nrow = n1.dose, ncol = n2.dose))
      misestimated_pTox[misestimated_pTox < 0] = 0
      misestimated_pTox[misestimated_pTox > 1] = 1
      a.pTox = var.ratio * misestimated_pTox
      b.pTox <- var.ratio * (1 - misestimated_pTox)
    }
    
    
    n.tox <- matrix(rep(0, n1.dose*n2.dose), nrow = n1.dose, ncol = n2.dose)  # total no. of toxicities at each dose combination in each simulation
    n.assign <- matrix(rep(0, n1.dose*n2.dose), nrow = n1.dose, ncol = n2.dose) # total no. of patients assigned to each dose combination in each simulation
    
    # random escalation to find starting dose
    # escalating from dose one until toxicity presents
    
    d = c(1,1)
    while (d[1] < n1.dose &
           d[2] < n2.dose) {
      # enrolled one patient at this dose
      sample.size[s] <- sample.size[s] + 1
      tot.assign[d[1], d[2]] <- tot.assign[d[1], d[2]] + 1
      n.assign[d[1], d[2]] <- n.assign[d[1], d[2]] + 1
      
      # generate the patient's responses
      # update prior parameters of toxicity rate using WORKING data
      # if toxicity presents, this dose is the starting dose
      # otherwise, move up to the next dose
      if (runif(1) < pTox[d[1], d[2]]) {
        tot.start[d[1], d[2]] <- tot.start[d[1], d[2]] + 1
        n.tox[d[1], d[2]] <- n.tox[d[1], d[2]] + 1
        a.pTox[d[1]:nrow(a.pTox), d[2]:ncol(a.pTox)] = a.pTox[d[1]:nrow(a.pTox), d[2]:ncol(a.pTox)] +1
        break;
      }
      else {
        b.pTox[1:d[1], 1:d[2]] <- b.pTox[1:d[1], 1:d[2]] + 1
        
        # random escalation to find starting dose
        if(runif(1) <= 0.5) {
          d = d + c(0,1)
        }
        else {
          d = d + c(1,0)
        }
      }
    }
    
    repeat {
      # enroll/apply FLW on at least n.min.mtd patients first
      # stop trial for overly toxic or running out of patients for stage 1
      if (sample.size[s] >= n.max) break
      
      is.stop.S3 <- (1 - pbeta(T.max, a.pTox[1, 1], b.pTox[1, 1])) > r1
      
      p.min = c()
      for (i in d[1] : n1.dose) {
        for (j in d[2] : n2.dose) {
          p.min = c(p.min, 1 - pbeta(T.max, a.pTox[i, j], b.pTox[i, j]))
        }
      }
      p.min = p.min[-1]
      is.stop.S4 = min(p.min) > r2
      
      is.stop = is.stop.S3 | is.stop.S4
      
      if (sample.size[s] >= n.min && is.stop) break
      
      # continue trial
      # assign one patient at the current mtd; update parameters at the dose
      d = dose_escalation(d, a.pTox, b.pTox, target, alpha, eta)
      
      sample.size[s] <- sample.size[s] + 1
      tot.assign[d[1], d[2]] <- tot.assign[d[1], d[2]] + 1
      n.assign[d[1], d[2]] <- n.assign[d[1], d[2]] + 1
      
      # generate responses
      is.tox <- runif(1) < pTox[d[1], d[2]]
      
      # update the real data: no. of efficacy, no. of toxicity at this dose
      # update prior parameters of toxicity rate using WORKING data
      
      if (is.tox) {
        n.tox[d[1], d[2]] <- n.tox[d[1], d[2]] + 1
        a.pTox[d[1]:nrow(a.pTox), d[2]:ncol(a.pTox)] = a.pTox[d[1]:nrow(a.pTox), d[2]:ncol(a.pTox)] +1
      } 
      else {
        b.pTox[1:d[1], 1:d[2]] <- b.pTox[1:d[1], 1:d[2]] + 1
      }
    }
    
    # update simulation parameters: no. of patients, toxicities, efficacies, and mtd selections
    tot.tox <- tot.tox + n.tox
    
    if (is.stop.S3 == FALSE) { # no mtd, all overly toxic
      is.MTD[s] = 1
      DLT.list[s] = pTox[d[1], d[2]]
      MTD.list[s, ] = d
      tot.select[d[1], d[2]] <- tot.select[d[1], d[2]] + 1
      if (pTox[d[1], d[2]] >= target & pTox[d[1], d[2]] <= T.max) {is.good[s] = 1}
    }
    
} # end of simulation
  
  #-------------------
  # summary statistics
  #-------------------
  n = mean(sample.size) # mean sample size
  percentFound = round(mean(is.MTD) * 100, digits = 1) # percent of conclusion of MTD found
  percentSuccess = sum(tot.select[which(pTox < round(target + 0.1, 1) & pTox > round(target - 0.1, 1))])/n.sim*100   # percent that pooled efficacy above E.min
  percentTox = round(100 * sum(tot.tox)/ sum(tot.assign), digits = 1)
  percentMTD = round(100 * tot.select / n.sim, 1) # distribution of mtd recommendation
  percentPatients = round(100 * tot.assign/ sum(tot.assign), digits = 1) # distributions of patients' assignment and their responses
  selection_prop = merge(as.data.frame(table(DLT.list[DLT.list != 0])/sum(is.MTD)), data.frame('Var1' = factor(unique(as.vector(pTox)))), by = 'Var1', all.y = TRUE)
  selection_prop[is.na(selection_prop)] = 0
  

  return (list(n = n,
               selection_prop = selection_prop,
               percentFound = percentFound,
               percentTox= percentTox,
               percentMTD = percentMTD,
               percentPatients = percentPatients))
}
