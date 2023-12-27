#'
#' Conduct sensitivity analysis of Curve-Free Design for Maximum Tolerated Dose in One-Agent Trials based on random errors
#'
#' @description Conduct sensitivity analysis by adding random errors to the means of the Beta prior distributions to investigate the robustness of proposed design. The random errors follow uniform distribution from \code{lower} to \code{upper}.
#'
#' @usage sensitivity_randomError_FLW(pTox, error.T, var.ratio, target, T.max, 
#' n.min.mtd, n.max.mtd, n.sim, seed, alpha = 1, eta = 1, p1 = 0.1, p2 = 0.1)
#'
#' @param pTox a list of true toxicity probabilities at all dose levels
#' @param error.T error size for toxicity rates; taken from -0.5 - 0.5
#' @param var.ratio an equivalent no. of patients contained in the prior information, chosen such that prior.n.mtd * target is larger than but close to 1
#' @param target the target toxicity rate
#' @param T.max the absolute overly toxic, usually .05 above target
#' @param n.min.mtd the minimal sample size for MTD identification (phase Ia)
#' @param n.max.mtd the maximal sample size for MTD identification (phase Ia)
#' @param n.sim the total number of trials to be simulated
#' @param seed the random seed for simulation
#' @param alpha the weight of penalty of toxicity rate below target level; 1 is suggested
#' @param eta the weight of penalty of toxicity rate above target level; 1 is suggested
#' @param p1 the error rate for concluding inadmissible; .10 is suggested
#' @param p2 the error rate for concluding admissible; .10 is suggested
#'
#' @return \code{sensitivity_randomError_FLW()} returns the operating characteristics of proposed design as a list,
#'         including:
#'         (1) \code{$n}: the average sample size over all simulated trials
#'         (2) \code{$percentFound}: the percentage of trials recommending BEDs
#'         (3) \code{$percentCorrect}: within those trials recommending BEDs, the percentage of trials of which all recommended doses are truly admissible and acceptable
#'         (4) \code{$percentTox}: the percentage of in-trial toxicities
#'         (5) \code{$percentMTD}: the selection percentages of MTD
#'         (6) \code{$percentPatients}: patient allocation for all doses under CFBD
#'
#' @seealso Fan SK, Lu Y, Wang YG. A simple Bayesian decision-theoretic design for dose-finding trials. \emph{Stat Med}, 2012 Dec 10;31(28):3719-30. doi: 10.1002/sim.5438. Epub 2012 Jul 5. PMID: 22763943. \url{https://doi.org/10.1002/sim.5438}
#'
#' @export


sensitivity_randomError_FLW <- function (pTox, error.T, var.ratio, target, T.max, n.min.mtd, n.max.mtd, n.sim, seed, alpha = 1, eta = 1, p1 = 0.1, p2 = 0.1) {
  
  if (target < 0.05) {
    stop("the target is too low!")
  }
  if (target > 0.6) {
    stop("the target is too high!")
  }
  if (min(pTox) < 0 | max(pTox) > 1) {
    stop("the true toxicity probabilities are not valid, need to be larger than 0 and less than 1")
  }
  
  set.seed(seed)
  n.dose <- length(pTox)
  
  #--------------------
  # initialize storage
  #--------------------
  sample.size <- numeric(n.sim)      # sample size for each simulation
  is.MTD      <- numeric(n.sim)      # finding a BED interval or not per simulation
  is.good     <- numeric(n.sim)      # Toxicity rate is below T.max or not per simulation
  tot.assign  <- numeric(n.dose)     # total no. of patients assigned to each dose in all simulations
  tot.start   <- numeric(n.dose)     # total no. of selections as starting dose in all simulation
  tot.select  <- numeric(n.dose + 1) # total no. of selections as mtd in all simulations
  tot.tox     <- numeric(n.dose)     # total no. of toxicities at each dose in all simulations
  n.good <- 0 # no. of good MTD recommendation in all simulations
  
  #-------------------
  # simulation begins
  #-------------------
  for (s in 1:n.sim) {
    
    mu.T <- numeric(n.dose)
    for (i in 1: n.dose) {
      mu.T[i]<-pTox[i] * (1+runif(1,-1,1) * error.T)
    }
    mu.T[mu.T > 1] = 1
    mu.T[mu.T < 0] = 0
    a.pTox <- var.ratio * mu.T
    b.pTox <- var.ratio * (1 - mu.T)
    
    n.tox <- numeric(n.dose)
    n.assign <- numeric(n.dose)
    
    # random escalation to find starting dose
    # escalating from dose one until toxicity presents
    
    for (dose in 1:n.dose) {
      # enroll one patient at this dose
      sample.size[s] <- sample.size[s] + 1
      tot.assign[dose] <- tot.assign[dose] + 1
      n.assign[dose] <- n.assign[dose] + 1
      
      # generate the patient's responses  
      is.tox <- runif(1) < pTox[dose]
      
      # update the real data: no. of efficacy, no. of toxicity at this dose
      # update prior parameters of toxicity rate using WORKING data
      
      # if toxicity presents, this dose is the starting dose
      # otherwise, move up to the next dose
      
      if (is.tox) {
        tot.start[dose] <- tot.start[dose] + 1
        n.tox[dose] <- n.tox[dose] + 1
        a.pTox[dose:n.dose] <- a.pTox[dose:n.dose] + 1
        break # starting dose is found, so quit the loop
      } else
        b.pTox[1:dose] <- b.pTox[1:dose] + 1
    }
    
    # stage 1: finding mtd with a stopping rule 
    osla <- findmtd(target, a.pTox, b.pTox, alpha, eta)
    mtd <- osla$mtd
    
    repeat {
      # enroll/apply FLW on at least n.min.mtd patients first
      # stop trial for overly toxic or running out of patients for stage 1
      if (sample.size[s] >= n.max.mtd) break
      is.stop <- stopmtd(a.pTox, b.pTox, mtd, target, T.max, p1, p2)
      if (sample.size[s] >= n.min.mtd && is.stop) break
      
      # continue trial
      # assign one patient at the current mtd; update parameters at the dose
      sample.size[s] <- sample.size[s] + 1
      tot.assign[mtd] <- tot.assign[mtd] + 1
      n.assign[mtd] <- n.assign[mtd] + 1
      
      # generate responses
      is.tox <- runif(1) < pTox[mtd]
      
      # update the real data: no. of efficacy, no. of toxicity at this dose
      # update prior parameters of toxicity rate using WORKING data
      
      if (is.tox) {
        n.tox[mtd] <- n.tox[mtd] + 1
        a.pTox[mtd:n.dose] <- a.pTox[mtd:n.dose] + 1
      } else b.pTox[1:mtd] <- b.pTox[1:mtd] + 1
      
      osla <- findmtd(target, a.pTox, b.pTox, alpha, eta)
      mtd <- osla$mtd
    }
    
    # update simulation parameters: no. of patients, toxicities, efficacies, and mtd selections
    tot.tox <- tot.tox + n.tox
    
    if (is.stop == 2) { # no mtd, all overly toxic
      tot.select[n.dose + 1] <- tot.select[n.dose + 1] + 1
    }
    else {
      is.MTD[s] <-1 # find a MTD in simulation s
      tot.select[mtd] <- tot.select[mtd] + 1
      if (pTox[mtd] <= T.max) n.good <- n.good + 1
      if (n.assign[mtd] != 0 & n.tox[mtd]/n.assign[mtd] <= T.max) is.good[s] <-1 # sample pooled eff above E.min
    }
  } # End of simulation
  
  # Simulation Results
  n.MTD <- sum(is.MTD) # no. of finding BED in all simulations
  
  #-------------------
  # summary statistics
  #-------------------
  n = mean(sample.size) # mean sample size
  percentFound = round((n.MTD / n.sim) * 100, digits = 1) # percent of conclusion of MTD found
  percentCorrect = round((n.good / n.MTD) * 100, digits = 1) # percent that all recommended MTDs are acceptable given a MTD is found 
  percentSuccess = round((sum(is.good)/n.MTD)*100, digits=1) # percent that toxicity rate below T.max
  percentTox = round(100 * sum(tot.tox)/ sum(tot.assign), digits = 1)
  percentMTD = round(100 * tot.select[1:n.dose] / n.sim, 1) # distribution of mtd recommendation
  percentPatients = round(100 * tot.assign/ sum(tot.assign), digits = 1) # distributions of patients' assignment and their responses
  
  
  return (list(n = n,
               percentFound = percentFound,
               percentCorrect = percentCorrect,
               percentTox= percentTox,
               percentMTD = percentMTD,
               percentPatients = percentPatients))
}
