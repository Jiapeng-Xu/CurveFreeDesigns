#'
#' Find intervals of biological efficacious doses (BEDs)
#'
#' @description  Find the intervals of biological efficacious doses (BEDs) based on prior information, current MTD, and current data. 
#'
#' @usage findbeds_CFBD(mtd, n.eff, n.assign, a.pEff, b.pEff, E.min, gain.A, gain.AC, phi, lo)
#'
#' @param mtd the current maximum tolerated dose (MTD)
#' @param n.eff a list of no. of efficacies at all dose levels
#' @param n.assign a list of no. of patients at all dose levels
#' @param a.pEff parameter alpha for current prior distribution
#' @param b.pEff parameter beta for current prior distribution
#' @param E.min the minimum acceptable efficacy rate
#' @param gain.A the weight of reward that a BED is selected correctly; 1 is suggested
#' @param gain.AC the weight of reward that a non-BED is not selected correctly; 1 is suggested
#' @param phi the weight for how much more efficacious P_{E,i} is than E.min; 1 is suggested
#' @param lo the weight for how much less efficacious P_{E,i} is than E.min; 1 is suggested
#' 
#' @details \code{findbeds()} implement the following procedure:
#'          (1) find the interval with maximum utility A among all the possible intervals
#'          (2) find the subinterval B that maximizes the posterior probability
#'          In this function, only doses at or below MTD are investigated and we use real data rather than working data to address the BED interval.
#'
#' @return \code{findbeds()} returns (1) iL.B: the lower boundary of the interval of BEDs and; (2) iU.B: the upper boundary of the interval of BEDs and; (3) p.best: the posterior probability that interval B is acceptable
#'
#' @seealso Fan, S., Lee, B. L., & Lu, Y. (2020). A curve-free bayesian decision-theoretic design for phase Ia/Ib trials considering both safety and efficacy outcomes. \emph{Statistics in Biosciences}, 12(2), 146–166. \url{https://doi.org/10.1007/s12561-020-09272-5}
#'
#' @export

findbeds_CFBD <- function(mtd, n.eff, n.assign, a.pEff, b.pEff, E.min, gain.A, gain.AC, phi, lo) {

  # In this function, utility (called Lik) is the total goodness of fit of this step interval (model)
  max.Lik <- - .Machine$double.xmax
  
  #----------------
  # case 1: the interval is 1 piece: A=[1, mtd]
  # mtd must be at least 1
  #----------------
  Lik <-0
  for (i in 1:mtd){
    if (n.assign[i] > 0) # otherwise no data on dose i
      Lik <- Lik + n.assign[i]*((gain.A-phi*E.min)*(1-zipfR::Rbeta(E.min, a.pEff[i], b.pEff[i]))
                                +phi*a.pEff[i]/(a.pEff[i]+b.pEff[i])*(1-zipfR::Rbeta(E.min, a.pEff[i]+1, b.pEff[i])))
  }
  
  if (Lik > max.Lik) {
    max.Lik <- Lik
    iL <- 1
    iU <- mtd
  }
  
  if (mtd == 1) { # can only be one piece
    p.best <- 1-zipfR::Rbeta(E.min,a.pEff[1], b.pEff[1]) # Pr[[1,1] is acceptable]
    return(c(iL, iU, p.best)) 
  }
  
  #-----------------------------
  # case 2a: 2 pieces (high-low): A=[a=1,b], b< mtd
  # mtd must be at least 2
  #-----------------------------
  a <- 1
  for (b in 1:(mtd - 1)) {
    Lik <- 0 # initial utility for different b 
    # utility for both pieces
    # in set A = [a=1,b]
    for (i in a:b){
      if (n.assign[i] > 0) # otherwise no data on dose i
        Lik <- Lik + n.assign[i]*((gain.A-phi*E.min)*(1-zipfR::Rbeta(E.min, a.pEff[i], b.pEff[i]))
                                  +phi*a.pEff[i]/(a.pEff[i]+b.pEff[i])*(1-zipfR::Rbeta(E.min, a.pEff[i]+1, b.pEff[i])))
    }
    # in set AC = [b+1,mtd]
    for (i in (b+1):mtd){
      if (n.assign[i] > 0) # otherwise no data on dose i
        Lik <- Lik + n.assign[i]*((gain.AC+lo*E.min)*(zipfR::Rbeta(E.min, a.pEff[i], b.pEff[i]))
                                  -lo*a.pEff[i]/(a.pEff[i]+b.pEff[i])*(zipfR::Rbeta(E.min, a.pEff[i]+1, b.pEff[i])))
    }
    
    if (Lik > max.Lik) {
      max.Lik <- Lik
      iL <- a
      iU <- b
    }
  }
  
  #-----------------------------
  # case 2b: 2 pieces (low-high) A = [a, b=mtd], a > 1
  # mtd must be at least 2
  #-----------------------------
  b <- mtd
  for (a in 2:mtd) {
    Lik <- 0
    # utility for both pieces
    # in set A = [a, b=mtd]
    for (i in a:b){
      if (n.assign[i] > 0) # otherwise no data on dose i
        Lik <- Lik + n.assign[i]*((gain.A-phi*E.min)*(1-zipfR::Rbeta(E.min, a.pEff[i], b.pEff[i]))
                                  +phi*a.pEff[i]/(a.pEff[i]+b.pEff[i])*(1-zipfR::Rbeta(E.min, a.pEff[i]+1, b.pEff[i])))
    }
    
    # in set AC = [1, a-1]
    for (i in 1:(a-1)){
      if (n.assign[i] > 0) # otherwise no data on dose i
        Lik <- Lik + n.assign[i]*((gain.AC+lo*E.min)*(zipfR::Rbeta(E.min, a.pEff[i], b.pEff[i]))
                                  -lo*a.pEff[i]/(a.pEff[i]+b.pEff[i])*(zipfR::Rbeta(E.min, a.pEff[i]+1, b.pEff[i])))
    }
    
    if (Lik > max.Lik) {
      max.Lik <- Lik
      iL <- a
      iU <- b
    }
  }
  
  if (mtd == 2) { # cannot be 3 pieces
    # find best acceptable subset B now
    iL.B <- iL
    iU.B <- iU
    p.best <- 1-zipfR::Rbeta(E.min,sum(a.pEff[iL:iU]), sum(b.pEff[iL:iU])) # Pr[[iL,iU] is acceptable]
    
    for (c in iL:iU){
      for (d in c: iU){
        p.B <-1-zipfR::Rbeta(E.min,sum(a.pEff[c:d]), sum(b.pEff[c:d]))
        if (p.B > p.best){
          iL.B <- c
          iU.B <- d
          p.best <- p.B
        }
      }
    }
    return(c(iL.B,iU.B, p.best))
  }
  
  #-----------------
  # case 3: 3 pieces [1, a-1], A=[a,b], [b+1,mtd]
  # mtd must be at least 3
  #-----------------
  for (a in 2:(mtd - 1)) {
    for (b in a:(mtd - 1)) {
      Lik <- 0   
      # utility for three pieces
      # in set A = [a, b=mtd]
      for (i in a:b){
        if (n.assign[i] > 0) # otherwise no data on dose i
          Lik <- Lik + n.assign[i]*((gain.A-phi*E.min)*(1-zipfR::Rbeta(E.min, a.pEff[i], b.pEff[i]))
                                    +phi*a.pEff[i]/(a.pEff[i]+b.pEff[i])*(1-zipfR::Rbeta(E.min, a.pEff[i]+1, b.pEff[i])))
      }
      
      # in set AC = [1, a-1] & [b+1, mtd]
      for (i in 1:(a-1)){
        if (n.assign[i] > 0) # otherwise no data on dose i
          Lik <- Lik + n.assign[i]*((gain.AC+lo*E.min)*(zipfR::Rbeta(E.min, a.pEff[i], b.pEff[i]))
                                    -lo*a.pEff[i]/(a.pEff[i]+b.pEff[i])*(zipfR::Rbeta(E.min, a.pEff[i]+1, b.pEff[i])))
      }
      
      for (i in (b+1):mtd){
        if (n.assign[i] > 0) # otherwise no data on dose i
          Lik <- Lik + n.assign[i]*((gain.AC+lo*E.min)*(zipfR::Rbeta(E.min, a.pEff[i], b.pEff[i]))
                                    -lo*a.pEff[i]/(a.pEff[i]+b.pEff[i])*(zipfR::Rbeta(E.min, a.pEff[i]+1, b.pEff[i])))
      }
      
      if (Lik > max.Lik) {
        max.logLik <- logLik
        iL <- a
        iU <- b
      }
    }
  }
  iL.B <- iL
  iU.B <- iU
  p.best <- 1-zipfR::Rbeta(E.min,sum(a.pEff[iL:iU]), sum(b.pEff[iL:iU])) # Pr[[iL,iU] is acceptable]
  
  for (c in iL:iU){
    for (d in c: iU){
      p.B <-1-zipfR::Rbeta(E.min,sum(a.pEff[c:d]), sum(b.pEff[c:d]))
      if (p.B > p.best){
        iL.B <- c
        iU.B <- d
        p.best <- p.B
      }
    }
  }
  c(iL.B, iU.B, p.best)
}

