#'
#' Generate operating characteristics of the Rule-Based Design for Agents with Same Dose-Limiting Toxicities
#'
#' @description Obtain the operating characteristics of the Rule-Based Design for Agents with Same Dose-Limiting Toxicities
#'
#' @usage get_oc_2stage(p, pmin, pmax, ns, seed, stages = 2, p.preDLT = NULL)
#'
#' @details This package is developed to generate the operating characteristics of the 2 stage 2+1+3 dose combination finding design
#'          under the prespecified simulation scenarios. The two stage 2+1+3 dose combination finding design allows the first stage
#'          (m + n) to search for the potential MTD and a total of six (=m + n + h) patients to confirm and identify the MTD.
#'          
#' @param p the target toxicity level
#' @param pmin the minimum toxicity level
#' @param pmax the maximum toxicity level
#' @param ns the total number of trials to be simulated
#' @param seed the random seed for simulation
#' @param stages a number indicate whether to include the stage 0, number 2 indicates not include, otherwise, the stage 0 (Finding the starting combination) will be included
#' @param p.preDLT a list of true pre-DLT probabilities at all dose levels, if param \code{$stages} not equal to 2, \code{$p.preDLT} must be specified. Otherwise, an error will be return.
#'
#' @return \code{get_oc_2stage()} returns:
#'         (1) \code{$recommended_prop}: the selection percentages of each dose combination
#'         (2) \code{$selection_prop}: the selection percentages of each toxicity level
#'         (3) \code{$mean_DLT}: the percentage of in-trial toxicities
#'         (4) \code{$mean_patient}: patient allocation for all doses under CFBD
#'         (5) \code{$prop_no_move}: the percentage of trials that stopped at starting dose combination
#'         (6) \code{$found}: the percentage of trials recommending BEDsDLT_list
#'         (6) \code{$DLT_list}: a list of recommended toxicity level of each trial
#'
#' @seealso Fan SK, Venook AP, Lu Y. Design issues in dose-finding Phase I trials for combinations of two agents. J Biopharm Stat. 2009;19(3):509-23. doi: 10.1080/10543400902802433. PMID: 19384692.
#' 
#' @export

get_oc_2stage = function(p, pmin, pmax, ns, seed, stages = 2, p.preDLT = NULL){
  
  set.seed(seed)
  # t[i, j]: the # of toxicities at row=i, column=j 
  # n[i, j]: the # of patients at row=i, column=j 
  # r[i, j]: 0=not tested, 1=reject, 2=accept for (i,j) 
  # p.true, p_hat: the true and estimated toxicity rates   
  # mtdp: the toxicity rate of MTD 
  n1.dose = nrow(p)
  n2.dose = ncol(p) 
  none=0
  found=0 
  mtd = matrix(0, nrow = n1.dose, ncol = n2.dose)
  DLT_rates = c()
  
  # ns: # of simulations 
  # diag: the highest diag order of test
  # found[i]: indicator of finding MTD at simulation i
  # ttotal[i]: total # of toxicities at simulation i 
  # ntotal[i]: total # of patients at simulation i  
  # mtd[i, j]: # of simulations recommending (i,j) as MTD 
  
  ttotal = rep(0, ns)
  ntotal = rep(0, ns)
  
    
    for (is in 1:ns) {
      # initiate values
      
      t = matrix(0, nrow = n1.dose, ncol = n2.dose)
      n = matrix(0, nrow = n1.dose, ncol = n2.dose)
      p_hat = matrix(0, nrow = n1.dose, ncol = n2.dose)
      r = matrix(0, nrow = n1.dose, ncol = n2.dose)
      diag=0
      mtdp=-0.5 #- .Machine$double.xmax
      # I set p_hat for no data cells as -1 
      nmtd=0 
      qbar=0
      d.starting=1
      
      if (stages == 3) {
        # Stage 0: Finding the starting combination
        for (i in 1:min(n1.dose, n2.dose)) {
          if (runif(1) < p.preDLT[i] |
              runif(1) < p[i,i]) {
            n[i,i] = n[i,i] + 1
            t[i,i] = t[i,i] + 1
            d.starting = i
            break;
          }
          else {
            n[i,i] = n[i,i] + 1
            if (i == min(n1.dose, n2.dose)) {
              none = none + 1
            }
          }
          
        }
      }
            
      # stage 1: search for MTD candidates 
      # 2+1 cohort for testing (0,0): 1=reject, 2=accept 
      res = test(n1.dose, n2.dose, n, p, t, r, d.starting, d.starting);
      t = res$t
      n = res$n
      r = res$r
      if (r[1, 1]==1) none=none+1 # (0,0) is rejected, cannot move at all 
      else {
        # (0,0) is accepted; move forward 
        res = two.dim.search(p_hat,n1.dose, n2.dose, n, p, t, r, 1, 1)
        t = res$t
        n = res$n
        r = res$r
                    
        # find diag, the highest sum of row and column of tested cells 
        for (i in 1:n1.dose) {
          for (j in 1:n2.dose){
            if (n[i, j]>0 & (i+j)>diag) {
              diag=(i+j)
            } 
          }
        }
                        
        p_hat = twodim.iso(n1.dose, n2.dose, p_hat, t, n)
          #twodim.iso(n1.dose, n2.dose, p_hat, t, n)
          #twodim.iso(p_hat, t, n, diag) # find estimated toxicity rates for tested cells and lower cells 
                      
        # cells with no data are not considered 
        for (i in 1:n1.dose) {
          for (j in 1:n2.dose){
            if (n[i, j]==0) {
              p_hat[i, j]=-1
            } 
          }
        }
                    
        # stage 2: confirm those MTD candiates 
        for (i in 1:n1.dose){
          for (j in 1:n2.dose){
            if ( pmin<=p_hat[i, j] & pmax>=p_hat[i, j]){	
              for (k in 1:(6-n[i, j])){
                obs = runif(1)
                if (obs <p[i, j]) t[i, j] = t[i, j]+1
              }
              n[i, j]=6
            }
          }
        }
        p_hat = twodim.iso(n1.dose, n2.dose, p_hat, t, n)
          #twodim.iso(n1.dose, n2.dose, p_hat, t, n)
          #twodim.iso(p_hat, t, n, diag) # update the estimates 
                      
        # cells with no data are not considered 
        for (i in 1:n1.dose) {
          for (j in 1:n2.dose){
            if (n[i, j]==0) {
              p_hat[i, j]=-1
            } 
          }
        }
                    
        # recommending MTD: the acceptable trtmnts with the largest toxicity 
        for (i in 1:n1.dose) {
          for (j in 1:n2.dose) {
            if (p_hat[i, j]<=pmax & p_hat[i, j]>mtdp) mtdp=p_hat[i, j]
          }
        }
                        
                    
        if (mtdp>=0) found=found+1; # MTD found 
                      
        for (i in 1:n1.dose) {
          for (j in 1:n2.dose) {
            if (p_hat[i, j]==mtdp & low(p_hat, n1.dose, n2.dose, mtdp, i, j)!= 1){ # (i,j) is recommended as MTD 
              mtd[i, j] = mtd[i, j] + 1
              qbar = qbar + p[i, j]
              nmtd = nmtd + 1;
              DLT_rates = c(DLT_rates, p[i, j])
            }
          }
        }
                        
        qbar=qbar/nmtd
        
        # calculate the total # of toxicities and patients 
        for (i in 1:n1.dose) {
          for (j in 1:n2.dose){
            ttotal[is]=ttotal[is]+t[i, j];
            ntotal[is]=ntotal[is]+n[i, j];
          }
        }
      }
      } # end of simulations 
      
      # report mean # of toxicities and mean # of patients 
      # % of recommendation for each trtmnt and % of no MTD (none)
      
      mtdmean = mtd/ns
      tmean=sum(ttotal)/ns
      nmean=sum(ntotal)/ns
    
      none=none/ns;
      found=found/ns;
      recommended.DLT.rates = table(DLT_rates)/(found*ns)
      
      return(list(recommended_prop = mtdmean,
                  selection_prop = recommended.DLT.rates,
                  mean_DLT = tmean,
                  mean_patient = nmean,
                  prop_no_move = none,
                  found = found,
                  DLT_list = DLT_rates))
}