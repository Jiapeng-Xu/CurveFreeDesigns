% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_oc_CFBD.R
\name{get_oc_CFBD}
\alias{get_oc_CFBD}
\title{Generate operating characteristics for Curve-Free Bayesian Design (CFBD)}
\usage{
get_oc_CFBD(pTox, pEff, var.ratio, var.ratio.E, target, T.max, E.min, 
n.min.mtd, n.max.mtd, n.min.int, n.max.int, n.cohort, n.sim, seed, alpha = 1, eta = 1,
p1 = 0.1, p2 = 0.1, gain.A = 1, gain.AC = 1, phi = 1, lo = 1, q1 = 0.1, q2 = 0.1)
}
\arguments{
\item{pTox}{a list of true toxicity probabilities at all dose levels}

\item{pEff}{a list of true efficacy probabilities at all dose levels}

\item{var.ratio}{an equivalent no. of patients contained in the prior information, chosen such that prior.n.mtd * target is larger than but close to 1}

\item{var.ratio.E}{an equivalent no. of patients contained in the prior information, chosen such that prior.n.beds * E.min is larger than but close to 1}

\item{target}{the target toxicity level}

\item{T.max}{the absolute overly toxic level, usually .05 above target}

\item{E.min}{the minimum acceptable efficacy rate}

\item{n.min.mtd}{the minimal sample size for MTD identification (phase Ia)}

\item{n.max.mtd}{the maximal sample size for MTD identification (phase Ia)}

\item{n.min.int}{the minimum sample size to stop a trial (phase Ia/Ib)}

\item{n.max.int}{the maximum sample size to stop a trial (phase Ia/Ib)}

\item{n.cohort}{the size of a confirmation cohort per BED}

\item{n.sim}{the total number of trials to be simulated}

\item{seed}{the random seed for simulation}

\item{alpha}{the weight of penalty of toxicity rate below target level; 1 is suggested}

\item{eta}{the weight of penalty of toxicity rate above target level; 1 is suggested}

\item{p1}{the error rate for concluding inadmissible; .10 is suggested}

\item{p2}{the error rate for concluding admissible; .10 is suggested}

\item{gain.A}{the weight of reward that a BED is selected correctly; 1 is suggested}

\item{gain.AC}{the weight of reward that a non-BED is not selected correctly; 1 is suggested}

\item{phi}{the weight for how much more efficacious P_{E,i} is than E.min; 1 is suggested}

\item{lo}{the weight for how much less efficacious P_{E,i} is than E.min; 1 is suggested}

\item{q1}{the threshold of the posterior probability of a dose being acceptable; 0.1 is suggested}

\item{q2}{the threshold of the posterior probability of a dose being unacceptable ; 0.1 is suggested}
}
\value{
\code{get_oc_CFBD()} returns the operating characteristics of the CFBD design as a list,
        including:
        (1) \code{$percentFound}: the percentage of trials recommending BEDs
        (2) \code{$percentCorrect}: within those trials recommending BEDs, the percentage of trials of which all recommended doses are truly admissible and acceptable
        (3) \code{$percentSuccess}: within those trials recommending BEDs, the percentage of trials of which the sample pooled efficacy rate being above \code{E.min}
        (4) \code{$percentToxicity}: the percentage of in-trial toxicities
        (5) \code{$percentEfficacy}: the percentage of in-trial efficacies
        (6) \code{$averageSampleSize}: the average sample size over all simulated trials
        (7) \code{$percentMTD}: the selection percentages of MTD
        (8) \code{$percentL}: the selection percentages of lower boundary of BEDs interval
        (9) \code{$percentU}: the selection percentages of upper boundary of BEDs interval
        (10) \code{$percentPatients}: patient allocation for all doses under CFBD
}
\description{
Obtain the operating characteristics of the CFBD design by simulation
}
\details{
The operating characteristics of the CFBD design are generated by simulating trials under the prespecified true toxicity and efficacy probabilities.
         We assume that the dose-toxicity curve is monotonic non-decreasing and the dose-efficacy curve is unimodal. 
         In phase Ia stage, a Bayesian model on the toxicity rates is used to locate the maximum tolerated doese.
         In phase Ib stage, we model the dose-efficacy curve using a step function while continuing to monitor the toxicity rates. 
         At the end of phase Ib, if some doses are recommended as BEDs, a cohort of confirmation is recruited and assigned at these doses to improve the precision of estimates at these doses.
}
\seealso{
Fan, S., Lee, B. L., & Lu, Y. (2020). A curve-free bayesian decision-theoretic design for phase Ia/Ib trials considering both safety and efficacy outcomes. \emph{Statistics in Biosciences}, 12(2), 146–166. \url{https://doi.org/10.1007/s12561-020-09272-5}
}
