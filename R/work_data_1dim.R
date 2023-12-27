#'
#' Calculate the working data for small sample to find MTD efficiently.
#'
#' @description Calculate the working data for small sample to find MTD efficiently.
#'
#' @usage work_data_1dim(n.tox, n.assign)
#'
#' @param n.tox a list of no. of toxicities at all dose levels
#' @param n.assign a list of no. of patients at all dose levels
#'
#' @return \code{work_data_1dim()} returns (1) tox: the working number of toxicities and; (2) n: the working number of patients
#'
#' @seealso Fan, S., Lee, B. L., & Lu, Y. (2020). A curve-free bayesian decision-theoretic design for Phase IA/IB trials considering both safety and efficacy outcomes. \emph{Statistics in Biosciences}, 12(2), 146–166. \url{https://doi.org/10.1007/s12561-020-09272-5}
#'          Fan SK, Lu Y, Wang YG. A simple Bayesian decision-theoretic design for dose-finding trials. \emph{Stat Med}, 2012 Dec 10;31(28):3719-30. doi: 10.1002/sim.5438. Epub 2012 Jul 5. PMID: 22763943. \url{https://doi.org/10.1002/sim.5438}
#'
#' @examples
#' work <- work_data_1dim(c(1,1),c(2,3))
#' @export

work_data_1dim <- function(n.tox, n.assign){
  k <- length(n.tox)
  work.tox <- rep(0,times=k)
  work.n <- rep(0, times=k)
  for(i in 1:k){
    work.tox[i] <- sum(n.tox[1:i])
    work.n[i] <- work.tox[i]+ sum(n.assign[i:k])-sum(n.tox[i:k])
  }
  return(list(tox = work.tox,n = work.n))
}
