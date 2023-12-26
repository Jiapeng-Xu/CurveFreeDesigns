#'
#' Calculate the 2 dimensional working data for small sample to find MTD efficiently.
#'
#' @description Calculate the 2 dimensional working data for small sample to find MTD efficiently.
#'
#' @usage workdata(n.tox, n.assign)
#'
#' @param n.tox a list of no. of toxicities at all dose levels
#' @param n.assign a list of no. of patients at all dose levels
#'
#' @return \code{workdata()} returns (1) tox: the working number of toxicities and; (2) n: the working number of patients
#'
#' @seealso Fan SK, Venook AP, Lu Y. Design issues in dose-finding Phase I trials for combinations of two agents. J Biopharm Stat. 2009;19(3):509-23. doi: 10.1080/10543400902802433. PMID: 19384692.
#'
#' @examples
#' work <- workdata(c(1,1),c(2,3))
#' @export

work_data_2dim <- function(n.tox, n.assign){
  k1 <- nrow(n.tox)
  k2 <- ncol(n.tox)
  work.tox <- matrix(rep(0,times=k1*k2), nrow = k1, ncol = k2)
  work.n <- matrix(rep(0,times=k1*k2), nrow = k1, ncol = k2)
  for(i in 1:k1){
    for (j in 1:k2) {
      work.tox[i,j] <- sum(n.tox[1:i, 1:j])
      work.n[i,j] <- work.tox[i,j]+ sum(n.assign[i:k1, j:k2])-sum(n.tox[i:k1, j:k2])
    }
  }
  return(list(tox = work.tox, n = work.n))
}
