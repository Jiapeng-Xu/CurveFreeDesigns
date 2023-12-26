#'
#' Update isotonic estimation of each dose combination
#'
#' @description Update isotonic estimation of each dose combination
#'
#' @usage twodim.iso(n1.dose, n2.dose, p_hat, t, n)
#'
#' @param n1.dose the number of dose level of agent 1
#' @param n2.dose the number of dose level of agent 2
#' @param p_hat a matrix containing the isotonic estimation of DLT probability of each dose level
#' @param t a matrix containing the number of DLT outcomes of each dose level
#' @param n a matrix containing the total number of patients of each dose level
#'
#' @return \code{two.dim.search()} returns: p_hat a matrix containing the isotonic estimation of DLT probability of each dose level
#'
#' @seealso Fan SK, Venook AP, Lu Y. Design issues in dose-finding Phase I trials for combinations of two agents. J Biopharm Stat. 2009;19(3):509-23. doi: 10.1080/10543400902802433. PMID: 19384692.
#' 
#' @export

twodim.iso = function(n1.dose, n2.dose, p_hat, t, n) {
  for (i in 1:(n1.dose+n2.dose-1)) {
    for (j in 1:i){
      if (j >=1 & 
          j <= n1.dose & 
          (i+1)-j >= 1 & 
          (i+1)-j <= n2.dose) {
        row = j
        col = (i+1)-j
        
        # submatrix 
        # t_sub = t[1:row, 1:col]
        # n_sub = n[1:row, 1:col]
        #n.sum = 0
        #t.sum = 0
        
        for (k in (row+col-1):1) {
          for (l in 1:k) {
            if (l >= 1 &
                l <= row &
                (k+1)-l >=1 &
                (k+1)-l <= col)
            {
              #n.sum = n.sum + n[l,(k+1)-l]
              #t.sum = t.sum + t[l,(k+1)-l]
              n.sum = sum(n[l:row, ((k+1)-l):col])
              t.sum = sum(t[l:row, ((k+1)-l):col])
              p_hat[row, col] = max(p_hat[row, col], ifelse(is.na(t.sum/n.sum), 0, t.sum/n.sum))
            }
          }
         
        }
        
        if (i>1) {
          for (a in (row+col-2):1) {
            for (b in 1:a) {
              if (b >= 1 &
                  b <= row &
                  (a+1)-b >=1 &
                  (a+1)-b <= col)
              {
                if((a+1)-b + 1<= col & b + 1 <= row) {
                  p_hat[b, (a+1)-b] = min(p_hat[b, (a+1)-b], p_hat[b+1, (a+1)-b], p_hat[b, (a+1)-b+1])
                }
                else if ((a+1)-b + 1 <= col) {
                  p_hat[b, (a+1)-b] = min(p_hat[b, (a+1)-b], p_hat[b, (a+1)-b+1])
                }
                else if (b + 1 <= row) {
                  p_hat[b, (a+1)-b] = min(p_hat[b, (a+1)-b], p_hat[b+1, (a+1)-b])
                }
              }
            }
          }
        }
      }
    }
  }
  return(p_hat)
}
