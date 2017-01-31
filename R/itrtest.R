#' Value function used for initial treatment heterogeneity assessment.  This is used 
#' inside the tree growing function. 
#' 
#' @param dat dataset being assessed 
#' @param z new (alternative) treatment assignment in the splitting procedure.  
#' @param n0 minimum number of observations needed to make a split. 
#' @return itr value from a defined split.
#' @export



itrtest <- function(dat,z,n0){
  y <- dat$y
  trt <- dat$trt
  prtx <- dat$prtx
  itr <- NA
  n <- nrow(dat)
  if (length(y)!=length(z)) stop("the vector z must have the same length as data.")
  if(n > n0) {
    itr <- mean(trt*y*z/prtx+(1-trt)*y*(1-z)/(1-prtx))/mean(trt*z/prtx+(1-trt)*(1-z)/(1-prtx))
  }
  itr
}
