#' @title Value function used for treatment assignment value calculation.  
#' 
#' @description This function is used inside the tree growing functions.  
#' 
#' @param dat dataset being assessed 
#' @param z new (alternative) treatment assignment in the splitting procedure.  
#' @param n0 minimum number of observations needed to make a split. 
#' @param aug logical indicator for AIPWE estimation (TRUE) or IPWE estimation (FALSE).
#' @return itr value from a defined split.
#' @examples 
#' dat <- gdataM(n=1000, depth=2, beta1=3, beta2=1)
#' #Assign every other patient to treatment
#' z <- rep(c(0,1), nrow(dat)/2)
#' IPWE value
#' itrtest(dat, z, 5, F)
#' AIPWE value
#' itrtest(dat, z, 5, T)
#' @export



itrtest <- function(dat,z,n0,aug){
  y <- dat$y
  trt <- dat$trt
  prtx <- dat$prtx
  
  if(aug){
    if(length(trt==1)>0 & length(trt==0)>0){
      m1 <- mean(y[trt==1])
      m0 <- mean(y[trt==0])
    }
  }
  itr <- NA
  n <- nrow(dat)
  if (length(y)!=length(z)) stop("the vector z must have the same length as data.")
  #if(min(sum(trt), sum(1-trt)) >= n0) {
  if(!aug) itr <- mean(trt*y*z/prtx+(1-trt)*y*(1-z)/(1-prtx))/mean(trt*z/prtx+(1-trt)*(1-z)/(1-prtx))
  if(aug){
    first <- (trt*z+(1-trt)*(1-z))*y/(prtx*trt+(1-prtx)*(1-trt))
    second <- (((trt*z+(1-trt)*(1-z)) - (prtx*trt+(1-prtx)*(1-trt)))/((prtx*trt+(1-prtx)*(1-trt))))*(m1*z+m0*(1-z))
    itr <- mean(first-second)
    #}
  }
  return(round(itr,4))
}