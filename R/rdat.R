#' @title Simulates data from an RCT design
#' @description Simulates data from an RCT according to the following model:
#' 2 + 2*sign(x1<cut2) + beta1*trt*subgrp + beta2*(1-trt)*subgrp + N(0,sigma)
#' If depth=1, then subgrp=(x1<=cut1)
#' If depth!=1 then subgrp=(x1>=0.3 & x3>=0.1)
#' 
#' @param n size of the dataset to be generate.  Defaults to 100.
#' @param depth gives the number of interacting covariates. If set to 1, then 
#'  then covariate X1 interacts with treatment. If set to another value, then 
#'  covariates X1 and X3 both interact with treatment effect (one-way interactions). Defaults to 1. 
#' @param beta1 controls the strength of the treatment effect. Defaults to 2. 
#' @param beta2 controls the strength of the noise. Defaults to 2. 
#' @param K internal variable used to generate the fineness of the unit interval on which covariates are generated. Defaults to 50.
#' @param sigma controls standard deviation of random variation.  Defaults to 1. 
#' @param cut1,cut2,cut3,cut4 controls where the cutpoints are to define subgroups. 
#' @return data frame containing y (outcome), x1-x4 (covariates), trt (treatment), prtx (probability of being in treatment group), and ID variable
#' @export
#' @examples
#' 
#' # This generates a dataframe with 500 observations, X1 as the only variable interacting with 
#' # the treatment, and a signal to noise ratio of 2/2=1.
#' data<-rdat(n=500)

rdat<-function(n = 100, 
               K = 50, 
               beta1 = 2, 
               beta2 = 2,
               sigma = 1, 
               cut1 = 0.5, 
               cut2 = 0.5, 
               cut3 = 0.3, 
               cut4 = 0.1, 
               depth = 1)
{
  trt <- sample(c(0,1), n, replace=T)
  #### Generate Covariates
  x1 <- x2 <- x3 <- x4 <- NULL
  for (j in 1:4) {
    assign(paste("x", j, sep=""),sample(1:K, n, replace=T)/K)
  }
  ### 
  if(depth==1){
    mean <- 2 + 2*sign(x2<cut2) + beta1*sign(x1<cut1)*trt + beta2*sign(x1>=cut1)*(1-trt)
    ##### Output
  }else{ 
    mean <- 2 + 2*sign(x2<cut2) + beta1*sign(x1>cut3 & x3>cut4)*trt + beta2*(1-sign(x1>cut3 & x3>cut4))*(1-trt)
  }
  y <- mean + rnorm(n, mean=0, sd=sigma)
  data.frame(x1=x1, x2=x2, x3=x3, x4=x4, y=y, id=1:n, trt=trt,prtx=rep(0.5,n))
}