#' Simulates data from an observational study according to the following model:
#' 1 + 2*X2 + 4*X4 + beta1*trt*subgrp + beta2*(1-trt)*(1-subgrp) + N(0,1)
#' where subgrp is the group of interacting variable(s).
#' If depth=1, then subgrp=(X1 < 0.5)
#' If depth!=1 then subgrp=(X1>0.3 & X3>0.1)
#' 
#' @param n size of the dataset to be generate.  Required input. 
#' @param depth gives the number of interacting covariates. If set to 1, then 
#'  then covariate X1 interacts with treatment. If set to another value, then 
#'  covariates X1 and X3 both interact with treatment effect (one-way interactions). Required input.
#' @param beta1 controls the strength of the treatment effect. Required input. 
#' @param beta2 controls the strength of the noise. Required input. 
#' @param cut1/cut2/cut3 controls the cutpoints for defining subgroups. 
#' @return dataframe containing y (outcome), X1-X4 (covariates), trt (treatment), prtx (probability of being in treatment group)
#' @export
#' @examples
#' data<-gdataM(n=500, depth=1, beta1=1, beta2=2)
#' This generates a dataframe with 500 observations, X1 as the only variable interacting with 
#' the treatment, and a signal to noise ratio of 1/2.


gdataM <- function(n,depth, beta1, beta2, 
                   cut1=0.5, cut2=0.3, cut3=0.1, K=50){
  NX  <- 4
  NPATIENT  <- n
  for (j in 1:4) {
    assign(paste("X", j, sep=""),sample(1:K, n, replace=T)/K)
  }
  covariatesX <- matrix(c(X1,X2,X3,X4), nrow = NPATIENT, ncol = NX)
  expLogit  <- exp(-4+3*covariatesX[,1]+5*covariatesX[,3])
  treatmentProbT  <- round(expLogit/(1+expLogit), 3)
  treatmentT  <- rbinom(NPATIENT,1,treatmentProbT)
  #fit <- glm(treatmentT~covariatesX, family = binomial(link = "logit"))
  #preds <- predict(fit, data.frame(covariatesX))
  #prtx <- round(exp(preds)/(1+exp(preds)),2)

  if(depth==1){
    subGroupIndex  <- ( covariatesX[,1] < cut1)
  }else {
    subGroupIndex  <- ( covariatesX[,1] > cut2 & covariatesX[,3] > cut3)
  }
  responseY1Mean  <- 1 + 2*covariatesX[,2] + 4*covariatesX[,4] + beta1*(subGroupIndex)*treatmentT
  responseY1  <- responseY1Mean  +  rnorm(NPATIENT)
  
  responseY0Mean  <- 1 + 2*covariatesX[,2] + 4*covariatesX[,4] + beta2*(1-subGroupIndex)*(1-treatmentT)
  responseY0  <- responseY0Mean + rnorm(NPATIENT)
  
  responseY  <- treatmentT*responseY1+(1-treatmentT)*responseY0

  dataM  <- as.data.frame(cbind(covariatesX,responseY,treatmentT,treatmentProbT, 1:n))
  names(dataM)  <- c(paste("X",c(1:4), sep=""),"y","trt","prtx", "id")
  dataM
}
