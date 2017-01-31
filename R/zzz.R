.onAttach <- function(libname, pkgname){
  options(stringsAsFactors = F)   
  install.packages("glmnet")
  install.packages("dplyr")
  install.packages("RColorBrewer")
  install.packages("e1071")
  
  library(glmnet)
  library(dplyr)
  library(RColorBrewer)
  library(e1071)
  }
