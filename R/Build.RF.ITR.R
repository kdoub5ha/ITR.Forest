#' Builds the random forest of interaction trees
#' 
#' @param dat the data set being used to grow the random forest. Required input. 
#' @param test indicates if testing data for each tree is to be run down the tree and included in the RF output.  Defaults to FALSE. 
#' @param col.y the response variable. Required input. 
#' @param col.trt the treatment indicator.  Must be binary. Required input.
#' @param col.prtx the probability of being assigned to treatment group. Required input. 
#' @param split.var vector of columns containing the desired splitting variables.  Required input. 
#' @param ctg identifies the categorical input columns.  Defaults to NA.  Not available yet. 
#' @param N0 minimum number of observations needed to call a node terminal.  Defaults to 20. 
#' @param n0 minimum number of treatment/control observations needed in a split to call a node terminal. Defaults to 5. 
#' @param max.depth controls the maximum depth of the tree. Defaults to 10. 
#' @param mtry sets the number of randomly selected splitting variables to be included. Defaults to max of length(split.var)/3 rounded down and 1.
#' @param ntree sets the number of trees to be generated. Defaults to 500.
#' @param avoid.nul.tree controls if trees with no splits (null trees) are allowed. Defaults to FALSE.
#' @param AIPWE logical for use of augmented robust estimator
#' @return summary of randomly generated trees (summary done by tree)
#' @export
#' @examples
#' forest<-Build.RF.ITR(dat=data, col.y="y", col.trt="trt", col.prtx="prtx", 
#'                      split.var=3:7)
#' This builds a forest of 500 trees using the dataset called 'data' with columns
#' 'y', 'trt', and 'prtx' for the outcome, treatement indicator, and probability of being
#' in treatment group, respectively.  The splitting variables are found in columns 3-7.


Build.RF.ITR<-function(dat, test=FALSE, col.y, col.trt, col.prtx=col.prtx, split.var, ctg=NULL,N0=20, n0=5,  max.depth=10,ntree=500, 
         mtry = max(floor(length(split.var)/3), 1),avoid.nul.tree=F, AIPWE=F)
{
  out <- as.list(NULL)
  out$ID.Boots.Samples  <- as.list(1:ntree)
  out$TREES <- as.list(1:ntree)
  b <- 1
  while (b <= ntree) {
    # TAKE BOOTSTRAP SAMPLES
    id.b <- sample(1:nrow(dat), size=nrow(dat), replace = T)
    dat.b <- dat[unique(id.b),]
    dat.test <- dat[-unique(id.b),]
    # Generate tree based on b-th bootstrap sample
    if(test==FALSE){
      tre.b <- grow.ITR(data=dat.b, test=NULL, min.ndsz=N0, n0=5, split.var=split.var, ctg=NULL, max.depth=15, mtry=mtry, AIPWE=AIPWE)
    } else {
      tre.b <- grow.ITR(data=dat.b, test=dat.test, min.ndsz=N0, n0=5, split.var=split.var, ctg=NULL, max.depth=15, mtry=mtry, AIPWE=AIPWE)
    } 
    if (avoid.nul.tree) {
      if (nrow(tre.b) > 1) {
        out$ID.Boots.Samples[[b]] <- id.b
        out$TREES[[b]] <- tre.b; 
        b <- b +1			
      }
    }
    else {
      out$ID.Boots.Samples[[b]] <- id.b
      out$TREES[[b]] <- tre.b; 
      b <- b +1		
    }
  }
  
  Model.Specification <- as.list(NULL)
  Model.Specification$data <- dat
  Model.Specification$split.var <- split.var
  Model.Specification$ctg <- ctg
  Model.Specification$col.y <- col.y
  Model.Specification$col.trt <- col.trt
  Model.Specification$col.prtx <- col.prtx
  out$Model.Specification <- Model.Specification
  return(out)
}
