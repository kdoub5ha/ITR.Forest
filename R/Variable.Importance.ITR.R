#' @title Calcuates variable importance measures for a random forest object.  
#' 
#' @description This function accepts a forest object from the `Build.RF.ITR` function and estimates the importance of each
#' predictor. This is accomplished by considering each tree in the forest, obtaining the out-of-bag value for each predictor in that tree, 
#' obtaining the permuted out-of-bag value for each predictor in the tree, and comparing the values. A larger discrepancy between 
#' the original value and permuted value indicates the predictor is more important in predicting treatment. The function
#' returns the variables in order of importance along with the importance measure, scaled to be out of 1.  
#' 
#' @param RF.fit forest object from Build.RF.ITR. Required input. 
#' @param n0 minimum number of treatment/control observations needed in a split to call a node terminal. Defaults to 2. 
#' @param sort sort the variable importance measure? Defaults to TRUE. 
#' @param N0 minimum number of observations needed in a split to call a node terminal. Defaults to 20. 
#' @param details print details of each tree as the function progresses. Defaults to FALSE.
#' @param truncate.zeros sets variable importances less than 0 to 0. Defaults to TRUE.
#' @param depth internal variable.
#' @param AIPWE indicator for AIPWE estimation.
#' @return Returns ordered variable importance measure calculated for each splitting variable. 
#' @import randomForest
#' @examples 
#' set.seed(1)
#' dat <- gdataM(n = 1000, depth = 2, beta1 = 3, beta2 = 1)
#' # Build a forest with 100 trees
#' forest <- Build.RF.ITR(dat, col.y="y", col.trt="trt", col.prtx="prtx", split.var=1:4, ntree=100)
#' # Calculate variable importance measures (X1 and X3 should be returned as the most important)
#'  Variable.Importance.ITR(forest)
#' # X1          X3          X4          X2 
#' # 0.727854671 0.260080262 0.009528276 0.002536791 
#' @export


Variable.Importance.ITR <- function(RF.fit, 
                                    n0 = 5, 
                                    N0 = 20, 
                                    sort = TRUE, 
                                    details = FALSE,
                                    truncate.zeros = TRUE,
                                    depth = 1, 
                                    AIPWE = FALSE)
  {
  trees <- RF.fit$TREES
  id.boots <- RF.fit$ID.Boots.Samples
  # ARGUMENTS FOR MODEL SPECIFICATION 
  Model.Specification <- RF.fit$Model.Specification
  dat0 <- Model.Specification$data
  col.y <- Model.Specification$col.y
  col.trt <- Model.Specification$col.trt
  col.prtx <- Model.Specification$col.prtx
  split.var <- Model.Specification$split.var
  ctg <- Model.Specification$ctg
  vnames <- colnames(dat0)[split.var]
  # 
  ntree <- length(trees)
  p <- length(split.var)
  VI <- rep(0, p)
  for (b in 1:ntree){
    id.b <- id.boots[[b]]
    dat.oob <- dat0[-sort(unique(id.b)), ] 
    n.oob <- nrow(dat.oob)	
    tre.b <- trees[[b]]
    ########## NOTE THAT revise.tree=T HERE! ##########
    out0.b <- send.down.VI.ITR(dat.new=dat.oob, tre=tre.b, col.y=col.y, col.trt=col.trt, 
                               col.prtx=col.prtx, ctg=ctg, n0=n0, N0=N0, revise.tree=T,depth=depth,AIPWE = AIPWE)  
    tre0.b <- out0.b$tre0				
    if (nrow(tre0.b) > 0) {					### AVOID NULL TREES	
      Xs.b <- sort(unique(na.omit(tre0.b$var))) 
      G.oob <- out0.b$score
      for (j in 1:p) {
        if (details) print(j)
        G.j <- G.oob
        col.xj <- split.var[j] 
        if (is.element(col.xj, Xs.b)){			
          x.j <- dat.oob[, col.xj]
          dat.permuted <- dat.oob
          dat.permuted[ , col.xj] <- x.j[sample(1:n.oob,n.oob, replace=F)]
          ########## NOTE THAT revise.tree=F HERE! ##########
          out0.bj <- send.down.VI.ITR(dat.new=dat.permuted, tre=tre0.b, col.y=col.y, col.trt=col.trt, 
                                      col.prtx=col.prtx, ctg=ctg, n0=n0, N0=N0, revise.tree=F,depth=1,AIPWE = AIPWE)
          tre0.bj <- out0.bj$tre0		
          G.j <- ifelse(nrow(tre0.bj) ==1, G.oob, out0.bj$score)
        }
        if (G.j > G.oob) G.j <- G.oob  		
        ##################### PREVENTS NEGATIVE IMPORTANCE VALUES 
        VI[j] <- VI[j] + (G.oob - G.j)/G.oob
      }
    }	
  }
  if (truncate.zeros) VI[VI <0] <- 0  		####### IS THIS STEP NECESSARY? NOPE. 
  names(VI) <- vnames
  if (sort) VI <- sort(VI, decreasing=T) 
  VI<-VI/sum(VI)
  return(VI)
}