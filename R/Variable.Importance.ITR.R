#' Calcuates variable importance measures for a random forest object.  Input must be an object
#'  from the random forest function Build.RF.ITR. 
#' 
#' @param RF.fit forest object from Build.RF.ITR. Required input. 
#' @param n0 minimum number of treatment/control observations needed in a split to call a node terminal. Defaults to 2. 
#' @param sort sort the variable importance measure? Defaults to TRUE. 
#' @param n0 minimum number of treatment/control observations needed in a split to call a node terminal. Defaults to 5. 
#' @param details print details. Defaults to FALSE.
#' @param truncate.zeros sets variable importances less than 0 to 0. Defaults to TRUE.
#' @param AIPWE logical for use of augmented robust estimator
#' @return summary of tree performance
#' @export


Variable.Importance.ITR<-function(RF.fit, n0=5, N0=20, sort=T, details=F, truncate.zeros=T,depth=1, AIPWE = F){
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
