#' @title Generates partition summary based on itr value. 
#' 
#' @description The partitioning function is used inside of the tree growing functions.  It 
#' selects the best split among a set of predictors based on the IPWE or AIPWE value. 
#' 
#' @param dat data set from which the partition is to be made.  Must contain outcome, binary 
#'  treatment indicator, columns of splitting covariates, and column of probability of being
#'  in treatment group.
#' @param split.var columns of potential spliting variables. Required input.
#' @param min.ndsz minimum number of observations required to call a node terminal. Defaults to 20.
#' @param ctg identifies the categorical input columns.  Defaults to NULL.  Not available yet. 
#' @param n0 minimum number of treatment/control observations needed in a split to call a node terminal. Defaults to 5. 
#' @param max.depth controls the maximum depth of the tree. Defaults to 15. 
#' @param mtry sets the number of randomly selected splitting variables to be included. Defaults to number of splitting variables.
#' @param max.score controls the minimum score required to make an additional split (internally controlled).
#' @param dat.rest data outside the current node being split
#' @param AIPWE indicator for AIPWE estimation
#' @return summary of the best split for a given data frame. 
#' @export



partition.ITR<-function(dat, test=NULL, name="0", min.ndsz=20, n0=5, split.var, ctg=ctg, 
                        max.depth=15, mtry=length(split.var), dat.rest=NULL, max.score=NULL, AIPWE = AIPWE)
{   
  # inialize various variable
  call <- match.call()
  out <- match.call(expand = F)
  out$info <- NULL
  out$name.l <- NULL
  out$name.r <- NULL
  out$left <- NULL
  out$right <- NULL
  out$... <- NULL
  # label the binary tree by 1 (left) and 2 (right).
  name.l <- paste(name, 1, sep="")
  name.r <- paste(name, 2, sep="")
  # sample size
  n <- nrow(dat)
  # check whether testing data is provided
  if (!is.null(test)) {
    n.test <- nrow(test)
    score.test <- NA
  }
  # prepare for the first cut these variable were used to store final cut information
  var <- NA
  vname <- NA
  cut <- NA
  
  if(name=="0") {
    dat.comb<-dat
  } else{
    dat.comb <- rbind(dat.rest[,which(is.element(colnames(dat.rest), c('y', 'trt', 'prtx')))], dat[,which(is.element(colnames(dat), c('y', 'trt', 'prtx')))])
  }
  
  # extract value from data
  trt <- dat$trt
  y <- dat$y
  vnames <- colnames(dat)
  # COMPUTE THE TREATMENT EFFECT IN CURRENT NODE
  trt.effect <- NA
  n.0 = length(y[trt==0])
  n.1 = n - n.0
  if (min(n.1, n.0) >0) {
    trt.effect <- mean(y[trt==1]) - mean(y[trt==0])
  }
  # CONTROL THE MAX TREE DEPTH
  # name is the current tree label.
  # only when currently depth < max.depth and n > min terminal node proceed.
  depth <- nchar(name) 
  if (depth <= max.depth && n >= min.ndsz) {
    if (is.null(mtry)) {
      m.try = length(split.var)
    }else{
      m.try = mtry
    }
    
    # Search across all covariates
    for(i in sample(split.var, size=m.try, replace=F)) {
      x <- dat[,i]
      v.name <- vnames[i]
      temp <- sort(unique(x))
      if(length(temp) > 1) { 
        # handle categorial variable first, otherwise take out the final value as no cut after it.
        if (is.element(i,ctg)){
          zcut <- power.set(temp)
        } else{
          zcut <- temp[-length(temp)]
        }
        # zcut are the values for all possible cuts 
        for(j in zcut) {
          
          score <- NA
          if (is.element(i,ctg)){
            grp.l <- sign(is.element(x, j))
            grp.r <- sign(!is.element(x, j))    #This is modified
            cut1.l <- cbind("l", as.character(j))    #This is modified
            cut1.r <- cbind("r", as.character(j))    #This is modified
          } else  {
            # define left and right groups
            grp.l <- sign(x <= j)
            cut1.l <-  cbind("l",as.character(j))
            grp.r <- sign(x > j)
            cut1.r <- cbind("r",as.character(j))
          }
          
          if(min(sum(grp.l), sum(grp.r)) >= min.ndsz & min(sum(trt*grp.l), sum((1-trt)*grp.l), sum((1-trt)*grp.r), sum(trt*grp.r)) >= n0){
            # use itr rule to calcuate measure of splitting
            #n.1 <- sum(grp.l==1)
            #n.0 <- n-n.1
            if(name=="0") score.l <- itrtest(dat.comb, z=grp.l, n0, AIPWE)
            if(name!="0") score.l <- itrtest(dat.comb, z=c(dat.rest$trt.new, grp.l), n0, AIPWE)
            #n.1 <- sum(grp.r==1)
            #n.0 <- n-n.1
            
            if(name=="0") score.r <- itrtest(dat.comb, z=1-grp.l, n0, AIPWE)
            if(name!="0") score.r <- itrtest(dat.comb, z=c(dat.rest$trt.new, 1-grp.l), n0, AIPWE)
            # record the one with improved utility
            if (!is.na(score.l) && !is.na(score.r)) {
              if(score.l>max.score & score.r>max.score){
                if(score.l>score.r) {
                  max.score <- score.l
                  var <- i
                  vname <- v.name
                  cut <- cut1.l
                  best.cut<-j
                }else{
                  max.score <- score.r
                  var <- i
                  vname <- v.name
                  cut <- cut1.r
                  best.cut<-j
                }
                
              }else if(score.l>max.score & score.r<max.score){
                max.score <- score.l
                var <- i
                vname <- v.name
                cut <- cut1.l
                best.cut<-j
              }else if(score.l<max.score & score.r>max.score){
                max.score <- score.r
                var <- i
                vname <- v.name
                cut <- cut1.r
                best.cut<-j
              }
            }
          }
        }
      }
    }
  }
  # when testing data is provided, assess new treatment assignment 
  # using testing sample and the rule calculated from training sample
  # var is the covariates calcualted before where spliting adopts. 
  # best.cut is the cutting point.
  if (!is.null(test)) {
    n.test <- nrow(test)
    score.test <- NA
    if (!is.na(var)) {
      if (is.element(var,ctg)) {
        grp.test <- sign(is.element(test[,var], best.cut))
      } else  {
        if(cut[1]=="l"){
          grp.test <- sign(test[,var] <= best.cut)
        } else{
          grp.test <- sign(test[,var] > best.cut)
        }
      }
      score.test <- itrtest(test, z=grp.test, n0=n0, aug = AIPWE)
      if (!is.na(score.test)){
        out$name.l <- name.l
        out$name.r <- name.r
        if(cut[1]=="l"){
          out$left.test <- test[grp.test==1,  ]
          out$right.test <- test[grp.test==0,  ]
        } else{
          out$left.test <- test[grp.test==0,  ]
          out$right.test <- test[grp.test==1,  ]
        }
        if (is.element(var,ctg)) {
          out$left  <- dat[is.element(dat[,var], best.cut),]
          out$right <- dat[!is.element(dat[,var], best.cut), ]}
        else {
          if(cut[1]=='l'){
            out$left  <- cbind(dat[dat[,var]<= best.cut,],new.trt=rep(1,n=sum(dat[,var]<= best.cut)))
            out$right <- cbind(dat[dat[,var]> best.cut, ],new.trt=rep(0,n=sum(dat[,var]> best.cut)))
          }else{
            out$left  <- cbind(dat[dat[,var]<= best.cut,],new.trt=rep(0,n=sum(dat[,var]<= best.cut)))
            out$right <- cbind(dat[dat[,var]> best.cut, ],new.trt=rep(1,n=sum(dat[,var]> best.cut)))
          }  
        }
      } else {
        var <- NA
        vname <- NA
        cut <- NA
        max.score <- NA
      }
      # output results from both testing and training data.
      if(!is.na(var)){  out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect,var = var, 
                                               vname=vname, cut.1 = cut[1], cut.2 = cut[2], score=ifelse(max.score==-1e20, NA, max.score), score.test=score.test, size.test=n.test)
      } else{
        out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect, var = NA, 
                               vname=NA, cut.1 = NA, cut.2 = NA, score=NA,score.test=NA, size.test=n.test)
      }
    } else {
      out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect, var = NA, 
                             vname=NA, cut.1=NA, cut.2=NA, score=NA,score.test=NA, size.test=n.test)
    }
  } else {
    # if no testing data output results from training data only.
    if (!is.na(var)) {
      out$name.l <- name.l
      out$name.r <- name.r
      if (is.element(var,ctg)) {                                                                               
        out$left  <- dat[is.element(dat[,var], best.cut),]
        out$right <- dat[!is.element(dat[,var], best.cut), ]
      } else {
        if(cut[1]=='l'){
          out$left  <- cbind(dat[dat[,var]<= best.cut,],new.trt=rep(1,n=sum(dat[,var]<= best.cut)))
          out$right <- cbind(dat[dat[,var]> best.cut, ],new.trt=rep(0,n=sum(dat[,var]> best.cut)))
        }else{
          out$left  <- cbind(dat[dat[,var]<= best.cut,],new.trt=rep(0,n=sum(dat[,var]<= best.cut)))
          out$right <- cbind(dat[dat[,var]> best.cut, ],new.trt=rep(1,n=sum(dat[,var]> best.cut)))
        }  
      }
      out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect, var = var, 
                             vname=vname, cut.1 = unique(cut[,1]), cut.2 = paste(cut[,2], collapse = ','), score=ifelse(max.score==-1e20, NA, max.score))
    } else {
      out$info <- data.frame(node=name, size = n, n.1=n.1, n.0=n.0, trt.effect=trt.effect,var=NA, 
                             vname=NA, cut.1= NA,cut.2=NA, score=NA)
    }
  }
  out 
}