#' @title Grows a large interaction tree
#' 
#' @description This function grows an interaction tree using either the IPWE or AIPWE method (AIPWE=F, T). 
#' 
#' 
#' @param data data set from which the tree is to be grown.  Must contain outcome, binary 
#'  treatment indicator, columns of splitting covariates, and column of probability of being
#'  in treatment group.
#' @param test testing data
#' @param split.var columns of potential spliting variables. Required input.
#' @param min.ndsz minimum number of observations required to call a node terminal. Defaults to 20.
#' @param ctg identifies the categorical input columns.  Defaults to NULL.  Not available yet. 
#' @param n0 minimum number of treatment/control observations needed in a split to call a node terminal. Defaults to 5. 
#' @param max.depth controls the maximum depth of the tree. Defaults to 15. 
#' @param mtry sets the number of randomly selected splitting variables to be included. Defaults to number of splitting variables.
#' @param AIPWE logical. Should AIPWE (TRUE) or IPWE (FALSE) be used. 
#' @param in.forest logical for if the tree is being constructed in a forest. Should not be changed from defaults.
#' @param stabilize.type gives the method used for calculating residuals. Current options are 'rf' for random forest and 'linear' for linear model. 
#' @return Summary of a single interaction tree. Each `node` begins with "0" indicating the root node, 
#' followed by a "1" or "2" indicating the less than (or left) child node or greater than (or right) child node. 
#' Additionally, the number of observations `size`, number treated `n.1`, number on control `n.0`, and treatment effect `trt.effect`
#' summaries are provided.  The splitting information includes the column of the chosen splitting variable `var`, the variable name 'vname',
#' the direction the treatment is sent `cut.1` ("r" for right child node, and "l" for left), the chosen split value `cut.2`, 
#' and the estimated value function `score`.
#' @import randomForest
#' @export
#' @examples
#' dat <- gdataM(n=1000, depth=2, beta1=3, beta2=1)
#' # Generates tree using simualated EMR data with splitting variables located in columns 1-4.
#' tree <- grow.ITR(data=dat, split.var=1:4)


grow.ITR <- function(data, 
                     test=NULL, 
                     min.ndsz=20, 
                     n0=5, 
                     split.var, 
                     ctg=NULL, 
                     max.depth=15, 
                     mtry=length(split.var), 
                     AIPWE=F, 
                     in.forest = FALSE, 
                     stabilize.type = 'rf')
{
  # initialize variables.
  out <- NULL
  list.nd <- NULL
  list.test <- NULL
  temp.list <- NULL
  temp.test <- NULL
  temp.name <- NULL
  # record total dataset for spliting 
  if(!in.forest){
    if(stabilize.type == "rf"){
      # Replace raw measures with residuals
      fit <- randomForest(y = data$y, x = as.data.frame(data[,split.var]))
      resids <- fit$y - fit$predicted
      data$y <- resids
    }
    if(stabilize.type == "linear"){
      fit <- lm(data$y~as.matrix(data[ , c(split.var, which(colnames(data) == "trt")) ]))
      data$y <- fit$residuals
    }
  }
  list.nd <- list(data)
  if (!is.null(test)) list.test <- list(test)
  name <- "0"
  full.set <- data
  max.score <- NULL
  # loop over dataset for spliting 
  while (length(list.nd)!=0) {
    for (i in sample(1:length(list.nd), size = length(list.nd), replace = F)){
      if (!is.null(dim(list.nd[[i]])) && nrow(list.nd[[i]]) > 1){
        test0 <- NULL
        if (!is.null(test)) test0 <- list.test[[i]]
        if(length(list.nd)>1) {temp.tree=temp.tree} else{temp.tree <- NULL}
        
        # inilize score statistics. 
        # at the initial stage, there is no subgroup. 
        # Inidividuals either assign to trt=1 (z=rep(1,dim(dat)[1]))  
        # or trt=0 (z=rep(0,dim(dat)[1])) depending on which one gives better utility.
        # This sets the default utility as the max of all patients in treatment or all in control
        if(name[i]=="0"){
          max.score <- max(itrtest(data, z=rep(0,nrow(data)), n0, AIPWE),itrtest(data, z=rep(1,nrow(data)), n0, AIPWE))
        }else{
          send<-send.down(data, temp.tree, ctgs = ctg)
          node<-substr(send$data$node,1,nchar(send$data$node)-1)
          direction<-substr(send$data$node,nchar(send$data$node),nchar(send$data$node))
          trt.dir <- temp.tree[match(node,temp.tree$node),]$cut.1
          
          trt.pred<-ifelse(trt.dir=="r" & direction=="1",0,
                           ifelse(trt.dir=="r" & direction=="2",1,
                                  ifelse(trt.dir=="l" & direction=="1",1,0)))
          
          max.score <- itrtest(dat = data, z=trt.pred, n0, AIPWE)
          
          #Obtain treatments for those not included in the node
          dat.rest <- data[!is.element(data$id, list.nd[[i]]$id),]
          dat.rest$trt.new <- trt.pred[!is.element(data$id, list.nd[[i]]$id)]
        }
        
        # Determine best split across all covariates
        split <- partition.ITR(dat = list.nd[[i]], test = test0, name = name[i], min.ndsz=min.ndsz,n0=n0, split.var=split.var, ctg=ctg, 
                               max.depth=max.depth, mtry=mtry, dat.rest=dat.rest, max.score=max.score, AIPWE = AIPWE)
        out <- rbind(out, split$info)
        if(!is.null(nrow(split$left))&&!is.null(nrow(split$right))){
          min.n <- min(nrow(split$left),nrow(split$right))
        }
        if (!is.null(split$left) && min.n>min.ndsz && is.null(test)) {
          temp.list <- c(temp.list, list(split$left, split$right))
          temp.test <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
        } else if (!is.null(split$left) && min.n>min.ndsz && !is.null(test) && !is.null(split$left.test)) {
          temp.list <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
          temp.test <- c(temp.test, list(split$left.test, split$right.test))
        }
      }
      temp.tree <- out
    }
    list.nd <- temp.list
    list.test <- temp.test
    name <- temp.name
    temp.tree<-out
    temp.list <- NULL
    temp.test <- NULL
    temp.name <- NULL
  } 
  out$node <- as.character(out$node)
  out <- out[order(out$node), ]
  for(p in 1:dim(out)[1]) {out$var[p]<-ifelse(is.na(de(out$node[p],out)[1]),NA,out$var[p])}
  for(p in 1:nrow(out)) {if(is.na(out$var[p])) out$vname[p] <- out$cut.1[p] <- out$cut.2[p] <- out$score[p] <- NA}
  if(!is.null(test)){
    for(p in 1:nrow(out)) {out$score.test[p]<-ifelse(is.na(out$var[p]),NA,out$score.test[p])}
  }
  if(!is.null(test)){
    pruned<-prune(tre = out, a=0, train=data, test=test, AIPWE = AIPWE, ctgs = ctg)
    out$score.test <- as.numeric(pruned$Va.test[match(out$node, pruned$node.rm)])
  }
  out
}