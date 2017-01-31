#' Grows a large interaction tree
#' 
#' @param data data set from which the tree is to be grown.  Must contain outcome, binary 
#'  treatment indicator, columns of splitting covariates, and column of probability of being
#'  in treatment group.
#' @param test provided testing data.  Defaults to NULL. 
#' @param split.var columns of potential spliting variables. Required input.
#' @param min.ndsz minimum number of observations required to call a node terminal. Defaults to 20.
#' @param ctg identifies the categorical input columns.  Defaults to NULL.  Not available yet. 
#' @param N0 minimum number of observations needed to call a node terminal.  Defaults to 20. 
#' @param n0 minimum number of treatment/control observations needed in a split to call a node terminal. Defaults to 5. 
#' @param max.depth controls the maximum depth of the tree. Defaults to 15. 
#' @param mtry sets the number of randomly selected splitting variables to be included. Defaults to number of splitting variables.
#' @return summary of single interaction tree 
#' @export
#' @examples
#' tree<-grow.ITR(data=rctdata, split.var=3:7)
#' Generates tree using rctdata with potential splitting variables located in columns 3-7.

grow.ITR<-function(data, test=NULL, min.ndsz=20, n0=5, split.var, ctg=NULL, max.depth=15, mtry=length(split.var))
{
  # initialize variables.
  out <- NULL
  list.nd <- NULL
  list.test <- NULL
  temp.list <- NULL
  temp.test <- NULL
  temp.name <- NULL
  # record total dataset for spliting 
  list.nd <- list(data)
  if (!is.null(test)) list.test <- list(test)
  name <- "0"
  max.score <- NULL
  full <- data
  # loop over dataset for spliting 
  while (length(list.nd)!=0) {
    for (i in 1:length(list.nd)){
      if (!is.null(dim(list.nd[[i]])) && nrow(list.nd[[i]]) > 1){
        test0 <- NULL
        if (!is.null(test)) test0 <- list.test[[i]]
        if(length(list.nd)>1) {temp.tree=temp.tree} else{temp.tree <- NULL}
        # Determine best split across all covariates
        split <- partition.ITR(dat = list.nd[[i]], test = test0, name = name[i], min.ndsz=min.ndsz,n0=n0, split.var=split.var, ctg=ctg, max.depth=max.depth, mtry=mtry, temp.tree=temp.tree, full=full)
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
  for(p in 1:dim(out)[1]) {out$vname[p]<-ifelse(is.na(de(out$node[p],out)[1]),NA,out$vname[p])}
  for(p in 1:dim(out)[1]) {out$cut.1[p]<-ifelse(is.na(de(out$node[p],out)[1]),NA,out$cut.1[p])}
  for(p in 1:dim(out)[1]) {out$cut.2[p]<-ifelse(is.na(de(out$node[p],out)[1]),NA,out$cut.2[p])}
  for(p in 1:dim(out)[1]) {out$score[p]<-ifelse(is.na(de(out$node[p],out)[1]),NA,out$score[p])}
  if(!is.null(test)){
    for(p in 1:dim(out)[1]) {out$score.test[p]<-ifelse(is.na(out$var[p]),NA,out$score.test[p])}
  }
  if(!is.null(test)){
  pruned<-prune(tre = out, a=0, train=data, test=test)
  out$score.test <- as.numeric(pruned$Va.test[match(out$node, pruned$node.rm)])
  }
  out
}
