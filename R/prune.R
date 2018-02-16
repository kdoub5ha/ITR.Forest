#' @title Prunes a tree for a given penalty value
#' 
#' @description The `prune` function allows the user to specify a value of the penalty for a given tree. 
#' This function uses the "weakest link" criteria in order to evaluate the order in which branches are pruned
#' and gives the penalized value along with the unpenalized value. If testing data are provided for validation, 
#' then the penalized and unpenalized values from the testing data run down the tree structure are also provided.  
#' 
#' @param tre sets the tree to be pruned 
#' @param a sets the value of the splitting penalty
#' @param train the training data used to create the tree
#' @param test the testing data to be used.  Defaults to NULL.
#' @param AIPWE indicator for AIPWE estimation.
#' @param ctgs columns of categorical variables.
#' @return summary of pruned branches and the associated value of the tree after pruning. 
#' @return \item{result}{contains columns: `node.rm` which is the weakest link at each
#' iteration of the pruning algorithm; `size.tree` and `n.tmnl` give number of total nodes and
#' terminal nodes in each subtree; `alpha` is the weakest link scoring criteria; `V` and `V.test`
#' are the overall value of the tree for the training and tesing samples; `V.a` and `Va.test`
#' give the penalized value for the training and testing samples.}
#' @export
#' 


prune <- function(tre, a, train, test=NULL, AIPWE = F, n0=5, ctgs = NULL){
  if(is.null(dim(tre))) stop("No Need to Prune Further.")
  result <- NULL
  n.tmnl <- sum(is.na(tre$var))
  subtree <- 1
  while (n.tmnl > 1 ) {
    #internal keeps track of all splits which are not terminal <NA> for score value
    internal <- tre$node[!is.na(tre$cut.1)]
    l <- length(internal)
    #r.value is the vector of mean score values across all splits
    r.value <- 1:l
    for(j in 1:l) {
      #branch keeps track of all splits (terminal or not)
      #branch is a single path which can be followed down a given tree
      sub.tree <- tre[!is.element(tre$node,de(internal[j], tree=tre)),]
      if(nrow(sub.tree)>1){
        if(!is.null(test)) sub.tree[sub.tree$node==internal[j],][6:11] <- NA
        if(is.null(test)) sub.tree[sub.tree$node==internal[j],][6:10] <- NA
        
        trt.pred <- predict.ITR(sub.tree, train, ctgs = ctgs)$trt.pred
        
        score <- itrtest(dat = train, z=trt.pred, n0=n0, AIPWE)
      }else{
        score <- max(itrtest(train, rep(0,nrow(train)), 0, AIPWE), itrtest(train, rep(1,nrow(train)), 0, AIPWE))
      }
      #r.value is a penalized mean value score across all internal nodes
      r.value[j] <- score #dim(branch)[1] - nchar(internal[j]) + score[1]/1000
    } 
    
    alpha <-max(r.value, na.rm = T)
    nod.rm <- internal[r.value==alpha]
    # if (length(nod.rm)>1) print("Multiple Nodes will be pruned. Check!")
    trt.pred <- predict.ITR(tre, train, ctgs = ctgs)$trt.pred
    
    V <- itrtest(dat = train, z=trt.pred, n0=n0, AIPWE)
    V.a <- V - a*sum(!is.na(tre$score))
    
    if(!is.null(test)){
      # Calculate value for the training set
      trt.pred <- predict.ITR(tre, test, ctgs = ctgs)$trt.pred
      V.test <- itrtest(dat = test, z=trt.pred, n0=-1, AIPWE)
      Va.test <- V.test - a*sum(!is.na(tre$score))
    }
    
    # Calculate value for testing data
    if(is.null(test)) {result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tree=nrow(tre), 
                                                     size.tmnl=nrow(tre)-l, alpha=alpha, V=V, V.a=V.a, V.test=NA, Va.test=NA))}
    if(!is.null(test)) {result <- rbind(result, cbind(subtree=subtree, node.rm=nod.rm, size.tree=nrow(tre), 
                                                      size.tmnl=nrow(tre)-l, alpha=alpha, V=V, V.a=V.a, V.test=V.test, Va.test=Va.test))}
    
    if(length(nod.rm)>1){
      for(k in 1:length(nod.rm)){
        tre <- tre[!is.element(tre$node, de(nod.rm[k],tre)),]
        if(is.null(test)) tre[match(nod.rm[k], tre$node), c("var", "vname", "cut.1", "cut.2", "score")] <- NA
        if(!is.null(test))  tre[match(nod.rm[k], tre$node), c("var", "vname", "cut.1", "cut.2", "score", "score.test")] <- NA
        n.tmnl <- sum(is.na(tre$cut.1))
        subtree <- subtree + 1  
      }
    } else{
      tre <- tre[!is.element(tre$node, de(nod.rm,tre)),]
      if(is.null(test)) tre[match(nod.rm, tre$node), c("var", "vname", "cut.1", "cut.2", "score")] <- NA
      if(!is.null(test))  tre[match(nod.rm, tre$node), c("var", "vname", "cut.1", "cut.2", "score", "score.test")] <- NA
      n.tmnl <- sum(is.na(tre$cut.1))
      subtree <- subtree + 1
    }
  }
  # HANDLE THE NULL TREE WITH THE ROOT NODE ONLY
  if(!is.null(test)){
    result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre), 
                                  size.tmnl=1, alpha=9999, 
                                  V=max(itrtest(train, rep(1,nrow(train)), 5, AIPWE), itrtest(train, rep(0, nrow(train)), 5, AIPWE)), 
                                  V.a=max(itrtest(train, rep(1,nrow(train)), 5, AIPWE), itrtest(train, rep(0, nrow(train)), 5, AIPWE)), 
                                  V.test=max(itrtest(test, rep(1,nrow(test)), 5, AIPWE), itrtest(test, rep(0, nrow(test)), 5, AIPWE)), 
                                  Va.test=max(itrtest(test, rep(1,nrow(test)), 5, AIPWE), itrtest(test, rep(0, nrow(test)), 5, AIPWE))))
  }else{
    result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre), 
                                  size.tmnl=1, alpha=9999, 
                                  V=max(itrtest(train, rep(1,nrow(train)), 5, AIPWE), itrtest(train, rep(0, nrow(train)), 5, AIPWE)), 
                                  V.a=max(itrtest(train, rep(1,nrow(train)), 5, AIPWE), itrtest(train, rep(0, nrow(train)), 5, AIPWE)), 
                                  V.test=NA, Va.test=NA))    
  }    
  result <- as.data.frame(result)
  result
}