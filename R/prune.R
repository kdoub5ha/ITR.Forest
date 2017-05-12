#' Prunes a bootstrap tree with penalty (a) for additional splitting
#' 
#' @param tre sets the tree to be pruned 
#' @param a sets the value of the splitting penalty
#' @param train the training data used to create the tree
#' @param test the testing data to be used.  Defaults to NULL.
#' @param AIPWE logical indicates use of the robust augmented estimator
#' @return summary of pruned branches and the associated value of the tree after pruning
#' @export
#' @examples
#' 

prune <- function(tre, a, train, test=NULL, AIPWE = F){
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
      sub.tree[sub.tree$node==internal[j],][6:11] <- NA
      
      send<-send.down(data = train, tre = sub.tree)
      node<-substr(send$data$node,1,nchar(send$data$node)-1)
      direction<-substr(send$data$node,nchar(send$data$node),nchar(send$data$node))
      trt.dir<-sub.tree[match(node, sub.tree$node),]$cut.1
      
      trt.pred<-ifelse(trt.dir=="r" & direction=="1",0,
                       ifelse(trt.dir=="r" & direction=="2",1,
                              ifelse(trt.dir=="l" & direction=="1",1,0)))
      
      
      score <- itrtest(dat = train, z=trt.pred, n0=2, AIPWE)
      
      #r.value is a penalized mean value score across all internal nodes
      r.value[j] <- score #dim(branch)[1] - nchar(internal[j]) + score[1]/1000
    } 
    
    alpha <-max(r.value, na.rm = T)
    nod.rm <- internal[r.value==alpha]
    # if (length(nod.rm)>1) print("Multiple Nodes will be pruned. Check!")
    send<-send.down(data = train, tre = tre)
    node<-substr(send$data$node,1,nchar(send$data$node)-1)
    direction<-substr(send$data$node,nchar(send$data$node),nchar(send$data$node))
    trt.dir<-tre[match(node,tre$node),]$cut.1
    
    trt.pred<-ifelse(trt.dir=="r" & direction=="1",0,
                     ifelse(trt.dir=="r" & direction=="2",1,
                            ifelse(trt.dir=="l" & direction=="1",1,0)))
    
    
    V <- itrtest(dat = train, z=trt.pred, n0=2, AIPWE)
    V.a <- V - a*sum(!is.na(tre$score))

    if(!is.null(test)){
    #Calculate value for the training set
    send<-send.down(data = test, tre = tre)
    node<-substr(send$data$node,1,nchar(send$data$node)-1)
    direction<-substr(send$data$node,nchar(send$data$node),nchar(send$data$node))
    trt.dir<-tre[match(node,tre$node),]$cut.1
    
    trt.pred<-ifelse(trt.dir=="r" & direction=="1",0,
                     ifelse(trt.dir=="r" & direction=="2",1,
                     ifelse(trt.dir=="l" & direction=="1",1,0)))
      
    V.test <- itrtest(dat = test, z=trt.pred, n0=2, AIPWE)
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
  result <- rbind(result, cbind(subtree=subtree, node.rm='NA', size.tree=nrow(tre), 
                                size.tmnl=1, alpha=9999, V=0, V.a=0, V.test=0, Va.test=0))    
  result <- as.data.frame(result)
  result
}
