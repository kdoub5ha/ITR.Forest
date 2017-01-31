#' Performs k-fold cross validation for choosing correctly sized tree. 
#' 
#' @param data data set used to create trees.  
#' @param param vector of pruning parameters to be considered.  
#' @param min.ndsz minimum number of obervations allowed in a terminal node.  Defaults to 20. 
#' @param n0 minimum number of observations from each treatment group allowed in a terminal node.  Defaults to 5.
#' @param split.var specifies columns of splitting variable in the input data set.
#' @param graphic indicates whether a cross validation plot is generated.  Defaults to FALSE.
#' @param graphicname name of output graphic when graphic=TRUE.  Must be specified if graphic=TRUE and will be output to working directory.
#' @param nfolds number of folds to be considered in the cross validation.  Defaults to 10. 
#' @return summary of pruned branches and the associated value of the tree after pruning
#' @export
#' @examples
#' 


PruneCV<-function(data, param, min.ndsz=20, n0=5, split.var, graphic=FALSE, nfolds=10, graphicname=NULL){
  tree<-prune.tree<-out<-NULL
  
  pb <- txtProgressBar(min = 0, max = length(param), style = 3)
  
  for(x in 1:length(param)){
    data<-data[sample(1:nrow(data),nrow(data)),]
    folds <- cut(seq(1,nrow(data)),breaks=nfolds,labels=FALSE)
    comb <- NULL
    for(i in 1:nfolds){
      #Segement data
      testIndexes <- which(folds==i,arr.ind=TRUE)
      test <- data[unique(testIndexes), ]
      train <- data[-unique(testIndexes), ]
      
      tree<-grow.ITR(data = train, test=test, split.var=split.var, min.ndsz = min.ndsz, n0=n0)
      
      prune.tree<-prune(tree, param[x], train=train, test=test)
      comb <- rbind(comb, prune.tree)
      
    }
    #obtain aggregated summaries
    col1 <- aggregate(as.numeric(comb[,9]), FUN=mean, by=list(comb[,4]), na.rm=T)
    col2 <- aggregate(as.numeric(comb[,9]), FUN=sd, by=list(comb[,4]), na.rm=T)
    col3 <- aggregate(as.numeric(comb[,9]), FUN=length, by=list(comb[,4]))
    
    #merge them together
    opt<-Reduce(function(x, y) merge(x, y, by="Group.1"), list(col1, col2, col3))
    opt$Parameter <- as.factor(param[x]) 
    
    #generate the output data frame
    out <- rbind(out, opt)
    setTxtProgressBar(pb, value=x)
  }
  close(pb)
  out<-data.frame(out)
  names(out) <- c("n.tmnl", "MeanValue", "SDValue","Length", "Parameter")
  #plot graphic if indicated
  if((graphic)){
    postscript(file = graphicname)
    p <- ggplot(data=out[which(out$Length>4),], aes(x=as.numeric(n.tmnl), y=log(MeanValue), group=Parameter))
    p <- p + geom_point() + geom_line(aes(linetype=Parameter))
    p <- p + xlab("Number of Terminal Nodes") + ylab(expression(log~V[lambda](T)))
    p <- p + theme(axis.title.x=element_text(size=24)) + theme(axis.title.y=element_text(size=24))
    p <- p + theme(axis.text.x=element_text(size=20)) + theme(axis.text.y=element_text(size=20))
    p <- p + scale_x_continuous(breaks=1:max(out$n.tmnl)) + theme(legend.position=c(0.9,0.5))
    p <- p + theme(legend.title=element_text(size=18) , legend.text=element_text(size=14))
    p <- p + theme(panel.background=element_blank()) + theme(axis.line = element_line(colour = "black"))
    print(p)
    dev.off()
  }
  
  return(out)
}
