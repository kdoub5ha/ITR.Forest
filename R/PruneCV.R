#' This function is outdated and it is recommended that users use treeCV instead.
#' Performs k-fold cross validation for selection optimal tuning parameter lambda. 
#' 
#' @param data input data 
#' @param param vector of pruning parameter values to be considered
#' @param min.ndsz sets the minimal number of observations allowed in terminal nodes. Defaults to 20. 
#' @param n0 sets the minimum number of observations from each treatment group to be in each child node.  Defaults to 5. 
#' @param split.var specifies the columns of splitting variables in the input data. 
#' @param nfolds defaults to 10.
#' @param graphic logical indicator for whether CV graphic should be generated.  Defaults to FALSE.
#' @param graphicname filename of the CV graphic.  Must be specified if graphic=TRUE.  Defaults to FALSE.  
#' @return summary of pruned branches and the associated value of the tree after pruning
#' @export
#' @examples
#' 


PruneCV<-function(data, param, min.ndsz=20, n0=5, split.var, nfolds=10, graphic=FALSE, graphicname=NULL){
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
    if(is.null(graphicname)) stop("Must specify name for output graphic")
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