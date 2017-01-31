#' Sends testing data down a tree to assess the performance of the tree. This is used inside
#' the variable importance and pruning functions. 
#' 
#' @param dat data to be run down the tree.  Required input. 
#' @param tre tree object from grow.ITR().  Required input.
#' @return summary of tree performance
#' @export

send.down <- function(data, tre, char.var=1000)
{
  call <- match.call()
  out <- match.call(expand = F)
  out$tree <- out$data <- out$... <- NULL
  dat <- cbind(data, node=0)
  tre <- cbind(tre, n.test=NA)
  cut.point <- as.vector(tre$cut.2) 
  split.var <- as.numeric(as.vector(tre$var)) 
  for (i in 1:nrow(tre)){
    in.node <- (dat$node)==(tre$node[i])
    tre$n.test[i] <- sum(in.node)                       
    if (!is.na(split.var[i])){
      # print(cbind(i, var=tre$var[i], cut=tre$cut[i]))
      var.split <- dat[,split.var[i]] 
      cut <- cut.point[i]
      if (!is.element(split.var[i], char.var)) { 
        cut1 <- as.numeric(cut)    
        l.nd <- dat$node[in.node & var.split <= cut1] 
        r.nd <- dat$node[in.node & var.split > cut1]
        dat$node[in.node & var.split <= cut1] <- paste(l.nd, 1, sep="")
        dat$node[in.node & var.split >  cut1] <- paste(r.nd, 2, sep="")  
      }
      else {
        var.split <- as.character(var.split)
        cut1 <- unlist(strsplit(as.character(cut), split=" ")) #####################
        l.nd <- dat$node[in.node & is.element(var.split, cut1)] 
        r.nd <- dat$node[in.node & !is.element(var.split, cut1)]                  
        dat$node[in.node & is.element(var.split, cut1)] <- paste(l.nd, 1, sep="")  
        dat$node[in.node & !is.element(var.split, cut1)] <- paste(r.nd, 2, sep="")}                   
    }
  }
  # print(tre)
  out$data <- dat
  out$tree <- tre
  out 
}
