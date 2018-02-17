#' @title Sends testing data down a tree to obtain terminal node assignments 
#' 
#' @description Sends dat.new down tree 'tre' to obtain node assignment. 
#' @param dat.new data to be run down the tree.  Required input. 
#' @param tre tree object from grow.ITR().  Required input.
#' @param ctgs categorical variables, entered as columns in `dat.new`
#' @param char.var internal variable.
#' @return \item{data}{input data with extra column of node assignments}
#' @return \item{tree}{input tree with extra column for number of observations in each node}
#' @export



send.down <- function(dat.new, 
                      tre, 
                      char.var = 1000, 
                      ctgs = NULL)
{
  call <- match.call()
  out <- match.call(expand = F)
  out$tree <- out$data <- out$... <- NULL
  dats <- cbind(dat.new, node=0)
  tre.new <- cbind(tre, n.test=NA)
  cut.point <- as.vector(tre$cut.2) 
  split.v <- as.numeric(as.vector(tre$var)) 
  for (i in 1:nrow(tre)){
    in.node <- (dats$node)==(tre.new$node[i])
    tre.new$n.test[i] <- sum(in.node)                       
    if (!is.na(split.v[i])){
      # print(cbind(i, var=tre$var[i], cut=tre$cut[i]))
      var.split <- dats[,split.v[i]] 
      cut <- cut.point[i]
      if (!is.element(split.v[i], char.var)) {
        if(is.element(split.v[i], ctgs)){
          cut1 <- as.character(strsplit(cut, split = ",")[[1]])
          l.nd <- dats$node[in.node & is.element(var.split, cut1)]
          r.nd <- dats$node[in.node & !is.element(var.split, cut1)]
          dats$node[in.node & is.element(var.split, cut1)] <- paste(l.nd, 1, sep="")
          dats$node[in.node & !is.element(var.split, cut1)] <- paste(r.nd, 2, sep="") 
        } else{
          cut1 <- as.numeric(cut)    
          l.nd <- dats$node[in.node & var.split <= cut1] 
          r.nd <- dats$node[in.node & var.split > cut1]
          dats$node[in.node & var.split <= cut1] <- paste(l.nd, 1, sep="")
          dats$node[in.node & var.split >  cut1] <- paste(r.nd, 2, sep="") 
        }
      } else {
        var.split <- as.character(var.split)
        cut1 <- unlist(strsplit(as.character(cut), split=" ")) #####################
        l.nd <- dats$node[in.node & is.element(var.split, cut1)] 
        r.nd <- dats$node[in.node & !is.element(var.split, cut1)]                  
        dats$node[in.node & is.element(var.split, cut1)] <- paste(l.nd, 1, sep="")  
        dats$node[in.node & !is.element(var.split, cut1)] <- paste(r.nd, 2, sep="")
      }                   
    }
  }
  # print(tre)
  out$data <- dats
  out$tree <- tre.new
  out 
}