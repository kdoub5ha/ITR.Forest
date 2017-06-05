#' @title Finds all decendants of a particular node in a tree.   
#' 
#' @description This function identifies all nodes in a given tree structure
#' which are descendants of the input node. 
#' 
#' @param x node 
#' @param tree interaction tree 
#' @return Returns all the descendant nodes from the node `x`
#' @examples
#' set.seed(10) 
#' dat<- gdataM(1000,2,3,1)
#' tre <- grow.ITR(dat, split.var = 1:4)
#' de('01', tree = tre)
#' "011"  "012"  "0121" "0122"

#' 
#' @export


de <- function(x, tree)
{
  if(length(x) != 1) stop("The length of x in function de must be 1.")    
  y <- tree$node
  de <- NA
  if(sum(match(x, y), na.rm = T) != 0) {
    temp <- 1:length(y)
    start <- match(x, y) + 1    
    end <- length(y)
    if(start <= length(y) & nchar(y[start]) > nchar(x)) {
      temp1 <- temp[temp >= start & nchar(y) <= nchar(x)][1] - 1
      if(!is.na(temp1)) end <- temp1
      de <- y[start:end]
    }}
  de
}
