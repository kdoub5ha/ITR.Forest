#' @title Treatment Prediction Function
#'
#' @description Used to make treatment prediction for a single tree or random forest. If the 
#' input is a forest, then the proportion of trees voting for treatment (`trt=1`) is returned. 
#' If the input is a single tree, then the function returns the vote (treatment or control) from the tree. 
#' 
#' @param input tree or forest object from `grow.ITR` or `Build.RF.ITR`.
#' @param new.dat data for which predictions are desired
#' @param ctgs columns of categorical variables. 
#' @return A summary list of the following elements:
#' @return \item{SummaryTreat}{proportion of trees voting for treatment (trt=1). 
#' If input is a single tree then SummaryTreat is a single number. 
#' If input is a forest then SummaryTreat is a vector equal to the length of the number of trees.}
#' @return \item{trt.pred}{vector of treatment assignments {0, 1} based on the tree vote (single tree) or majority of tree votes (forest). This vector has length equal to the number of rows in `new.dat`.}
#' @return \item{n.trees}{number of tree in `input`}
#' @return \item{tree.votes}{matrix of votes for each tree for each subject in `new.dat`. Rows correspond to trees in `input` and columns correspond to subjects in `new.dat`.}
#' @return \item{data}{input data frame `new.dat`}
#' @return \item{NA.trees}{number of trees returning no votes. In a forest, this is the number of null trees.}
#' @export
#' @examples
#' dat <- gdataM(n=1000, depth=2, beta1=3, beta2=1)
#' # Build a forest with 100 trees
#' forest <- Build.RF.ITR(dat, col.y="y", col.trt="trt", col.prtx="prtx", split.var=1:4, ntree=100)
#' # Predict treatment assignments for 1000 observations in `dat`
#' predict.ITR(forest, dat)



predict.ITR <- function(input, new.dat, ctgs = NULL){
  if(is.null(dim(input))) trees <- input$TREES
  if(!is.null(dim(input))) trees <- input
  dat <- new.dat
  n <- nrow(dat)
  if(is.null(dim(input))) n.trees <- length(trees)
  if(!is.null(dim(input))) n.trees <- 1
  out <- NULL
  
  result <- NULL
  for(i in 1:n.trees){
    if(is.null(dim(input))) tre <- trees[[i]]
    if(!is.null(dim(input))) tre <- trees
    
    if(nrow(tre)>1){
      send<-send.down(dat.new = dat, tre = tre, ctgs = ctgs)
      node<-substr(send$data$node,1,nchar(send$data$node)-1)
      direction<-substr(send$data$node,nchar(send$data$node),nchar(send$data$node))
      trt.dir<-tre[match(node,tre$node),]$cut.1
      
      trt.pred<-ifelse(trt.dir=="r" & direction=="1",0,
                       ifelse(trt.dir=="r" & direction=="2",1,
                              ifelse(trt.dir=="l" & direction=="1",1,0)))
    }else{
      if(is.null(dim(input))) trt.pred <- rep(NA, n)
      if(!is.null(dim(input))) trt.pred <- NA
    }
    result <- rbind(result, trt.pred) 
  }
  out$SummaryTreat <- apply(result, 2, FUN = mean, na.rm=T)
  if(is.null(dim(input))) out$trt.pred <- ifelse(out$SummaryTreat < 0.5, 0, 1)
  if(!is.null(dim(input))) out$trt.pred <- out$SummaryTreat
  out$n.trees <- n.trees
  out$tree.votes <- result
  out$data <- new.dat
  out$NA.trees <- sum(is.na(result[,1]))
  return(out)
}