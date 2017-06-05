#' Treatment Prediction Function
#'
#' Used to make treatment prediction for a single tree or random forest
#' @param input tree or forest object from `grow.ITR` or `Build.RF.ITR`.
#' @param new.dat data for which predictions are desired
#' @export
#' @examples
#' dat <- gdataM(n=1000, depth=2, beta1=3, beta2=1)
#' #Build a forest with 100 trees
#' forest <- Build.RF.ITR(dat, col.y="y", col.trt="trt", col.prtx="prtx", split.var=1:4, ntree=100)
#' #Predict treatment assignments for 1000 observations in `dat`
#' predict.ITR(forest, dat)



predict.ITR <- function(input, new.dat){
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
      send<-send.down(dat.new = dat, tre = tre)
      node<-substr(send$data$node,1,nchar(send$data$node)-1)
      direction<-substr(send$data$node,nchar(send$data$node),nchar(send$data$node))
      trt.dir<-tre[match(node,tre$node),]$cut.1
      
      trt.pred<-ifelse(trt.dir=="r" & direction=="1",0,
                       ifelse(trt.dir=="r" & direction=="2",1,
                              ifelse(trt.dir=="l" & direction=="1",1,0)))
    }else{
      trt.pred <- NA
    }
    result <- rbind(result, trt.pred) 
  }
  out$SummaryTreat <- apply(result, 2, FUN = mean, na.rm=T)
  out$n.trees <- n.trees
  out$tree.votes <- result
  out$data <- new.dat
  out$NA.trees <- sum(is.na(result[,1]))
  return(out)
}