#' Performs k-fold cross validation for selection optimal tuning parameter lambda. 
#' 
#' @param tre original large tree from grow.ITR() 
#' @param dat data frame used to grow tre
#' @param N0 sets the minimal number of observations allowed in terminal nodes. Defaults to 20. 
#' @param n0 sets the minimum number of observations from each treatment group to be in each child node.  Defaults to 5. 
#' @param sp.var specifies the columns of splitting variables in the input data to grow the original tree. 
#' @param nfolds number of folds in the cross validation. Defaults to 10.
#' @param param vector of numeric values to be considered as the tuning parameter. Defaults to seq(0, 0.15, 0.1) but should be modified for each specific problem. 
#' @param AIPWE logical indicating if the augmented estimator should be used. Defaults to FALSE.
#' @param sort logical indicating if data should be sorted before folds are created. Defaults to FALSE. 
#' @param stabilize logical. Should stabilization be used?
#' @param stabilize.type Type of stabilization. 'rf' for random forests (default) or 'linear' for linear model.  
#' @param ctgs columns of categorical variables. 
#' @return A summary of the cross validation including optimal penalty parameter and the optimal model. 
#' @return \item{best.tree}{optimal ITR tree model}
#' @return \item{best.lambda}{optimal lambda}
#' @return \item{full.tree}{original input tree}
#' @return \item{pruned.tree}{pruning output given the optimal lambda}
#' @return \item{data}{input data}
#' @return \item{details}{summary of model performance for each lambda under consideration}
#' @export
#' @examples
#' 
#' # Grow large tree
#' set.seed(1)
#' dat <- rdat(1000, depth = 2)
#' tre <- grow.ITR(dat, split.var = 1:4, min.ndsz = 5, n0 = 2)
#' 
#' # This tre should have 3 terminal nodes (correct size), but has 4 as grown.
#' 
#' cv.prune <- treeCV(tre, dat, nfolds = 5, param = seq(0, 0.15, 0.01), sp.var = 1:4)
#' 
#' # The best tree returned has the correct number of nodes
#' # cv.prune$best.tree
#' #  node size n.1 n.0 trt.effect var vname cut.1 cut.2  score
#' # 1    0 1000 480 520  0.4264782   1    x1     r   0.3 4.8839
#' # 3   01  317 162 155 -1.8782577  NA  <NA>  <NA>  <NA>     NA
#' # 2   02  683 318 365  1.4908419   3    x3     r   0.1 5.0174
#' # 4  021   85  36  49 -1.7575391  NA  <NA>  <NA>  <NA>     NA
#' # 5  022  598 282 316  1.9529338  NA  <NA>  <NA>  <NA>     NA
#' 
#' 


treeCV <- function(tre, dat, nfolds = 5, param = seq(0, 0.15, 0.01), 
                   AIPWE = FALSE, N0=20, n0=5, sp.var, sort = FALSE, ctgs = NA, 
                   stabilize.type = 'rf', stabilize = TRUE){
  input.tre <- tre
  input.dat <- dat
  
  if(stabilize){
    if(stabilize.type == 'rf'){
      fit <- randomForest(y = input.dat$y, as.data.frame(input.dat[,sp.var]))
      resids <- fit$y - fit$predicted
      input.dat$y <- resids
    }
    if(stabilize.type == "linear"){
      fit <- lm(dat$y~as.matrix(dat[ , c(split.var, which(colnames(dat) == "trt")) ]))
      dat$y <- fit$residuals
    }
  }  
  # Shuffle data
  if(sort) input.dat <- input.dat[sample(1:nrow(input.dat), size = nrow(input.dat), replace = FALSE),]
  folds <- cut(seq(1,nrow(input.dat)), breaks = nfolds, labels = FALSE)
  
  in.train <- in.test <- trees <- list()
  result <- NULL
  
  for(k in 1:nfolds){
    in.train[[k]] <- input.dat[-which(folds==k,arr.ind=TRUE),]
    in.test[[k]]  <- input.dat[which(folds==k,arr.ind=TRUE),]
    trees[[k]] <- grow.ITR(in.train[[k]], in.test[[k]], split.var = sp.var, 
                           min.ndsz = N0, n0=n0, AIPWE = AIPWE, ctg = ctgs, 
                           in.forest = TRUE, stabilize.type = stabilize.type)
  }  
  
  out <- matrix(0, ncol = length(trees), nrow = length(param))
  for(t in 1:length(trees)){
    
    pru <- prune(trees[[t]], 0, in.train[[t]], in.test[[t]], AIPWE = FALSE, ctgs = ctgs)
    
    for(j in 1:length(param)){
      row <- which.max(as.numeric(pru$V) - (param[j])*(as.numeric(pru$size.tree)-as.numeric(pru$size.tmnl)))
      out[j, t] <- as.numeric(pru$V.test)[row]
    }
  }
  result <- cbind(apply(out, 1, mean), apply(out, 1, sd))
  
  result2 <- data.frame(Parameter=param, m=result[,1], SD=result[,2], lower=result[,1]-result[,2]/sqrt(length(trees)), 
                        upper=result[,1]+result[,2]/sqrt(length(trees)))
  
  best.lam <- result2$Parameter[which(as.numeric(as.vector(result2$m))==max(as.numeric(as.vector(result2$m))))][1]
  pruned <- prune(input.tre, best.lam, input.dat, AIPWE = FALSE, ctgs = ctgs)
  
  row.prune <- as.numeric(pruned$subtree[as.numeric(pruned$V.a)==max(as.numeric(as.vector(pruned$V.a)))])[1]
  
  best.tree <- input.tre
  if(row.prune > 1){
    for(i in 1:(row.prune-1)){
      best.tree <- best.tree[!is.element(best.tree$node, de(pruned$node.rm[i], best.tree)),]
      if(dim(best.tree)[2]==12) best.tree[which(best.tree$node==pruned$node.rm[i]), 6:11] <- NA
      if(dim(best.tree)[2]==10) best.tree[which(best.tree$node==pruned$node.rm[i]), 6:10] <- NA
    }
  }
  out <- list()
  out$best.tree <- best.tree
  out$best.lambda <- best.lam
  out$full.tree <- input.tre
  out$pruned.tree <- pruned
  out$data <- input.dat
  out$details <- result2
  return(out)
}