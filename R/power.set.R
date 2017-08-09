#' @title Categorical predictor subsetting in the tree growing function. 
#' 
#' @param x a categorical predictor. 
#' @return List of unique permutations of predictor `x` 
#' @export

power.set <- function(x) {
  if(length(x) == 0)
    return(vector(mode(x), 0))
  x <- sort(unique(x))
  n <- length(x)
  K <- NULL
  for(m in x) K <- rbind(cbind(K, F), cbind(K, T))
  out <- apply(K, 1, function(x, s) s[x], s = x)
  out <- out[-c(1, length(out))]
  l <- length(out)
  i <- 1
  while (i <= l) {
    if (length(out[[i]]) >= ceiling(n/2+.5)) {
      out[[i]] <- NULL
      i <- i-1
    }
    i <- i+1 
    l <- length(out) 
  }
  out
}  