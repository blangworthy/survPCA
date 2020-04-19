#'Output the entire martingale and counting process covariance and correlation matrices for p different variables at time t.
#'
#'This takes survival times and censoring/causue indicators for p different variables and caclulates the p times p martingale and counting process covariance matrices for all variables.
#'@param data A matrix where the first p columns are the event times for the p variables and the next p columns are the censoring indicators for the p variables. The order for the event times and censoring indicators needs to be the same.
#'@param p The number of variables to calculate the matrix for, if NA ncol(data)/2 will be used
#'@param t The timepoint to estimate the covariance and correlation matrices
#'@param mineigen The minimum eigenvalue for the martingale and counting process covariance matrices if not already positive semidefinite
#'@return A list with the following elements
#'\itemize{
#'\item{CovMart: The full martingale covariance matrix}
#'\item{CorMart: The full martingale correlation matrix}
#'\item{CovCount: The full counting process covariance matrix}
#'\item{CorCount: The full counting process correlation matrix}
#'}
#'@export
CovMatNoComp <- function(data,p=NA,t,mineigen = 0.001){
  if(is.na(p)){p=ncol(data)/2}
  data <- as.matrix(data[,1:(2*p)])
  if(sum(is.na(data)) > 0){stop("Data includes missing values")}
  if(is.numeric(data) == FALSE){stop("Event times and censoring/cause indicator should be numeric")}
  
  CovMart <- matrix(nrow = p,ncol = p)
  CovCount <- matrix(nrow = p,ncol = p)
  
  for(x in 1:p){
    for(y in x:p){
      if(x == y){
        survobj1 <- survival::Surv(time = data[,x],event=ifelse(data[,(x+p)]>0,1,0))
        km1 <- survival::survfit(survobj1~1)
        survest1 <- stepfun(km1$time, c(1, km1$surv))
        CovMart[x,x] <- 1-survest1(t)
        CovCount[x,x] <- survest1(t)*(1-survest1(t))
      }
      else{
        covout <- CovNoComp(data[,c(x,y,(p+x),(p+y))],t)
        CovMart[x,y] <- covout$MartCov
        CovCount[x,y] <- covout$CountCov
      }
    }
  }
  
  
  CovMart[lower.tri(CovMart)] <- t(CovMart)[lower.tri(CovMart)]
  CovMart <- PSDmat(CovMart,mineigen)
  CorMart <- cov2cor(CovMart)
  
  CovCount[lower.tri(CovCount)] <- t(CovCount)[lower.tri(CovCount)]
  CovCount <- PSDmat(CovCount,mineigen)
  if(0 %in% diag(CovCount)){
    CorCount <- NA
    warning("Variance of one or more of counting processes estimated to be 0, correlation matrix not well defined")
  }else{
  CorCount <- cov2cor(CovCount)
  }
  output <- list("CovMart" = CovMart,
                 "CorMart" = CorMart,
                 "CovCount" = CovCount,
                 "CorCount" = CorCount)
  
  return(output)
}