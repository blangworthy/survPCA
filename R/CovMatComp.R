#'Output the entire martingale and counting process covariance and correlation matrices for p different variables at time t in the presence of competing risks
#'
#'This takes survival times and censoring/causue indicators for p different variables and caclulates the p times p martingale and counting process covariance matrices for all variables in the presence of competing risks.
#'@param data A matrix where the first p columns are the event times for the p variables and the next p columns are the censoring/cause indicators for the p variables. The order for the event times and censoring/cause indicators needs to be the same.
#'@param causes A p-dimensional vector that has the causes for each of the p different variables.
#'@param p The number of variables to calculate the matrix for, if NA ncol(data)/2 will be used
#'@param t The timepoint to estimate the covariance and correlation matrices
#'@param mineigen The minimum eigenvalue for the martingale and counting process covariance matrices if not already positive semidefinite
#'@return A list with the following elements
#'\itemize{
#'\item{CovMartComp: The full martingale covariance matrix}
#'\item{CorMartComp: The full martingale correlation matrix}
#'\item{CovCountComp: The full counting process covariance matrix}
#'\item{CorCountComp: The full counting process correlation matrix}
#'}
#'@export
CovMatComp <- function(data,causes,p=NA,t,mineigen = 0.001){
  if(is.na(p)){p=ncol(data)/2}
  data <- as.matrix(data[,1:(2*p)])
  if(sum(is.na(data)) > 0){stop("Data includes missing values")}
  if(is.numeric(data) == FALSE){stop("Event times and censoring/cause indicator should be numeric")}
  
  CovMartComp <- matrix(nrow = p,ncol = p)
  CovCountComp <- matrix(nrow = p,ncol = p)
  for(x in 1:(p)){
    for(y in (x):p){
      if(x==y){
        cumincout1 <- cmprsk::cuminc(data[,x],data[,(x+p)])
        cif1 <- data.frame(cbind(cumincout1[[causes[p]]]$time,cumincout1[[causes[p]]]$est))
        cif1 <- dplyr::filter(dplyr::group_by(cif1,X1),X2 == max(X2))
        cif1fun <- stepfun(cif1$X1[-1],cif1$X2)
        CovMartComp[x,x] <- cif1fun(t)
        CovCountComp[x,x] <- cif1fun(t)*(1-cif1fun(t))
      }
      else{
        covout <- CovComp(data[,c(x,y,(p+x),(p+y))],t,causes[x],causes[y])
        CovMartComp[x,y] <- covout$MartCov
        CovCountComp[x,y] <- covout$CountCov
      }
    }
  }
  
  CovMartComp[lower.tri(CovMartComp)] <- t(CovMartComp)[lower.tri(CovMartComp)]
  CovMartComp <- PSDmat(CovMartComp,mineigen)
  CorMartComp <- cov2cor(CovMartComp)
  
  CovCountComp[lower.tri(CovCountComp)] <- t(CovCountComp)[lower.tri(CovCountComp)]
  CovCountComp <- PSDmat(CovCountComp,mineigen)
  CorCountComp <- cov2cor(CovCountComp)
  
  output <- list("CovMartComp" = CovMartComp,
                 "CorMartComp" = CorMartComp,
                 "CovCountComp" = CovCountComp,
                 "CorCountComp" = CorCountComp)
  return(output)
  
  
}