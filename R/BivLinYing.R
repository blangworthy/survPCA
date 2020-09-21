#'Output the Lin Ying bivariate survival function estimate
#'
#'This takes survival times and censoring indicators for two variables and outputs the Dabrowska bivariate survival function estimate. Modified from code originally written by Yu Cheng at University of Pittsburgh.
#'@param data An n by 4 matrix where the first column is the observed event time of the first variable, the second column is the observed event time of the second variable, the third column is the censoring indicator for the first variable (0 indicates censored) and the fourth column is the censoring indicator for the second variable (0 indicates censored)
#'@return The Lin Ying bivariate survival function estimate with times for variable 1 along the rows and times for variable 2 along the columns
#'@export
BivLinYing <- function(data){
  data <- as.matrix(data)
  if(sum(is.na(data)) >0){stop("Data includes missing values")}
  ndim = dim(data)[1]
  uniquetime1 = unique(sort(data[,1]))
  uniquetime2 = unique(sort(data[,2]))
  
  npts1 = length(uniquetime1)
  npts2 = length(uniquetime2)
  
  
  tmp1 = matrix(rep(data[,1],each = npts1),npts1,ndim)
  tmp2 = matrix(rep(data[,2],each = npts2),npts2,ndim)
  
  esurv <-  ifelse(tmp1>uniquetime1,1,0)%*%t(ifelse(tmp2>uniquetime2,1,0))/ndim
  
  censortime <- apply(data,1,function(x) max(x[1],x[2]))
  censorind <- apply(data,1,function(x) x[3]*x[4])
  
  censkm <- survival::survfit(survival::Surv(censortime,(1-censorind))~1)
  censsurv <- stepfun(censkm$time, c(1, censkm$surv))
  
  surv1 <- censsurv(uniquetime1)
  surv2 <- censsurv(uniquetime2)
  censsurvmat <- sapply(surv2,function(x) sapply(surv1,function(y) min(x,y)))
  
  mat <- esurv/censsurvmat
  
  rownames(mat) <- uniquetime1
  colnames(mat) <- uniquetime2
  return(mat)
  
}
  