#'Output the Lin Ying bivariate survival function estimate
#'
#'This takes survival times and censoring indicators for two variables and outputs the Dabrowska bivariate survival function estimate. Modified from code originally written by Yu Cheng at University of Pittsburgh.
#'@param data An n by 4 matrix where the first column is the observed event time of the first variable, the second column is the observed event time of the second variable, the third column is the censoring indicator for the first variable (0 indicates censored) and the fourth column is the censoring indicator for the second variable (0 indicates censored)
#'@return The Lin Ying bivariate survival function estimate with times for variable 1 along the rows and times for variable 2 along the columns
#'@export
BivLinYing <- function(data){
  data <- as.matrix(data)
  if(sum(is.na(data)) >0){stop("Data includes missing values")}
  
  uniquetime1 = unique(sort(data[,1]))
  uniquetime2 = unique(sort(data[,2]))
  
  censortime <- apply(data,1,function(x) max(x[1],x[2]))
  censorind <- apply(data,1,function(x) x[3]*x[4])
  
  censkm <- survival::survfit(survival::Surv(censortime,(1-censorind))~1)
  censsurv <- stepfun(censkm$time, c(1, censkm$surv))

  
  
  mat <- sapply(uniquetime2,function(x) sapply(uniquetime1, function(y) ifelse(censsurv(max(x,y))>0,mean(data[,2]>x & data[,1] > y)/censsurv(max(x,y)),NA)))
  rownames(mat) <- uniquetime1
  colnames(mat) <- uniquetime2
  return(mat)
  
}
  