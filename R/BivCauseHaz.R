#'Output the bivariate cause specific cumulative hazard and jumps in the bivariate cause specific cumulative hazard
#'
#'This takes survival times and censoring/cause indicators for two variables and outputs the bivariate cause specific cumulative hazard and jumps in the bivariate cause specific cumulative hazard
#'@param data An n by 4 matrix where the first column is the observed event time of the first variable, the second column is the observed event time of the second variable, the third column is the censoring/cause indicator for the first variable (0 indicates censored) and the fourth column is the censoring/cause indicator for the second variable (0 indicates censored)
#'@param cause1 An indicator of which cause to calculate the cause specific hazard for for the first variable, should be a non-zero value that appears in the censoring/cause indicator column for the first variable
#'@param cause1 An indicator of which cause to calculate the cause specific hazard for for the second variable, should be a non-zero value that appears in the cencoring/cause indicator column for the second variable
#'@return A list with the following elements:
#'\itemize{
#'\item{hazard: Jumps in the bivariate cause specific cumulative hazard}
#'\item{cumhazard: The bivariate cause specific cumulative hazard}
#'\item{cause1: The cause for the first variable}
#'\item{cause2: The cause for the second variable}
#'}
#'
#'@export
BivCauseHaz <- function(data,cause1,cause2){
  data <- as.matrix(data)
  if(sum(is.na(data)) > 0){stop("Data includes missing values")}
  if(cause1 ==0){stop("Zero should be for censoring not a cause")}
  if(cause2 ==0){stop("Zero should be for censoring not a cause")}
  if(sum(data[,3] == cause1) == 0){warning("Cause for variable 1 not found in censoring/cause indicator")}
  if(sum(data[,4] == cause2) == 0){warning("Cause for variable 2 not found in censoring/cause indicator")}
  uniquetime1 <- unique(sort(data[,1][which(data[,3] == cause1)]))
  uniquetime2 <- unique(sort(data[,2][which(data[,4] == cause2)]))
  
  n <- length(data[,1])
  
  npts1 <- length(uniquetime1)
  npts2 <- length(uniquetime2)
  
  tmp1 <- matrix(rep(data[,1],each = npts1),npts1,n)
  tmp2 <- matrix(rep(data[,3],each = npts1),npts1,n)
  tmp3 <- matrix(rep(data[,2],each = npts2),npts2,n)
  tmp4 <- matrix(rep(data[,4],each = npts2),npts2,n)
  
  numrisk <- ifelse(tmp1>=uniquetime1,1,0)%*%t(ifelse(tmp3>=uniquetime2,1,0))
  dhaz <-  ifelse(tmp1==uniquetime1 & tmp2==cause1,1,0)%*%t(ifelse(tmp3==uniquetime2 & tmp4 ==cause2,1,0))/numrisk
  dhaz <- ifelse(is.na(dhaz),0,dhaz)
  cumhaz <- t(apply(apply(dhaz, 2, cumsum), 1, cumsum))
  rownames(cumhaz) <- uniquetime1
  colnames(cumhaz) <- uniquetime2
  rownames(dhaz) <- uniquetime1
  colnames(dhaz) <- uniquetime2
  return(list("hazard" =dhaz,"cumhazard" =cumhaz,
              "cause1" = cause1,"cause2" = cause2))
}
