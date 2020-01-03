#'Output the cause specific hazard and cumulative hazard
#'
#'This takes a survival time and censoring/cause indicator and calculates the hazard and cumulative hazard 
#'@param data An n by 2 matrix where the first column is the observed event time and the second column is the censoring/cause indicator (0 indicates censored)
#'@param cause An indicator of which cause to calculate the cause specific hazard for, should be a non-zero value that appears in the cencoring/cause indicator column
#'@return A list with the following elements
#'\itemize{
#'\item{hazard: The cause specific hazard estimate}
#'\item{cumhazard: The cause specific cumulative hazard estimate}
#'\item{cause: The cause value}
#'}
#'@export
CauseHaz<- function(data,cause){
  if(sum(is.na(data)) > 0){stop("Data includes missing values")}
  if(cause ==0){stop("Zero should be for censoring not a cause")}
  if(sum(data[,2] == cause) == 0){warning("Cause not found in censoring/cause indicator")}
  uniquetime <- unique(sort(data[,1][which(data[,2] == cause)]))
  n <- length(data[,1])
  npts <- length(uniquetime)
  tmp1 <- matrix(rep(data[,1],each = npts),npts,n)
  tmp2 <- matrix(rep(data[,2],each = npts),npts,n)
  atrisk <- ifelse(tmp1>=uniquetime,1,0)
  failure <- ifelse(tmp1==uniquetime & tmp2==cause,1,0)
  
  numrisk <- apply(atrisk,1,sum)
  numfailure <- apply(failure,1,sum)
  
  if(min(uniquetime) > 0){
  dhazard <- data.frame(cbind(c(0,uniquetime),c(0,as.matrix(numfailure/numrisk))))
  cumhazard <- data.frame(cbind(c(0,uniquetime),as.matrix(cumsum(dhazard[,2]))))
  }
  else{
    dhazard <- data.frame(cbind(uniquetime,as.matrix(numfailure/numrisk)))
    cumhazard <- data.frame(cbind(uniquetime,as.matrix(cumsum(dhazard[,2]))))
  }
  colnames(dhazard) <- c("time","hazard")
  colnames(cumhazard) <- c("time","cumhazard")
  output <- list("hazard" = dhazard,
                 "cumhazard" = cumhazard,
                 "cause"= cause)
  return(output)
}
