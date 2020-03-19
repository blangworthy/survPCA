#'Output the Nelson Aalen hazard and cumulative hazard
#'
#'This takes a survival time and censoring indicator and calculates the hazard and cumulative hazard using the Nelson Aalen Estimator
#'@param data An n by 2 matrix where the first column is the observed event time and the second column is the censoring indicator (0 indicates censored)
#'@return A list with the following elements
#'\itemize{
#'\item{hazard: The Nelson-Aalen estimate for the hazard}
#'\item{cumhazard: The Nelson-Aalen estimate for the cumulative hazard}
#'}
#'@export
NelsonAalen <- function(data){
  data <- as.matrix(data)
  if(sum(is.na(data)) > 0){stop("Data includes missing values")}
  uniquetime <- unique(sort(data[,1]))
  n <- dim(data)[1]  
  npts <- length(uniquetime)
  tmp1 <- matrix(rep(data[,1],each = npts),npts,n)
  tmp2 <- matrix(rep(data[,2],each = npts),npts,n)
  atrisk <- ifelse(tmp1>=uniquetime,1,0)
  failure <- ifelse(tmp1==uniquetime & tmp2!=0,1,0)
  
  numrisk <- apply(atrisk,1,sum)
  numfailure <- apply(failure,1,sum)
  if(min(uniquetime) >0){
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
                 "cumhazard" = cumhazard)
  return(output)
}
