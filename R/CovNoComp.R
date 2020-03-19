#'Output the martingale and counting process covariance for two variables without presence of competing risks
#'
#'This takes a two survival time and  two censoring indicators and calculates the martingale and counting process covariance at time t
#'@param data An n by 4 matrix where the first column is the observed event time of the first variable, the second column is the observed event time of the second variable, the third column is the censoring indicator for the first variable (0 indicates censored) and the fourth column is the censoring indicator for the second variable (0 indicates censored)
#'@param t The time to evaluate the martingale and counting process covariance at
#'@return A list with the following elements
#'\itemize{
#'\item{MartCov: The martingale covariance between variable 1 and variable 2 at time t}
#'\item{CountCov: The counting process covariance between variable 1 and variable 2 at timet}
#'\item{SurVar1: The Kaplan-Meier survival function estimate for variable 1}
#'\item{Survar2: The Kaplan-Meier survival funciton estimate for variable 2}
#'}
#'@importFrom cmprsk cuminc
#'@importFrom dplyr group_by
#'@importFrom dplyr filter
#'@importFrom survival Surv
#'@importFrom survival survfit
#'
#'@export
CovNoComp <- CovNoComp <- function(data,t){
  if(min(data[,1]<=0)){stop("Observed times should be positive")}
  if(sum(is.na(data)) > 0){stop("Data includes missing values")}
  if(is.numeric(data) == FALSE){stop("Event times and censoring/cause indicator should be numeric")}
  
  data <- as.matrix(data)
  BivSurv <- BivDabrowska(data)
  
  survobj1 <- survival::Surv(time = data[,1],event=ifelse(data[,3]>0,1,0))
  km1 <- survival::survfit(survobj1~1)
  survest1 <- stepfun(km1$time, c(1, km1$surv))
  survobj2 <- survival::Surv(time = data[,2],event=ifelse(data[,4]>0,1,0))
  km2 <- survival::survfit(survobj2~1)
  survest2 <- stepfun(km2$time, c(1, km2$surv))
  tmp <- cbind(survest1(unique(sort(data[,1]))),BivSurv) 
  BivSurvFull <- rbind(c(1,survest2(unique(sort(data[,2])))),tmp)
  # to handle non-monotonicity of Dabrowska's estimator, use cummin
  BivSurvFull <- apply(BivSurvFull,2,cummin)
  BivSurvFull <- t(apply(BivSurvFull,1,cummin))
  rownames(BivSurvFull)[1] <- 0
  colnames(BivSurvFull)[1] <- 0
  
  dhazard1 <- NelsonAalen(data[,c(1,3)])$hazard
  dhazard1 <- dhazard1[dhazard1$time<=t & dhazard1$time > 0,]
  dhazard2 <- NelsonAalen(data[,c(2,4)])$hazard
  dhazard2 <- dhazard2[dhazard2$time<=t & dhazard2$time > 0,]
  prodhaz <- dhazard1$hazard%*%t(dhazard2$hazard)
  rownames(prodhaz) <- dhazard1$time
  colnames(prodhaz) <- dhazard2$time
  
  maxrowbiv <- max(as.numeric(rownames(BivSurvFull))[as.numeric(rownames(BivSurvFull))<=t])
  maxcolbiv <- max(as.numeric(colnames(BivSurvFull))[as.numeric(colnames(BivSurvFull))<=t])
  
  BivSurvFullMinusRow <- BivSurvFull[1:(nrow(BivSurvFull)-1),as.character(maxcolbiv)]
  names(BivSurvFullMinusRow) <- rownames(BivSurvFull)[-1]
  BivSurvFullMinusCol <- BivSurvFull[as.character(maxrowbiv),1:(ncol(BivSurvFull)-1)]
  names(BivSurvFullMinusCol) <- colnames(BivSurvFull)[-1]
  
  BivSurvMinus <- BivSurvFull[1:(nrow(BivSurvFull)-1),1:(ncol(BivSurvFull)-1)]
  colnames(BivSurvMinus) <- colnames(BivSurvFull)[-1] 
  rownames(BivSurvMinus) <- rownames(BivSurvFull)[-1] 
  
  cm <- BivSurvFull[as.character(maxrowbiv),as.character(maxcolbiv)] - 1 +
  sum(dhazard1$hazard * BivSurvFullMinusRow[names(BivSurvFullMinusRow) %in% dhazard1$time]) + 
  sum(dhazard2$hazard * BivSurvFullMinusCol[names(BivSurvFullMinusCol) %in% dhazard2$time]) + 
  sum(prodhaz*BivSurvMinus[rownames(BivSurvMinus)%in%rownames(prodhaz),colnames(BivSurvMinus)%in%colnames(prodhaz)])
  
  cn <- BivSurvFull[as.character(maxrowbiv),as.character(maxcolbiv)] - survest1(t)*survest2(t)
  
  output <- list("MartCov" = cm,"CountCov" = cn,"SurvVar1"= survest1,"SurvVar2" = survest2)
}
