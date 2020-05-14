#'Output the martingale and counting process covariance for two variables in presence of competing risks
#'
#'This takes survival times and censoring/cause indicators for two variables and outputs the martingale and counting process covariance at time t
#'@param data An n by 4 matrix where the first column is the observed event time of the first variable, the second column is the observed event time of the second variable, the third column is the censoring/cause indicator for the first variable (0 indicates censored) and the fourth column is the censoring/cause indicator for the second variable (0 indicates censored)
#'@param t The time to evaluate the martingale and counting process covariance at
#'@param cause1 An indicator of which cause to calculate the cause specific hazard for for the first variable, should be a non-zero value that appears in the censoring/cause indicator column for the first variable
#'@param cause2 An indicator of which cause to calculate the cause specific hazard for for the second variable, should be a non-zero value that appears in the censoring/cause indicator column for the second variable
#'@return A list with the following elements
#'\itemize{
#'\item{MartCov: The martingale covariance between variable 1 and variable 2 at time t}
#'\item{CountCov: The counting process covariance between variable 1 and variable 2 at timet}
#'\item{CIFVar1: The cause specifice cumulative incidence function estimate for variable 1}
#'\item{Survar2: The cause specifice cumulative incidence function estimate for variable 2}
#'}
#'@importFrom cmprsk cuminc
#'@importFrom dplyr group_by
#'@importFrom dplyr filter
#'@importFrom survival Surv
#'@importFrom survival survfit
#'
#'@export
CovComp <- function(data,t,cause1,cause2){
  
  data <- as.matrix(data)
  if(min(data[,1]<=0)){stop("Observed times should be positive")}
  if(sum(is.na(data)) > 0){stop("Data includes missing values")}
  if(is.numeric(data) == FALSE){stop("Event times and censoring/cause indicator should be numeric")}
  if(cause1 ==0){stop("Zero should be for censoring not a cause")}
  if(cause2 ==0){stop("Zero should be for censoring not a cause")}
  if(sum(data[,3] == cause1) == 0){warning("Cause for variable 1 not found in censoring/cause indicator")}
  if(sum(data[,4] == cause2) == 0){warning("Cause for variable 2 not found in censoring/cause indicator")}
  haz1 <- CauseHaz(data[,c(1,3)],cause1)
  haz2 <- CauseHaz(data[,c(2,4)],cause2)
  
  cumincout1 <- cmprsk::cuminc(data[,1],data[,3])
  cumincout2 <- cmprsk::cuminc(data[,2],data[,4])
  
  cif1 <- data.frame(cbind(cumincout1[[cause1]]$time,cumincout1[[cause1]]$est))
  cif1 <- dplyr::filter(dplyr::group_by(cif1,X1),X2 == max(X2))
  cif2 <- data.frame(cbind(cumincout2[[cause2]]$time,cumincout2[[cause2]]$est))
  cif2 <- dplyr::filter(dplyr::group_by(cif2,X1),X2 == max(X2))
 
  
  cif1fun <- stepfun(cif1$X1[-1],cif1$X2)
  cif2fun <- stepfun(cif2$X1[-1],cif2$X2)
  
  if1 <- diff(cif1$X2)
  if2 <- diff(cif2$X2)
  
  names(if1) <- cif1$X1[-1]
  names(if2) <- cif2$X1[-1]
  
  
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
  
  BivHaz11 <- BivCauseHaz(data,1,1)
  BivHaz12 <- BivCauseHaz(data,1,2)
  BivHaz21 <- BivCauseHaz(data,2,1)
  BivHaz22 <- BivCauseHaz(data,2,2)
  
  BivCIF11 <- BivCIF(BivHaz11$hazard,BivSurvFull)
  BivCIF12 <- BivCIF(BivHaz12$hazard,BivSurvFull)
  BivCIF21 <- BivCIF(BivHaz21$hazard,BivSurvFull)
  BivCIF22 <- BivCIF(BivHaz22$hazard,BivSurvFull)
  
  SurvTimesRow <-as.numeric(rownames(BivSurvFull))[as.numeric(rownames(BivSurvFull))<=t]
  SurvTimesCol <-as.numeric(colnames(BivSurvFull))[as.numeric(colnames(BivSurvFull))<=t]
  Haz11TimesRow <-as.numeric(rownames(BivHaz11$hazard))[as.numeric(rownames(BivHaz11$hazard))<=t]
  Haz11TimesCol <-as.numeric(colnames(BivHaz11$hazard))[as.numeric(colnames(BivHaz11$hazard))<=t]
  Haz12TimesRow <-as.numeric(rownames(BivHaz12$hazard))[as.numeric(rownames(BivHaz12$hazard))<=t]
  Haz12TimesCol <-as.numeric(colnames(BivHaz12$hazard))[as.numeric(colnames(BivHaz12$hazard))<=t]
  Haz21TimesRow <-as.numeric(rownames(BivHaz21$hazard))[as.numeric(rownames(BivHaz21$hazard))<=t]
  Haz21TimesCol <-as.numeric(colnames(BivHaz21$hazard))[as.numeric(colnames(BivHaz21$hazard))<=t]
  Haz22TimesRow <-as.numeric(rownames(BivHaz22$hazard))[as.numeric(rownames(BivHaz22$hazard))<=t]
  Haz22TimesCol <-as.numeric(colnames(BivHaz22$hazard))[as.numeric(colnames(BivHaz22$hazard))<=t]
  Haz1TimesRow <- haz1$cumhazard$time[haz1$cumhazard$time <= t]
  Haz2TimesRow <- haz2$cumhazard$time[haz2$cumhazard$time <= t]
  IF1Times <- as.numeric(names(if1))[as.numeric(names(if1))<= t & if1 > 0]
  IF2Times <- as.numeric(names(if2))[as.numeric(names(if2))<= t & if2 > 0]
  
  term1 <- haz1$cumhazard$cumhaz[haz1$cumhazard$time ==max(Haz1TimesRow)]*haz2$cumhazard$cumhaz[haz2$cumhazard$time ==max(Haz2TimesRow)]*BivSurvFull[as.character(max(SurvTimesRow)),as.character(max(SurvTimesCol))]
  term2 <- sum(colSums(BivCIF22$BivIF[as.character(Haz22TimesRow),as.character(Haz22TimesCol)])*cumhaz_time(haz1,Haz22TimesCol)*cumhaz_time(haz2,Haz22TimesCol))
  term3 <- sum(t(BivCIF11$BivIF[as.character(Haz11TimesRow),as.character(Haz11TimesCol)]*(1-cumhaz_time(haz1,Haz11TimesRow)))*(1-cumhaz_time(haz2,Haz11TimesCol)))
  term4 <- sum((if1[as.character(IF1Times)] - rowSums( BivCIF11$BivIF[as.character(Haz11TimesRow),as.character(Haz11TimesCol)]) - rowSums( BivCIF12$BivIF[as.character(Haz12TimesRow),as.character(Haz12TimesCol)]))*(1-cumhaz_time(haz1,IF1Times))*(-1*haz2$cumhazard$cumhaz[haz2$cumhazard$time ==max(Haz2TimesRow)]))
  term5 <- sum((if2[as.character(IF2Times)] - colSums( BivCIF11$BivIF[as.character(Haz11TimesRow),as.character(Haz11TimesCol)]) - colSums( BivCIF21$BivIF[as.character(Haz21TimesRow),as.character(Haz21TimesCol)]))*(1-cumhaz_time(haz2,IF2Times))*(-1*haz1$cumhazard$cumhaz[haz1$cumhazard$time ==max(Haz1TimesRow)]))
  term6 <- sum(t(BivCIF12$BivIF[as.character(Haz12TimesRow),as.character(Haz12TimesCol)]*(1-cumhaz_time(haz1,Haz12TimesRow)))*(-1*cumhaz_time(haz2,Haz12TimesCol)))
  term7 <- sum(t(BivCIF21$BivIF[as.character(Haz21TimesRow),as.character(Haz21TimesCol)]*(-1*cumhaz_time(haz1,Haz21TimesRow)))*(1-cumhaz_time(haz2,Haz21TimesCol)))
  
  cm <- sum(term1,term2,term3,term4,term5,term6,term7)
  
  bivcif <- get(paste0("BivCIF",cause1,cause2))
  maxrow <- max(as.numeric(rownames(bivcif$BivCIF))[as.numeric(rownames(bivcif$BivCIF))<=t])
  maxcol <- max(as.numeric(colnames(bivcif$BivCIF))[as.numeric(colnames(bivcif$BivCIF))<=t])
  
  cn <- bivcif$BivCIF[as.character(maxrow),as.character(maxcol)] - cif2fun(t)*cif2fun(t)
  
  output <- list("MartCov" = cm,"CountCov" = cn,"CIFVar1" = cif1fun,"CIFVar2" = cif2fun)
  
  
 return(output)
}



