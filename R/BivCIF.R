#'Output estimate for the bivariate cumulative incidence function
#'
#'Takes bivariate survival and bivariate hazard function and outputs estimate for bivariate cumulative incidence function
#'@param BivHaz A two dimensional matrix with estimates of bivariate cause specific hazard function. Times for variable 1 are along the rows, times for variable 2 along columns.
#'@param BivSurv A two dimensional matrix with estimates of bivariate survival times. Times for variable 1 are along the rows, times for variable 2 along columns.
#'@return A list with the following elements:
#'\itemize{
#'\item{BivIF: The bivariate cause specific incidence function estimate}
#'\item{BivCIF: The bivariate cause specific cumulative incidence function estimate}
#'}
#'@export
BivCIF<- function(BivHaz,BivSurv){
  BivSurvMinus <- BivSurv[1:(nrow(BivSurv)-1),1:(ncol(BivSurv)-1)]
  colnames(BivSurvMinus) <- colnames(BivSurv)[-1] 
  rownames(BivSurvMinus) <- rownames(BivSurv)[-1] 
  
  BivSurvMinus <- BivSurvMinus[rownames(BivHaz),colnames(BivHaz)]
  BivIF <- BivSurvMinus*BivHaz
  BivCIF <- t(apply(apply(BivIF , 2, cumsum), 1, cumsum))
  rownames(BivCIF) <- rownames(BivHaz)
  colnames(BivCIF) <- colnames(BivHaz)
  rownames(BivIF) <- rownames(BivHaz)
  colnames(BivIF) <- colnames(BivHaz)
  return(list("BivIF" = BivIF,"BivCIF" = BivCIF))
}
