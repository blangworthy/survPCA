#'Output the Dabrowska bivariate survival function estimate
#'
#'This takes two survival times and censoring indicators and outputs the Dabrowska bivariate survival function estimate
#'@param data An n by 4 matrix where the first column is the observed event time of the first variable, the second column is the observed event time of the second variable, the third column is the censoring indicator for the first variable (0 indicates censored) and the fourth column is the censoring indicator for the second variable (0 indicates censored)
#'@return The Debrowska bivariate survival function estimate with times for variable 1 on the rows and times for variable 2 on the columns
#'@export
BivDabrowska= function(data){
  if(sum(is.na(data)) >0){stop("Data includes missing values")}
  ###########################################  
  ## bivariate event and at-risk processes ##
  ###########################################  
  ndim = dim(data)[1]        
  uniquetime1 = unique(sort(data[,1]))
  uniquetime2 = unique(sort(data[,2]))
  
  npts1 = length(uniquetime1)
  npts2 = length(uniquetime2)
  
  tmp1 = matrix(rep(data[,1],each = npts1),npts1,ndim)
  tmp3 = matrix(rep(data[,3],each = npts1),npts1,ndim)
  tmp2 = matrix(rep(data[,2],each = npts2),npts2,ndim)
  tmp4 = matrix(rep(data[,4],each = npts2),npts2,ndim)
  
  numrisk = ifelse(tmp1>=uniquetime1,1,0)%*%t(ifelse(tmp2>=uniquetime2,1,0))
  lambda11 = ifelse(tmp1==uniquetime1 & tmp3!=0,1,0)%*%t(ifelse(tmp2==uniquetime2 & tmp4 !=0,1,0))/numrisk
  lstar = ifelse(tmp1==uniquetime1 & tmp3 ==1,1,0)%*%t(ifelse(tmp2==uniquetime2 & tmp4 ==1,1,0))/numrisk
  
  # the factor is evaluated at lambda_{10} (x,y-) and lambda_{01}(x-,y)
  lambda10.left =ifelse(tmp1==uniquetime1 & tmp3 !=0,1,0)%*%t(ifelse(tmp2>=uniquetime2,1,0))/numrisk
  lambda01.up = ifelse(tmp1>=uniquetime1,1,0)%*%t(ifelse(tmp2==uniquetime2 & tmp4 !=0,1,0))/numrisk
  
  numrisk = NULL
  
  
  lambda10.left = ifelse(lambda10.left == "NaN", 0, lambda10.left)
  lambda01.up = ifelse(lambda01.up == "NaN", 0, lambda01.up)
  lambda11 = ifelse(lambda11 == "NaN", 0, lambda11)   
  lstar = ifelse(lstar == "NaN", 0, lstar)
  
  
  lambda1 = apply(ifelse(tmp1==uniquetime1 & tmp3 != 0,1,0),1,sum)/apply(ifelse(tmp1>=uniquetime1,1,0),1,sum)   
  S1 = cumprod(1-lambda1)
  lambda2 = apply(ifelse(tmp2==uniquetime2 & tmp4 != 0,1,0),1,sum)/apply(ifelse(tmp2>=uniquetime2,1,0),1,sum)   
  S2 = cumprod(1-lambda2)
  
  lstar2 = apply(ifelse(tmp2==uniquetime2 & tmp4 == 1,1,0),1,sum)/apply(ifelse(tmp2>=uniquetime2,1,0),1,sum)
  
  
  ## Dabrowska estimator for the bivariate survival function
  tmp1 = tmp2 = tmp3 = tmp4 = tmp5 = tmp6  = NULL
  # factor
  factor = ifelse(lambda10.left!=1 & lambda01.up!=1, 1-(lambda10.left*lambda01.up-lambda11)/(1-lambda10.left)/(1-lambda01.up), 1)
  
  
  factor = cumprodmatrix(factor)   
  factor = ifelse(is.na(factor),0,factor)
  
  SS = S1%*%t(S2)*factor
  
  factor = NULL
  
  # to handle non-monotonicity of Dabrowska's estimator, use cummin
  SS = apply(SS,2,cummin)
  SS = t(apply(SS,1,cummin))
  rownames(SS) <- uniquetime1
  colnames(SS) <- uniquetime2
  
  return(SS)
}
