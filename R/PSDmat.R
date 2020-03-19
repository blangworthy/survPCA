###Make a matrix positive semidfinite if it is not already
PSDmat <- function(mat,mineigen){
if(min(eigen(mat)$values)>0){
  return(mat)
}else{
  return(eigen(mat)$vectors%*%diag(pmax(mineigen,eigen(mat)$values))%*%t(eigen(mat)$vectors))
}
}
