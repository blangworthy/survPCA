##Function to get cumulative hazard at a vector of times
cumhaz_time <- function(haz,s){
  sapply(s,function(x)haz$cumhazard$cumhaz[haz$cumhazard$time ==max(haz$cumhazard$time[haz$cumhazard$time <=x])])
}
