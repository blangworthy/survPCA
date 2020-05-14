#Modified from code received from Yu Cheng at University of Pittsburgh.
# subfunctions to calculate cumulative sum function of a matrix
cummatrix = function(ee)
{n2 = dim(ee)[1]
n4 = dim(ee)[2]

for (i in 1:n4)
  ee[,i] = cumsum(ee[,i])
for (i in 1:n2)
  ee[i,] = cumsum(ee[i,])
return(ee)
}

# subfunctions to calculate cumulative prod function of a matrix
cumprodmatrix = function(ee)
{n2 = dim(ee)[1]
n4 = dim(ee)[2]

for (i in 1:n4)
  ee[,i] = cumprod(ee[,i])

for (i in 1:n2)
  ee[i,] = cumprod(ee[i,])
return(ee)
}


