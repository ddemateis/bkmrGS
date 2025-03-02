#mod_vec1 is the first modifier contrast matrix used to construct the covariance matrix
#mod_vec1 is the second modifier contrast matrix used to construct the covariance matrix
#K is the covariance matrix that needs to be blocked by modifier group

block_kernel <- function(mod_vec1, mod_vec2, K){
  
  #convert the contrast modifier matrices to unique modifier groups 
  m1 <- apply(mod_vec1, 1, paste, collapse = "")
  m2 <- apply(mod_vec2, 1, paste, collapse = "")

  #construct a logical matrix to block K based on groups
  block_idx <- t(sapply(m1, function(x) x == m2))
  
  #return the original K matrix with elements in different groups set to 0
  return(K * block_idx)
}

#mod_vec1 is the first modifier contrast matrix used to construct the covariance matrix
#mod_vec1 is the second modifier contrast matrix used to construct the covariance matrix
#lambda is a row of the posterior samples, vector of each group's lambda
#data.comps is a list of info, need gs.tau and randint
lambda_scalar <- function(mod1, mod2, lambda, data.comps){
  
  if(data.comps$gs.tau){#return scalar matrix of lambda if group-specific tau method
    
    #don't use the last element of lambda if data.comps$randint is true
    if(data.comps$randint){
      lambda <- lambda[1:(length(lambda)-1)]
    }
    
    #get unique levels from modifier contrast matrices
    m1 <- apply(mod1, 1, paste, collapse = "")
    
    #match the lambdas to the modifier levels
    lambda_col <- lambda[m1] 
    
    #take a row/col and repeat it into a matrix
    #this will be multiplied (element-wise) by a blocked K matrix, so it doesn't matter that this doesn't have the block structure
    scalar <- matrix(rep(lambda_col, nrow(mod2)), ncol=nrow(mod2))
    
  }else{#return lambda as a scalar
    scalar <- lambda[1]
  }
  
  return(scalar)
}