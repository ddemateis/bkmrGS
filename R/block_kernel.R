#mod_vec1 is the first modifier contrast matrix used to construct the covariance matrix
#mod_vec1 is the second modifier contrast matrix used to construct the covariance matrix
#K is the covariance matrix that needs to be blocked by modifier group

block_kernel <- function(mod_vec1, mod_vec2, K){
  
  #convert the contrast modifier matrices to unique modifier groups 
  m1 <- apply(mod_vec1, 1, paste, collapse = "")
  m2 <- apply(mod_vec2, 1, paste, collapse = "")

  #construct a logical matrix to block K based on groups
  block_idx <- matrix(apply(expand.grid(m1,m2), 
                            1, 
                            function(x)x[1]==x[2]), 
                      ncol = length(m2))
  
  #return the original K matrix with elements in different groups set to 0
  return(K * block_idx)
}

