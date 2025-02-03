#pred is a posterior matrix (iter x hgrid)
#mod_vec is a vector (with the original modifier categorical values), not the contrast matrix
center_ERF <- function(pred, mod_vec){
  
  groups <- unique(mod_vec)
  pred_centered <- c()
  
  for(group in groups){
    group_idx <- which(mod_vec == group)
    pred_group <- pred[,group_idx]
    group_mean <- rowMeans(pred_group)
    group_center <- pred_group - group_mean
    pred_centered <- cbind(pred_centered, group_center)
  }
  
  return(pred_centered)
}