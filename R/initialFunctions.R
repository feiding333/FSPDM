# write the function get_index_interval to divid the covariate into small bins.
get_index_interval = function(y,num_bin){
  cutinterval_inter = cut(y,breaks = num_bin,labels = FALSE)
  group_interval = list()
  # get the index of each bin
  for(i in 1:num_bin){
    index_interval = which(cutinterval_inter == i )
    group_interval = c(group_interval,list(index_interval))
  }
  return(group_interval)
}
