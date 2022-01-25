vertex_mat_make = function(interval_min,interval_max){
  library("tidyverse")

  n_obs=nrow(interval_min)
  
  single_valued_mat =  rbind(interval_min,interval_max) %>%  as_tibble() %>% mutate(id=rep(1:n_obs,2)) %>%
    pivot_longer(cols=1:(ncol(interval_min)-1),names_to = "var_names",values_to="value")
  
  combo_mat = single_valued_mat %>% expand(what,id,var_names)
  
  
  
  min_max_mat = combo_mat %>% dplyr::select(-id) %>%  pivot_wider(values_from = what,names_from = var_names,values_fn=list) %>%
    unnest(cols=everything()) %>% distinct() %>% expand.grid() %>% 
    slice(rep(1:n(),n_obs)) %>% 
    mutate(id=rep(1:n_obs,each=2^(ncol(interval_min)-1)))
  
  
  single_min_max = min_max_mat %>% pivot_longer(cols = 1:(ncol(interval_min)-1),names_to="var_names",values_to="what")%>% 
    left_join(single_valued_mat,by=c("what","var_names","id")) %>% 
    dplyr::select(var_names,value) %>% 
    pivot_wider(names_from=var_names,values_from = value,values_fn=list) %>%
    unnest(cols=everything())
  
  out=list()
  out$vx_mat = single_min_max
  out$ids = rep(1:n_obs,each=2^(ncol(interval_min)-1))
  return(out) 
  
}
