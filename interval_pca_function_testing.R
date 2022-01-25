rm(list=ls())
library(tidyverse)
source("interval_PCA.r")
source("plot.interval_PCA_cpca.R")
source("plot.interval_PCA_vpca.R")
 F_Rec = read_delim(file="./Dataset/Face_recognition.csv",delim=";") 
# F_Rec %>% filter(str_detect(`...1`,pattern="ISA"))
 F_Rec = Face_recognition

obs_names=F_Rec %>% dplyr::select(`...1`)
F_Rec=F_Rec %>% dplyr::select(-`...1`)
var_names = F_Rec %>% names %>% .[seq(from=2,to=ncol(F_Rec),by=2)]
# F_Rec=read_csv2(file="./Dataset/peperoni4vars.csv") 
# obs_names=F_Rec %>% select(Id)
# F_Rec=F_Rec %>% select(-Id)

out_cpca=interval_PCA(x=F_Rec,method = "cpca",std_meth="centers",obs_names = obs_names %>% pull,var_names = var_names)
cpca_plot=plot(out_cpca)
cpca_plot$rect_plot
ggsave("C_PCA_f_r_sol.pdf", cpca_plot$rect_plot, width=8, height=6)

out_vpca=interval_PCA(x=F_Rec,method = "vpca",std_meth="centers",obs_names = obs_names %>% pull,var_names = var_names)
vpca_plot=plot(out_vpca)
vpca_plot$vpca_plot
ggsave("V_PCA_f_r_sol.pdf", vpca_plot$vpca_plot, width=8, height=6)

out_mrpca=interval_PCA(x=F_Rec,method = "mrpca",std_meth = "symb",obs_names = obs_names)


out_mrpca$observations_names=(obs_names) %>% pull(`...1`)
mrpca_plot = plot(out_mrpca)
out_vpca$last_rot_rad_coord2
out_mrpca$last_rot_rad_coord

  out_mrpca$mrpca_inertia


save(file="./Dataset/output_MRPCA_face_recognition.RData",out_mrpca)  

out_mrpca$
out_mrpca$last_rot_rad_coord
