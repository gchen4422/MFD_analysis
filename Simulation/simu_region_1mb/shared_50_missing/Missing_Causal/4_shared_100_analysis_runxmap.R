args<-as.numeric(commandArgs(TRUE))
LD_BLOCK = args[1]
h2_num = args[2]
num_causal = args[3]
causal_index = args[4]
External_index = args[5]
causal_index_name = c("Both","One")[causal_index]
External_index_name = c("","External_")[External_index]

library(data.table)
library(dplyr)
num_BB = 300000
num_EU = 300000

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Causal/",External_index_name,causal_index_name,"/causal_num_",num_causal,"/")
system(paste0("mkdir -p ",wrk_dir))

data_dir<-paste0(wrk_dir,"summary_data/")
system(paste0("mkdir -p ",data_dir))

result_dir<-paste0(wrk_dir,"result/")
system(paste0("mkdir -p ",result_dir))

paintor_dir<-data_dir

zfile<-read.table(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/causal_num_",num_causal,"/summary_data/CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),header=T)

SuSiE_res_name<-paste0(result_dir,"MESuSiE_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".RData")
load(SuSiE_res_name)
zfile_MESuSiE_SuSiE_Paintor = zfile_MESuSiE_SuSiE_Paintor %>% select(RSID,EU_missing,BB_missing)
zfile = zfile %>% left_join(zfile_MESuSiE_SuSiE_Paintor, by = "RSID")


EU_cov<-as.matrix(fread(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/causal_num_",num_causal,"/summary_data/CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")))
colnames(EU_cov)<-zfile$RSID
BB_cov<-as.matrix(fread(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/causal_num_",num_causal,"/summary_data/CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")))
colnames(BB_cov)<-zfile$RSID

if(causal_index==1){
  missing_index = which(zfile$EU_missing!=0|zfile$BB_missing!=0)
  zfile_subset<-zfile%>%filter(EU_missing==0,BB_missing==0)
  EU_cov_subset<-EU_cov[-missing_index,-missing_index]
  BB_cov_subset<-BB_cov[-missing_index,-missing_index]
}else if(causal_index==2){
  missing_index = which(zfile$EU_missing!=0)
  zfile_subset<-zfile%>%filter(EU_missing==0)
  EU_cov_subset<-EU_cov[-missing_index,-missing_index]
  BB_cov_subset<-BB_cov[-missing_index,-missing_index]
}




#################Run XMAP########################

library(XMAP)
z_afr <- zfile_subset$zscore_2
z_eur <- zfile_subset$zscore_1

# OmegaHat <- matrix(c(0.01,0,0,0.01), nrow = 2, ncol = 2)

#load("OmegaHat.Rdata")


start.time<-Sys.time()

xmap <- XMAP(simplify2array(list(EU_cov_subset, BB_cov_subset)), cbind(z_eur, z_afr), n = c(300000, 300000),
             K = 10, tol = 1e-6,
             maxIter = 300, estimate_residual_variance = F, estimate_prior_variance = T,
             estimate_background_variance = F)

xmap.time<-Sys.time()-start.time
xmap.time<-as.numeric(xmap.time,units = "mins")
#cs1 <- get_CS(xmap, Xcorr = EU_cov, coverage = 0.95, min_abs_corr = 0.5)
#cs2 <- get_CS(xmap, Xcorr = BB_cov, coverage = 0.95, min_abs_corr = 0.5)
#cs_xmap <- cs1$cs[intersect(names(cs1$cs), names(cs2$cs))]
#cs_xmap
#pip_xmap <- get_pip(xmap$gamma)
#plot_CS(pip_xmap, cs_xmap, main = "XMAP")




add_result_name<-paste0(result_dir,"XMAP_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_additional.RData")

save(xmap, file = add_result_name)


#system(paste0("rm ",ld_EU_paintor_name))
#system(paste0("rm ",ld_BB_paintor_name))
#system(paste0("rm ",ld_EU_name))
#system(paste0("rm ",ld_BB_name))

