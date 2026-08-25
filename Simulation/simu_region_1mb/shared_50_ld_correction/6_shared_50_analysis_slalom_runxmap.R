library(Rcpp)
#sourceCpp("/net/fantasia/home/borang/MESuSiE_inf/software/develop_mode/MESuSiE_inf.cpp")
#source("/net/fantasia/home/borang/MESuSiE_inf/software/develop_mode/MESuSiE_inf.R")
library(dplyr)
#library(susieR)
library(data.table)
#source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/modified_kriging_rss.R")


args <- as.numeric(commandArgs(TRUE))

LD_BLOCK <- args[1]
h2_num <- args[2]
num_causal <- args[3]
flip_prop <- args[4]
rho_num <- args[5]

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/causal_num_",num_causal,"/")
system(paste0("mkdir -p ",wrk_dir))
data_dir<-paste0(wrk_dir,"summary_data/")
system(paste0("mkdir -p ",data_dir))
data_dir_flipped = paste0(data_dir,"flipped_",flip_prop,"_rho_",rho_num,"/")
data_dir_flipped_slalom = paste0(data_dir,"flipped_",flip_prop,"_rho_",rho_num,"/slalom/")
result_dir<-paste0(wrk_dir,"result/flipped_",flip_prop,"_rho_",rho_num,"_slalom/")
system(paste0("mkdir -p ",result_dir))


geno_dir_eur <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur")
LD_dir_eur <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_",num_causal,"/summary_data/")

geno_dir_afr <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr")
LD_dir_afr <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_",num_causal,"/summary_data/")



##############################################
#
#		Format MESuSiE Input
#
##############################################
zfile<-read.table(paste0(
  data_dir_flipped,
  "CAUSAL_", num_causal,
  "_LOCI_",  LD_BLOCK,
  "_h2_",    h2_num,
  "_flipped_",flip_prop,
  "_rho_",   rho_num
),header=T)


EU_cov<-as.matrix(fread(paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")))
rownames(EU_cov)<-colnames(EU_cov)<-zfile$RSID
BB_cov<-as.matrix(fread(paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")))
rownames(BB_cov)<-colnames(BB_cov)<-zfile$RSID


##Check for LD mismatch, remove mismatched SNPs
slalom_eur_result = fread(paste0(data_dir_flipped_slalom,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_eur.slalom.txt")) 
slalom_afr_result = fread(paste0(data_dir_flipped_slalom,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_afr.slalom.txt"))

slalom_eur_index <- slalom_eur_result %>%  mutate(is_outlier = r^2 > 0.6 & nlog10p_dentist_s > 4,is_outlier = coalesce(is_outlier, FALSE)  ) %>%  pull(is_outlier)
slalom_afr_index <- slalom_afr_result %>%  mutate(is_outlier = r^2 > 0.6 & nlog10p_dentist_s > 4,is_outlier = coalesce(is_outlier, FALSE)  ) %>%  pull(is_outlier)

slalom_index = slalom_eur_index | slalom_afr_index


zfile_sub = zfile %>% filter(!slalom_index)
EU_cov_sub = EU_cov[!slalom_index,!slalom_index]
BB_cov_sub = BB_cov[!slalom_index,!slalom_index]




MESuSiE_name<-paste0(result_dir,"XMAP_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".RData")




#################Run XMAP########################

library(XMAP)
z_afr <- zfile_sub$zscore_2
z_eur <- zfile_sub$zscore_1

# OmegaHat <- matrix(c(0.01,0,0,0.01), nrow = 2, ncol = 2)

#load("OmegaHat.Rdata")


start.time<-Sys.time()

xmap <- XMAP(simplify2array(list(EU_cov_sub, BB_cov_sub)), cbind(z_eur, z_afr), n = c(300000, 300000),
             K = 10, tol = 1e-6,
             maxIter = 300, estimate_residual_variance = F, estimate_prior_variance = T,
             estimate_background_variance = F)

#xmap.time<-Sys.time()-start.time
#xmap.time<-as.numeric(xmap.time,units = "mins")
#cs1 <- get_CS(xmap, Xcorr = EU_cov, coverage = 0.95, min_abs_corr = 0.5)
#cs2 <- get_CS(xmap, Xcorr = BB_cov, coverage = 0.95, min_abs_corr = 0.5)
#cs_xmap <- cs1$cs[intersect(names(cs1$cs), names(cs2$cs))]
#cs_xmap
#pip_xmap <- get_pip(xmap$gamma)
#plot_CS(pip_xmap, cs_xmap, main = "XMAP")








save(xmap, file = MESuSiE_name)


