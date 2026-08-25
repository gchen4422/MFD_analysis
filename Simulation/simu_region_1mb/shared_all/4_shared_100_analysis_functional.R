args<-as.numeric(commandArgs(TRUE))
LD_BLOCK = args[1]
h2_num = args[2]
num_causal = args[3]

library(data.table)
library(dplyr)

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/causal_num_",num_causal,"_annot/")
system(paste0("mkdir -p ",wrk_dir))
data_dir<-paste0(wrk_dir,"summary_data/")
system(paste0("mkdir -p ",data_dir))
result_dir<-paste0(wrk_dir,"result/")
system(paste0("mkdir -p ",result_dir))
paintor_dir<-data_dir


##############################################
#
#		Format MESuSiE Input
#
##############################################
zfile<-read.table(paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),header=T)
summary_stat_1 = data.frame("SNP" = zfile$RSID, "Beta"=zfile$zscore_1/sqrt(zfile$N_1),"Se"=1/sqrt(zfile$N_1), "Z" =zfile$zscore_1 ,  "N" = zfile$N_1 )
summary_stat_2 = data.frame("SNP" = zfile$RSID, "Beta"=zfile$zscore_2/sqrt(zfile$N_2),"Se"=1/sqrt(zfile$N_2), "Z" =zfile$zscore_2 ,  "N" = zfile$N_2 )

summary_stat_1_bim = read.table(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur/loci_",LD_BLOCK,".bim"))
summary_stat_2_bim = read.table(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_",LD_BLOCK,".bim"))
summary_stat_2_bim$V2 = summary_stat_1_bim$V2


summary_stat_1_allele = data.frame("A1" = summary_stat_1_bim$V5, "A2" = summary_stat_1_bim$V6)
summary_stat_2_allele = data.frame("A1" = summary_stat_2_bim$V5, "A2" = summary_stat_2_bim$V6)

summary_stat_1_susiex = cbind(zfile[,c(1,2,3)],summary_stat_1_allele,summary_stat_1[,c(2:5)])
summary_stat_2_susiex = cbind(zfile[,c(1,2,3)],summary_stat_2_allele,summary_stat_2[,c(2:5)])

summary_stat_1_susiex$P = 2*pnorm(abs(summary_stat_1_susiex$Z),lower.tail = F)
summary_stat_2_susiex$P = 2*pnorm(abs(summary_stat_2_susiex$Z),lower.tail = F)

# colnames(summary_stat_1_susiex)[c(3,4,5)] = c("SNP","REF","ALT")
# colnames(summary_stat_2_susiex)[c(3,4,5)] = c("SNP","REF","ALT")



prior_weights = fread(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/polyfun_doc/polyfun_output_eur/loci_",LD_BLOCK,"_eur_output.txt"))$SNPVAR
prior_weights = prior_weights/sum(prior_weights)

if(length(prior_weights) != nrow(zfile)){
  print("Missing functional annotation for SNPs")
}


write.table(summary_stat_1_susiex,file = paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_eur"),col.names = T,row.names = F,quote=F, sep ="\t")
write.table(summary_stat_2_susiex,file = paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_afr"),col.names = T,row.names = F,quote=F, sep ="\t")
write.table(summary_stat_2_bim,file = paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_",LD_BLOCK,".bim"),col.names = F,row.names = F,quote=F, sep ="\t")



summary_stat_sd_list = list("EUR" = summary_stat_1,"AFR"=summary_stat_2 ) 

EU_cov<-as.matrix(fread(paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")))
colnames(EU_cov)<-zfile$RSID
BB_cov<-as.matrix(fread(paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")))
colnames(BB_cov)<-zfile$RSID
R_mat_list=list("EUR" = EU_cov,"AFR" = BB_cov)
##############################################
#
#		Run MESuSiE
#
##############################################

library(MESuSiE)
bayes_fac=3
ancestry_weight=c(bayes_fac/(bayes_fac*2+1),bayes_fac/(bayes_fac*2+1),1/(bayes_fac*2+1))
start.time<-Sys.time()
MESuSiE_res<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10,residual_variance=NULL,prior_weights=prior_weights,ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =T,max_iter =100)
MESuSiE.time<-Sys.time()-start.time
MESuSiE.time<-as.numeric(MESuSiE.time,units = "mins")
MESuSiE_name<-paste0(result_dir,"MESuSiE_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".RData")

#################Run SuSiEx########################

# loci_info = fread("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/simulation_locus_info.txt")
# 
# chr = loci_info$CHR[LD_BLOCK]
# bp_start = loci_info$MinBP[LD_BLOCK]
# bp_end = loci_info$MaxBP[LD_BLOCK]
# 
# 
# # Construct the sst_file argument
# sst_file <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/causal_num_", num_causal,
#                    "/summary_data/CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_eur,/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/causal_num_",
#                    num_causal, "/summary_data/CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_afr")
# 
# # Construct the ref_file argument
# ref_file <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur/loci_", LD_BLOCK,
#                    ",/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_", LD_BLOCK)
# 
# # Construct the ld_file argument
# ld_file <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/susiex_ld_eur/loci_", LD_BLOCK,
#                   ".ld,/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/susiex_ld_afr/loci_", LD_BLOCK, ".ld")
# 
# # Define the out_name variable
# out_name <- paste0("SuSiEx_CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_output_cs95")
# 
# 
# # Construct the command
# command <- paste(
#   "~/gwas_software/SuSiEx/bin/SuSiEx",
#   paste0("--sst_file=", sst_file),
#   paste0("--ref_file=", ref_file),
#   paste0("--ld_file=", ld_file),
#   "--n_gwas=300000,300000",
#   "--out_dir=./susiex_result",
#   paste0("--out_name=", out_name),
#   paste0("--chr=", chr),
#   paste0("--bp=", bp_start, ",", bp_end),
#   "--chr_col=1,1",
#   "--bp_col=2,2",
#   "--snp_col=3,3",
#   "--a1_col=4,4",
#   "--a2_col=5,5",
#   "--eff_col=6,6",
#   "--se_col=7,7",
#   "--pval_col=10,10",
#   "--plink=/home/chen4422/gwas_software/plink2/plink",
#   "--keep-ambig=True",
#   "--mult-step=True",
#   "--maf=0.001",
#   "--level=0.95",
#   "--min_purity=0.5",
#   "--pval_thresh=1e-05",
#   "--tol=0.0001",
#   "--n_sig=10",
#   "--max_iter=200",
#   "--threads=16",
#   sep = " "
# )
# 
# # Run the command in the system
# start.time<-Sys.time()
# system(command)
# susiex.time<-Sys.time()-start.time
# susiex.time<-as.numeric(susiex.time,units = "mins")

##############################################
#
#		Run SuSiE
#
##############################################
library(susieR)
start.time<-Sys.time()
susie_EU<-susie_rss(zfile$zscore_1,EU_cov,check_prior=F,prior_weights = prior_weights, n = 300000)
susie_BB<-susie_rss(zfile$zscore_2,BB_cov,check_prior=F,prior_weights = prior_weights, n = 300000)



# susie_EU_test<-susie_rss(zfile$zscore_1,EU_cov,check_prior=F, n = 300000)
# susie_BB_test<-susie_rss(zfile$zscore_2,BB_cov,check_prior=F, n = 300000)

# $L1
# [1] 476
# 
# $L2
# [1] 982
# 
# $L4
# [1] 1557
# 
# $L3
# [1] 437 461 507 508



# susie_EU_test$sets$cs
# susie_BB_test$sets$cs


# $L1
# [1] 476
# 
# $L4
# [1] 1557
# 
# $L3
# [1] 437 461 507 508
# 
# $L2
# [1]  982 1004 1005

## Functional annotation priors can reduce credible-set size.
##prioritizing the true causal snp


susie.time<-Sys.time()-start.time
susie.time<-as.numeric(susie.time,units = "mins")

#kriging_rss(zfile$zscore_1,EU_cov, n = 300000)
#kriging_rss(zfile$zscore_2,BB_cov, n = 300000)


##############################################
#
#		Run Paintor
#
##############################################

paintor_out_name<-paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num)
input_file_name<-paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".input")
ld_EU_paintor_name<-paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")
ld_BB_paintor_name<-paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")
paintor_suffix<-paste0("CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num)
input_dir<-paintor_dir

load("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/PAINTOR_top_annotation.Rdata")
loci_annot = fread(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_annot/loci_",LD_BLOCK,".annot")) %>% select(-c(1:4))
loci_annot = loci_annot %>% select(top_annotations_list[[LD_BLOCK]]) %>% mutate(base = 1)
file_path <- paste0(paintor_out_name,".annotations")
write.table(loci_annot, file = file_path, row.names = FALSE, col.names = T, quote = FALSE)


start.time<-Sys.time()
system(paste0("~/gwas_software/PAINTOR_V3.0/PAINTOR -input ",input_file_name," -Zhead  zscore_1,zscore_2 -LDname LD1,LD2 -in ",input_dir," -out ",result_dir," -mcmc ",paste0(" -annotations ",paste(colnames(loci_annot), collapse = ",")," -Gname "),paintor_suffix," -Lname ",paste0(paintor_suffix,"_bayes_factor")," -RESname ",paste0("mcmc.paintor")))

Paintor.time<-Sys.time()-start.time
Paintor.time<-as.numeric(Paintor.time,units = "mins")
print(paste0("Paintor time is ",Paintor.time))


#Remove LD files
#system(paste0("rm ",ld_EU_paintor_name))
#system(paste0("rm ",ld_BB_paintor_name))


#################Run XMAP########################

library(XMAP)
z_eur <- zfile$zscore_1
z_afr <- zfile$zscore_2

# OmegaHat <- matrix(c(0.01,0,0,0.01), nrow = 2, ncol = 2)

#load("OmegaHat.Rdata")


start.time<-Sys.time()

xmap <- XMAP(simplify2array(list(EU_cov, BB_cov)), cbind(z_eur, z_afr), n = c(300000, 300000),
             K = 10, tol = 1e-6, prior_weights = prior_weights,
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



save(susie_EU,susie_BB,MESuSiE_res,xmap, file = MESuSiE_name)




#################Run MultiSuSiE########################





##############################################
#
#		Save time used
#
##############################################
time_out_name<-paste0(result_dir,"time_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".txt")
time_vec<-c(MESuSiE.time,susie.time,Paintor.time, xmap.time)
write.table(time_vec,time_out_name,col.names = F,row.names = F,quote=F)
