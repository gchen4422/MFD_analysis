args<-as.numeric(commandArgs(TRUE))
LD_BLOCK = args[1]
h2_num = args[2]
num_causal = args[3]

library(data.table)
library(dplyr)

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/functional_annotation_sensitivity/causal_num_",num_causal,"/")
system(paste0("mkdir -p ",wrk_dir))
data_dir<-paste0(wrk_dir,"summary_data/w1_startover/")
system(paste0("mkdir -p ",data_dir))
result_dir<-paste0(wrk_dir,"result/w1_startover/")
system(paste0("mkdir -p ",result_dir))

# LD files are from the original shared_50 simulation (same genotypes, same LD)
LD_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_",num_causal,"/summary_data/")

# Functional annotation directory
annot_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/functional_annotation_sensitivity/annot/")

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

write.table(summary_stat_1_susiex,file = paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_eur"),col.names = T,row.names = F,quote=F, sep ="\t")
write.table(summary_stat_2_susiex,file = paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_afr"),col.names = T,row.names = F,quote=F, sep ="\t")
#write.table(summary_stat_2_bim,file = paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_",LD_BLOCK,".bim"),col.names = F,row.names = F,quote=F, sep ="\t")

summary_stat_sd_list = list("EUR" = summary_stat_1,"AFR"=summary_stat_2 )

EU_cov<-as.matrix(fread(paste0(LD_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")))
colnames(EU_cov)<-zfile$RSID
BB_cov<-as.matrix(fread(paste0(LD_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")))
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
MESuSiE_res<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10,residual_variance=NULL,prior_weights=NULL,ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =T,max_iter =100)
MESuSiE.time<-Sys.time()-start.time
MESuSiE.time<-as.numeric(MESuSiE.time,units = "mins")
MESuSiE_name<-paste0(result_dir,"MESuSiE_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".RData")

#################Run SuSiEx########################

loci_info = fread("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/simulation_locus_info.txt")

chr = loci_info$CHR[LD_BLOCK]
bp_start = loci_info$MinBP[LD_BLOCK]
bp_end = loci_info$MaxBP[LD_BLOCK]

sst_file <- paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_eur,",
                   data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_afr")

ref_file <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur/loci_", LD_BLOCK,
                   ",/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_", LD_BLOCK)

ld_file <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/susiex_ld_eur/loci_", LD_BLOCK,
                  ".ld,/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/susiex_ld_afr/loci_", LD_BLOCK, ".ld")

out_name <- paste0("SuSiEx_CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_output_cs95")

command <- paste(
  "~/gwas_software/SuSiEx/bin/SuSiEx",
  paste0("--sst_file=", sst_file),
  paste0("--ref_file=", ref_file),
  paste0("--ld_file=", ld_file),
  "--n_gwas=300000,300000",
  paste0("--out_dir=", result_dir),
  paste0("--out_name=", out_name),
  paste0("--chr=", chr),
  paste0("--bp=", bp_start, ",", bp_end),
  "--chr_col=1,1",
  "--bp_col=2,2",
  "--snp_col=3,3",
  "--a1_col=4,4",
  "--a2_col=5,5",
  "--eff_col=6,6",
  "--se_col=7,7",
  "--pval_col=10,10",
  "--plink=/home/chen4422/gwas_software/plink2/plink",
  "--keep-ambig=True",
  "--mult-step=True",
  "--maf=0.001",
  "--level=0.95",
  "--min_purity=0.5",
  "--pval_thresh=1e-05",
  "--tol=0.0001",
  "--n_sig=10",
  "--max_iter=200",
  "--threads=16",
  sep = " "
)

start.time<-Sys.time()
system(command)
susiex.time<-Sys.time()-start.time
susiex.time<-as.numeric(susiex.time,units = "mins")

##############################################
#
#		Run SuSiE
#
##############################################
library(susieR)
start.time<-Sys.time()
susie_EU<-susie_rss(zfile$zscore_1,EU_cov,check_prior=F, n = 300000)
susie_BB<-susie_rss(zfile$zscore_2,BB_cov,check_prior=F, n = 300000)

susie.time<-Sys.time()-start.time
susie.time<-as.numeric(susie.time,units = "mins")

##############################################
#
#		Run Paintor (without annotations)
#
##############################################
paintor_dir<-data_dir
paintor_suffix<-paste0("CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num)
input_file_name<-paste0(paintor_dir,paintor_suffix,".input")
input_dir<-paintor_dir

# Write PAINTOR input/LD files to working directory
write.table(paintor_suffix,input_file_name,col.names = F,row.names = F,quote=F)

ld_EU_paintor_name<-paste0(paintor_dir,paintor_suffix,".LD1")
ld_BB_paintor_name<-paste0(paintor_dir,paintor_suffix,".LD2")
fwrite(EU_cov,ld_EU_paintor_name,col.names = F,row.names = F,quote=F,sep=" ")
fwrite(BB_cov,ld_BB_paintor_name,col.names = F,row.names = F,quote=F,sep=" ")

# Write dummy annotation for baseline PAINTOR
annotation_baseline<-matrix(rep(1,nrow(zfile)),ncol=1)
colnames(annotation_baseline)<-"coding"
annotation_baseline_name<-paste0(paintor_dir,paintor_suffix,".annotations")
write.table(annotation_baseline,annotation_baseline_name,col.names = T,row.names = F,quote=F, sep=" ")

setwd(input_dir)
start.time<-Sys.time()
system(paste0("~/gwas_software/PAINTOR_V3.0/PAINTOR -input ",input_file_name," -Zhead  zscore_1,zscore_2 -LDname LD1,LD2 -in ",input_dir," -out ",result_dir," -mcmc -annotations coding -Gname ",paintor_suffix," -Lname ",paste0(paintor_suffix,"_bayes_factor")," -RESname ",paste0("mcmc.paintor")))
Paintor.time<-Sys.time()-start.time
Paintor.time<-as.numeric(Paintor.time,units = "mins")
print(paste0("Paintor time is ",Paintor.time))

##############################################
#
#		Run Paintor + functional annotations
#
##############################################
annot_file_subset<-fread(paste0(annot_dir,"Region_",LD_BLOCK,".annot"))

# Write functional annotation file for PAINTOR
annotation_fun_name<-paste0(paintor_dir,paintor_suffix,".annotations")
write.table(annot_file_subset,annotation_fun_name,col.names = T,row.names = F,quote=F, sep=" ")

annotations = colnames(annot_file_subset)
annotations_part <- paste0("-annotations ", paste(annotations, collapse = ","))

paintor_fun_log <- paste0(result_dir,"Paintor_fun_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_mem.txt")
start.time<-Sys.time()
system(paste0(
  "/usr/bin/time -v -o ", paintor_fun_log, " ",
  "~/gwas_software/PAINTOR_V3.0/PAINTOR -input ",input_file_name,
  " -Zhead zscore_1,zscore_2 -LDname LD1,LD2 -in ",input_dir,
  " -out ",result_dir," -mcmc ", annotations_part,
  " -Gname ", paste0(paintor_suffix, ".enrich"),
  " -Lname ", paste0(paintor_suffix, "_fun_bayes_factor"),
  " -RESname ", paste0("fun.mcmc.paintor")
))
Paintor_fun.time<-Sys.time()-start.time
Paintor_fun.time<-as.numeric(Paintor_fun.time,units = "mins")
print(paste0("Paintor+annot time is ",Paintor_fun.time))

log_content_fun <- readLines(paintor_fun_log)
paintor_fun_ram_kb <- as.numeric(gsub("[^0-9]", "", grep("Maximum resident set size", log_content_fun, value = TRUE)))
paintor_fun_ram_mib <- paintor_fun_ram_kb / 1024

# Remove LD files
#system(paste0("rm ",ld_EU_paintor_name))
#system(paste0("rm ",ld_BB_paintor_name))

##############################################
#
#		Run Paintor + all functional annotations
#
##############################################
annotation_file_name<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_annot/loci_",LD_BLOCK,".annot")
annot_file_full<-fread(annotation_file_name)
evo2_file<-fread(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_annot/CAUSAL_1_LOCI_",LD_BLOCK,"_h2_1_eur_scored.csv"))

# Drop metadata columns, keep only annotation columns + add evo2
annot_file_all <- annot_file_full %>% select(-c("CHR","BP","SNP","CM"))
annot_file_all$evo2_score <- evo2_file$evo2_delta_score
annot_file_all$evo2_score <- rank(-annot_file_all$evo2_score) / nrow(annot_file_all)

# Overwrite .annotations file with all annotations for PAINTOR
write.table(annot_file_all,paste0(paintor_dir,paintor_suffix,".annotations"),col.names = T,row.names = F,quote=F, sep=" ")

annotations_all <- colnames(annot_file_all)
annotations_all_part <- paste0("-annotations ", paste(annotations_all, collapse = ","))

paintor_all_fun_log <- paste0(result_dir,"Paintor_all_fun_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_mem.txt")
start.time<-Sys.time()
system(paste0(
  "/usr/bin/time -v -o ", paintor_all_fun_log, " ",
  "~/gwas_software/PAINTOR_V3.0/PAINTOR -input ",input_file_name,
  " -Zhead zscore_1,zscore_2 -LDname LD1,LD2 -in ",input_dir,
  " -out ",result_dir," -mcmc ", annotations_all_part,
  " -Gname ", paste0(paintor_suffix, ".all_enrich"),
  " -Lname ", paste0(paintor_suffix, "_all_fun_bayes_factor"),
  " -RESname ", paste0("all.fun.mcmc.paintor")
))
Paintor_all_fun.time<-Sys.time()-start.time
Paintor_all_fun.time<-as.numeric(Paintor_all_fun.time,units = "mins")
print(paste0("Paintor+all_annot time is ",Paintor_all_fun.time))

log_content_all_fun <- readLines(paintor_all_fun_log)
paintor_all_fun_ram_kb <- as.numeric(gsub("[^0-9]", "", grep("Maximum resident set size", log_content_all_fun, value = TRUE)))
paintor_all_fun_ram_mib <- paintor_all_fun_ram_kb / 1024

# Remove LD files (after all PAINTOR runs are done)
#system(paste0("rm ",ld_EU_paintor_name))
#system(paste0("rm ",ld_BB_paintor_name))

##############################################
#
#		Run XMAP
#
##############################################

library(XMAP)
z_eur <- zfile$zscore_1
z_afr <- zfile$zscore_2

start.time<-Sys.time()

xmap <-   XMAP(simplify2array(list(EU_cov, BB_cov)), cbind(z_eur, z_afr), n = c(300000, 300000),
               K = 10, tol = 1e-6,
               maxIter = 300, estimate_residual_variance = F, estimate_prior_variance = T,
               estimate_background_variance = F)

xmap.time<-Sys.time()-start.time
xmap.time<-as.numeric(xmap.time,units = "mins")

##############################################
#
#		Run CARMA-X
#
##############################################
gsl_lib <- "/apps/spack/negishi/apps/gsl/2.7.1-gcc-12.2.0-oxu7x5x/lib"
dyn.load(file.path(gsl_lib, "libgslcblas.so.0.0.0"))
dyn.load(file.path(gsl_lib, "libgsl.so.27.0.0"))

library(CARMAX)
library(Matrix)

lambda.power <- 1
z.list <- list()
ld.list <- list()
lambda.list <- c()

g <- 1
z.list[[g]] <- as.matrix(summary_stat_1$`Z`)
lambda.list[[g]] <- lambda.power

g <- 2
z.list[[g]] <- as.matrix(summary_stat_2$`Z`)
lambda.list[[g]] <- lambda.power

ld.list[[1]] <- as.matrix(EU_cov)
ld.list[[2]] <- as.matrix(BB_cov)
ld.list[[3]] <- as.matrix(bdiag(EU_cov, BB_cov))

prior.name <- 'Spike-slab'
set.seed(123)

start.time<-Sys.time()
carmax_results <- CARMAX(z.list, ld.list, lambda.list = lambda.list, effect.size.prior = prior.name, LD.estimation = FALSE, outlier.switch = FALSE)
carmax.time<-Sys.time()-start.time
carmax.time<-as.numeric(carmax.time,units = "mins")

##############################################
#
#		Save results/time used
#
##############################################
save(susie_EU,susie_BB,MESuSiE_res,xmap,carmax_results, file = MESuSiE_name)

time_out_name<-paste0(result_dir,"time_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".txt")
time_vec<-c(MESuSiE.time,susie.time,Paintor.time,Paintor_fun.time,Paintor_all_fun.time,xmap.time,susiex.time,carmax.time)
write.table(time_vec,time_out_name,col.names = F,row.names = F,quote=F)

# Save runtime + memory for Paintor_fun and Paintor_all_fun
perf_df <- data.frame(
  Method = c("Paintor_fun", "Paintor_all_fun"),
  Time_Min = c(Paintor_fun.time, Paintor_all_fun.time),
  Elapsed_Time_sec = c(Paintor_fun.time * 60, Paintor_all_fun.time * 60),
  Peak_RAM_Used_MiB = c(paintor_fun_ram_mib, paintor_all_fun_ram_mib)
)
print(perf_df)
perf_out_name <- paste0(result_dir,"Runtime_memory_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".txt")
write.table(perf_df, perf_out_name, col.names = T, row.names = F, quote = F)
