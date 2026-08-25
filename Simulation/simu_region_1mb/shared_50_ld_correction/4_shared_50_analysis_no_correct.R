library(Rcpp)
#sourceCpp("/net/fantasia/home/borang/MESuSiE_inf/software/develop_mode/MESuSiE_inf.cpp")
#source("/net/fantasia/home/borang/MESuSiE_inf/software/develop_mode/MESuSiE_inf.R")
library(dplyr)
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
result_dir<-paste0(wrk_dir,"result/flipped_",flip_prop,"_rho_",rho_num,"/")
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

summary_stat_1 = data.frame("SNP" = zfile$RSID, "Beta"=zfile$zscore_1/sqrt(zfile$N_1),"Se"=1/sqrt(zfile$N_1), "Z" =zfile$zscore_1 ,  "N" = zfile$N_1 )
summary_stat_2 = data.frame("SNP" = zfile$RSID, "Beta"=zfile$zscore_2/sqrt(zfile$N_2),"Se"=1/sqrt(zfile$N_2), "Z" =zfile$zscore_2 ,  "N" = zfile$N_2 )
summary_stat_sd_list = list("EUR" = summary_stat_1,"AFR"=summary_stat_2 ) 

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



write.table(summary_stat_1_susiex,file = paste0(data_dir_flipped,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_eur"),col.names = T,row.names = F,quote=F, sep ="\t")
write.table(summary_stat_2_susiex,file = paste0(data_dir_flipped,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_afr"),col.names = T,row.names = F,quote=F, sep ="\t")
#write.table(summary_stat_2_bim,file = paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_",LD_BLOCK,".bim"),col.names = F,row.names = F,quote=F, sep ="\t")



summary_stat_sd_list = list("EUR" = summary_stat_1,"AFR"=summary_stat_2 ) 

EU_cov<-as.matrix(fread(paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")))
colnames(EU_cov)<-zfile$RSID
BB_cov<-as.matrix(fread(paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")))
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
#start.time<-Sys.time()
MESuSiE_res<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10)
#MESuSiE.time<-Sys.time()-start.time
#MESuSiE.time<-as.numeric(MESuSiE.time,units = "mins")
MESuSiE_name<-paste0(result_dir,"MESuSiE_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".RData")

#################Run SuSiEx########################

chr = median(zfile$CHR)
bp_start = min(as.numeric(zfile$POS))
bp_end = max(as.numeric(zfile$POS))


# Construct the sst_file argument
sst_file <- paste0(data_dir_flipped,
                   "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_eur,",data_dir_flipped, 
                   "/CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_afr")

# Construct the ref_file argument
ref_file <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur/loci_", LD_BLOCK, 
                   ",/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_", LD_BLOCK)

# Construct the ld_file argument
ld_file <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/susiex_ld_eur/loci_", LD_BLOCK, 
                  ".ld,/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/susiex_ld_afr/loci_", LD_BLOCK, ".ld")



# Define the out_name variable
out_name <- paste0("SuSiEx_CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_output_cs95")

dir.create(file.path(result_dir, "susiex_result"), recursive = TRUE, showWarnings = FALSE)


# Construct the command
command <- paste(
  "~/gwas_software/SuSiEx/bin/SuSiEx",
  paste0("--sst_file=", sst_file),
  paste0("--ref_file=", ref_file),
  paste0("--ld_file=", ld_file),
  "--n_gwas=300000,300000",
  paste0("--out_dir=",result_dir,"/susiex_result"),
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
  "--pval_thresh=1e-04",
  "--tol=0.0001",
  "--n_sig=10",
  "--max_iter=200",
  "--threads=16",
  sep = " "
)

# Run the command in the system
start.time<-Sys.time()
system(command)
#susiex.time<-Sys.time()-start.time
#susiex.time<-as.numeric(susiex.time,units = "mins")


##############################################
#
#		Run SuSiE
#
##############################################
library(susieR)
start.time<-Sys.time()
susie_EU<-susie_rss(zfile$zscore_1,EU_cov,check_prior=F, n = 300000)
susie_BB<-susie_rss(zfile$zscore_2,BB_cov,check_prior=F, n = 300000)

#susie.time<-Sys.time()-start.time
#susie.time<-as.numeric(susie.time,units = "mins")

#kriging_rss(zfile$zscore_1,EU_cov, n = 300000)
#kriging_rss(zfile$zscore_1,EU_cov, n = 300000)


##############################################
#
#		Run Paintor
#
##############################################

annotation_file<-matrix(rep(1,nrow(zfile)),ncol=1)
colnames(annotation_file)<-"coding"
annotation_paintor_name<-paste0(data_dir_flipped, "CAUSAL_", num_causal,
                                "_LOCI_",  LD_BLOCK,
                                "_h2_",    h2_num,
                                "_flipped_",flip_prop,
                                "_rho_",   rho_num,".annotations")


write.table(annotation_file,annotation_paintor_name,col.names = T,row.names = F,quote=F)

input_file_name<-paste0(data_dir_flipped,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".input")
write.table(paste0("CAUSAL_", num_causal,
                     "_LOCI_",  LD_BLOCK,
                     "_h2_",    h2_num,
                     "_flipped_",flip_prop,
                     "_rho_",   rho_num),input_file_name,col.names = F,row.names = F,quote=F)




paintor_dir = data_dir_flipped
paintor_out_name<-paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num)
input_file_name<-paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".input")
ld_EU_paintor_name<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")
ld_BB_paintor_name<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")
paintor_suffix<-paste0("CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num)
input_dir<-paintor_dir
start.time<-Sys.time()
system(paste0("~/gwas_software/PAINTOR_V3.0/PAINTOR -input ",input_file_name," -Zhead  zscore_1,zscore_2 -LDname LD1,LD2 -in ",input_dir," -out ",result_dir," -mcmc "," -annotations coding -Gname ",paintor_suffix," -Lname ",paste0(paintor_suffix,"_bayes_factor")," -RESname ",paste0("mcmc.paintor")))
#Paintor.time<-Sys.time()-start.time
#Paintor.time<-as.numeric(Paintor.time,units = "mins")
#print(paste0("Paintor time is ",Paintor.time))

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




#################Run MultiSuSiE########################

## Use python scripts in the home directory ~



#################Run CARMA-X########################
gsl_lib <- "/apps/spack/negishi/apps/gsl/2.7.1-gcc-12.2.0-oxu7x5x/lib"
# Preload GSL dependencies (order matters: cblas first, then gsl)
dyn.load(file.path(gsl_lib, "libgslcblas.so.0.0.0"))
dyn.load(file.path(gsl_lib, "libgsl.so.27.0.0"))

library(CARMAX)
library(Matrix)



# Initialize lists
lambda.power <- 1
z.list <- list()
ld.list <- list()
full.w.list <- list()
lambda.list <- c()
label.list <- list()
input.prior.list <- list()
# Read the Z-scores 
# AFR ancestry
g <- 1
z.list[[g]] <- as.matrix(summary_stat_1$`Z`)
lambda.list[[g]] <- lambda.power

# EUR ancestry
g <- 2
z.list[[g]] <- as.matrix(summary_stat_2$`Z`)
lambda.list[[g]] <- lambda.power


ld.list[[1]] <- as.matrix(EU_cov)
ld.list[[2]] <- as.matrix(BB_cov)
ld.list[[3]] <- as.matrix(bdiag(EU_cov, BB_cov))


# Define the prior distribution
prior.name <- 'Spike-slab'

# Set a random seed
set.seed(123)
# Run CARMA-X
start.time<-Sys.time()

carmax_results <- CARMAX(z.list, ld.list, lambda.list = lambda.list, effect.size.prior = prior.name, LD.estimation = FALSE, outlier.switch = FALSE)

#carmax.time<-Sys.time()-start.time
#carmax.time<-as.numeric(carmax.time,units = "mins")




save(susie_EU,susie_BB,MESuSiE_res,xmap,carmax_results, file = MESuSiE_name)
