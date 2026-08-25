library(Rcpp)
#sourceCpp("/net/fantasia/home/borang/MESuSiE_inf/software/develop_mode/MESuSiE_inf.cpp")
#source("/net/fantasia/home/borang/MESuSiE_inf/software/develop_mode/MESuSiE_inf.R")
library(dplyr)
library(susieR)
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
result_dir<-paste0(wrk_dir,"result/flipped_",flip_prop,"_rho_",rho_num,"_susie_rss/")
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
snplist_filtered <- zfile$RSID
idx_to_remove = c()
while(!(identical(idx_to_remove, integer(0)))){
  
  zfile_sub = zfile %>% filter(RSID %in% snplist_filtered)
  EU_cov_sub = EU_cov[snplist_filtered,snplist_filtered]
  BB_cov_sub = BB_cov[snplist_filtered,snplist_filtered]
  
  
  EU_diagnostic <- kriging_rss(zfile_sub$zscore_1, EU_cov_sub, n = median(zfile_sub$N_1))
  idx_to_remove_eur = EU_diagnostic$plot$plot_env$idx
  BB_diagnostic <- kriging_rss(zfile_sub$zscore_2, BB_cov_sub, n = median(zfile_sub$N_2))
  idx_to_remove_afr = BB_diagnostic$plot$plot_env$idx
  idx_to_remove = sort(union(idx_to_remove_eur,idx_to_remove_afr))
  
  if (!(identical(idx_to_remove, integer(0)))){
    snplist_filtered = snplist_filtered[-idx_to_remove]
  }
  
}



# helper to make per-ancestry summary stats
make_sumstats <- function(df, z_col, n_col){
  df %>%
    transmute(
      SNP  = RSID,
      Beta = .data[[z_col]] / sqrt(.data[[n_col]]),
      Se   = 1 / sqrt(.data[[n_col]]),
      Z    = .data[[z_col]],
      N    = .data[[n_col]]
    )
}

# 1) summary stats
summary_stat_1 <- make_sumstats(zfile_sub, "zscore_1", "N_1")
summary_stat_2 <- make_sumstats(zfile_sub, "zscore_2", "N_2")

summary_stat_sd_list <- list(EUR = summary_stat_1, AFR = summary_stat_2)

# 2) allele info from BIM (only SNP + A1/A2)
bim_eur <- fread(
  paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur/loci_", LD_BLOCK, ".bim"),
  header = FALSE
) %>%
  transmute(SNP = V2, A1 = V5, A2 = V6)

bim_afr <- fread(
  paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_", LD_BLOCK, ".bim"),
  header = FALSE
) %>%
  transmute(SNP = V2, A1 = V5, A2 = V6)

# 3) base columns to carry over from zfile_sub
#    (replaces zfile_sub[, c(1,2,3)] but keeps names)
base_cols <- zfile_sub %>%
  select(1, 2, 3) %>%      # retain the first three columns
  mutate(SNP = RSID) %>%
  select(-RSID)

# 4) build SuSiEx input tables with joins
summary_stat_1_susiex <- base_cols %>%
  left_join(bim_eur, by = "SNP") %>%
  left_join(summary_stat_1, by = "SNP") %>%
  mutate(P = 2 * pnorm(abs(Z), lower.tail = FALSE)) %>% dplyr::rename(RSID = "SNP")

summary_stat_2_susiex <- base_cols %>%
  left_join(bim_afr, by = "SNP") %>%
  left_join(summary_stat_2, by = "SNP") %>%
  mutate(P = 2 * pnorm(abs(Z), lower.tail = FALSE)) %>% dplyr::rename(RSID = "SNP")

# colnames(summary_stat_1_susiex)[c(3,4,5)] = c("SNP","REF","ALT")
# colnames(summary_stat_2_susiex)[c(3,4,5)] = c("SNP","REF","ALT")



write.table(summary_stat_1_susiex,file = paste0(data_dir_flipped,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_susie_rss_eur"),col.names = T,row.names = F,quote=F, sep ="\t")
write.table(summary_stat_2_susiex,file = paste0(data_dir_flipped,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_susie_rss_afr"),col.names = T,row.names = F,quote=F, sep ="\t")
#write.table(summary_stat_2_bim,file = paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_",LD_BLOCK,".bim"),col.names = F,row.names = F,quote=F, sep ="\t")
z_file_out <- paste0(data_dir_flipped, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_flipped_", flip_prop, "_rho_", rho_num, "_susie_rss")
write.table(zfile_sub,z_file_out,col.names = T,row.names = F,quote=F,sep=" ")  
ld_EU_paintor_name<-paste0(data_dir_flipped,"CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_flipped_", flip_prop, "_rho_", rho_num,"_susie_rss.LD1")
ld_BB_paintor_name<-paste0(data_dir_flipped,"CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_flipped_", flip_prop, "_rho_", rho_num,"_susie_rss.LD2")
fwrite(EU_cov_sub,ld_EU_paintor_name,sep=" ",col.names=F)
fwrite(BB_cov_sub,ld_BB_paintor_name,sep=" ",col.names=F)


summary_stat_sd_list = list("EUR" = summary_stat_1,"AFR"=summary_stat_2 ) 
R_mat_list=list("EUR" = EU_cov_sub,"AFR" = BB_cov_sub)



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

chr = median(zfile_sub$CHR)
bp_start = min(as.numeric(zfile_sub$POS))
bp_end = max(as.numeric(zfile_sub$POS))


# Construct the sst_file argument
sst_file <- paste0(data_dir_flipped,
                   "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_susie_rss_eur,",data_dir_flipped, 
                   "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_susie_rss_afr")

# Construct the ref_file argument
ref_file <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur/loci_", LD_BLOCK, 
                   ",/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_", LD_BLOCK)

# Construct the ld_file argument
ld_file <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/susiex_ld_eur/loci_", LD_BLOCK, 
                  ".ld,/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/susiex_ld_afr/loci_", LD_BLOCK, ".ld")



# Define the out_name variable
out_name <- paste0("SuSiEx_CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_output_cs95_susie_rss")

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
  "--pval_thresh=1e-05",
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
susie_EU<-susie_rss(zfile_sub$zscore_1,EU_cov_sub,check_prior=F, n = 300000)
susie_BB<-susie_rss(zfile_sub$zscore_2,BB_cov_sub,check_prior=F, n = 300000)

#susie.time<-Sys.time()-start.time
#susie.time<-as.numeric(susie.time,units = "mins")

#kriging_rss(zfile$zscore_1,EU_cov, n = 300000)
#kriging_rss(zfile$zscore_1,EU_cov, n = 300000)


##############################################
#
#		Run Paintor
#
##############################################

annotation_file<-matrix(rep(1,nrow(zfile_sub)),ncol=1)
colnames(annotation_file)<-"coding"
annotation_paintor_name<-paste0(data_dir_flipped, "CAUSAL_", num_causal,
                                "_LOCI_",  LD_BLOCK,
                                "_h2_",    h2_num,
                                "_flipped_",flip_prop,
                                "_rho_",   rho_num,"_susie_rss.annotations")


write.table(annotation_file,annotation_paintor_name,col.names = T,row.names = F,quote=F)

input_file_name<-paste0(data_dir_flipped,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_susie_rss.input")
write.table(paste0("CAUSAL_", num_causal,
                   "_LOCI_",  LD_BLOCK,
                   "_h2_",    h2_num,
                   "_flipped_",flip_prop,
                   "_rho_",   rho_num,"_susie_rss"),input_file_name,col.names = F,row.names = F,quote=F)




paintor_dir = data_dir_flipped
paintor_out_name<-paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_susie_rss")
input_file_name<-paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_susie_rss.input")
ld_EU_paintor_name<-paste0(data_dir_flipped,"CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_flipped_", flip_prop, "_rho_", rho_num,"_susie_rss.LD1")
ld_BB_paintor_name<-paste0(data_dir_flipped,"CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_flipped_", flip_prop, "_rho_", rho_num,"_susie_rss.LD2")
paintor_suffix<-paste0("CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num)
input_dir<-data_dir_flipped
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
z_eur <- zfile_sub$zscore_1
z_afr <- zfile_sub$zscore_2

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


ld.list[[1]] <- as.matrix(EU_cov_sub)
ld.list[[2]] <- as.matrix(BB_cov_sub)
ld.list[[3]] <- as.matrix(bdiag(EU_cov_sub, BB_cov_sub))


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
