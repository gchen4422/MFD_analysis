args <- as.numeric(commandArgs(TRUE))
if(length(args) < 2) stop("Please provide Block ID and Trait Index arguments.")
current_block_id <- args[1]
trait_idx <- args[2]

# Preload GSL dependencies
gsl_lib <- "/apps/spack/negishi/apps/gsl/2.7.1-gcc-12.2.0-oxu7x5x/lib"
dyn.load(file.path(gsl_lib, "libgslcblas.so.0.0.0"))
dyn.load(file.path(gsl_lib, "libgsl.so.27.0.0"))




library(dplyr)
library(data.table)
library(MESuSiE)
library(susieR)
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Causal/MFD_utility.R")

traits_list <- c("BMI", "DBP", "SBP")
trait <- traits_list[trait_idx] # Select trait based on argument


# Directories
wrk_dir <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/results_", trait, "/")
data_dir <- paste0(wrk_dir, "summary_data/")
result_dir <- paste0(wrk_dir, "result/")
out_dir <- paste0(wrk_dir, "out/")

# Read Candidate Regions
candidate_region_file <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/", trait, "_risk_loci")
candidate_region <- fread(candidate_region_file)


## LOCI 31 (rs13400734),148,190,200 (rs7279010) for BMI is an example of how susie-post hoc can find AS-V but mesusie cannot, it is filtered out by taking common subset
## LOCI 11 (rs9821489) for DBP is an example of how susie-post hoc can find AS-V but mesusie cannot, it is filtered out by taking common subset

## LOCI 190 is an example of how susie-post hoc can find AS-V but mesusie find potential tagged snps, it is filtered out by taking common subset
## LOCI 115 (rs146787607)



# Use tryCatch to skip regions where files might be missing without stopping the whole script

cat(paste0("Processing Block: ", current_block_id, "... "))



# Define file paths for this iteration
eur_sum_file <- paste0(data_dir, "eur_loci_", current_block_id, "_summ")
afr_sum_file <- paste0(data_dir, "afr_loci_", current_block_id, "_summ")
eur_ld_file <- paste0(data_dir, "eur_loci_", current_block_id, ".ld")
afr_ld_file <- paste0(data_dir, "afr_loci_", current_block_id, ".ld")

# Check if files exist
if(!file.exists(eur_sum_file) | !file.exists(afr_sum_file)) {
  stop("Summary statistics file missing")
}

# Load Summary Stats
summary_stat_eur <- fread(eur_sum_file)
colnames(summary_stat_eur)[which(colnames(summary_stat_eur) == "BETA")] <- "Beta"
colnames(summary_stat_eur)[which(colnames(summary_stat_eur) == "SE")] <- "Se"

summary_stat_afr <- fread(afr_sum_file)
colnames(summary_stat_afr)[which(colnames(summary_stat_afr) == "BETA")] <- "Beta"
colnames(summary_stat_afr)[which(colnames(summary_stat_afr) == "SE")] <- "Se"

# Load LD Matrices
EU_cov <- as.matrix(fread(eur_ld_file))
BB_cov <- as.matrix(fread(afr_ld_file))

# Ensure Matrix Dimensions match Summary Stats (Basic error checking)
# SNP order is established by the upstream preparation step
rownames(EU_cov) <- summary_stat_eur$SNP
colnames(EU_cov) <- summary_stat_eur$SNP
rownames(BB_cov) <- summary_stat_afr$SNP
colnames(BB_cov) <- summary_stat_afr$SNP



##############################################
#
#		Run MF Decision tree
#
##############################################


mfd_result = run_mf_decision_fm(summary_stat_eur,summary_stat_afr,EU_cov,BB_cov, pop_names = c("EUR", "AFR"))

mfd_result = run_mf_decision(summary_stat_eur,summary_stat_afr,EU_cov,BB_cov, pop_names = c("EUR", "AFR"))

mfd_result_ori = mfd_result
##############################################
#
#		Run SuSiE
#
##############################################
susie_EU<-susie_rss(summary_stat_eur$Z,EU_cov, n = median(summary_stat_eur$N))
susie_BB<-susie_rss(summary_stat_afr$Z,BB_cov, n = median(summary_stat_afr$N))


susie_BB_cs<-rep(0,length(susie_BB$pip))
susie_BB_cs[unlist(susie_BB$sets$cs)]<-1
susie_EU_cs<-rep(0,length(susie_EU$pip))
susie_EU_cs[unlist(susie_EU$sets$cs)]<-1

SuSiE_EU_res_dataframe<-data.frame(SNP = summary_stat_eur$SNP,SuSiE_PIP_EU = susie_EU$pip,susie_EU_cs = susie_EU_cs )
SuSiE_BB_res_dataframe<-data.frame(SNP = summary_stat_afr$SNP,SuSiE_PIP_BB = susie_BB$pip,susie_BB_cs = susie_BB_cs )


mfd_result_all <- mfd_result %>%
  left_join(SuSiE_EU_res_dataframe, by = "SNP")%>%
  left_join(SuSiE_BB_res_dataframe, by = "SNP")

mfd_result_all <- mfd_result_all%>% replace(is.na(.), 0)
mfd_result_all <- mfd_result_all%>%mutate(SuSiE_cs = ifelse(susie_EU_cs+susie_BB_cs==0,0,1),SuSiE_PIP_Either = pmax(SuSiE_PIP_EU,SuSiE_PIP_BB),SuSiE_PIP_Shared = pmin(SuSiE_PIP_EU,SuSiE_PIP_BB))




##############################################
#
#		Run MESuSiE
#
##############################################

# MESuSiE requires variants present in both datasets
common_snps <- intersect(summary_stat_eur$SNP, summary_stat_afr$SNP)

g1_sub <- summary_stat_eur[SNP %in% common_snps][order(match(SNP, common_snps))]
g2_sub <- summary_stat_afr[SNP %in% common_snps][order(match(SNP, common_snps))]

# Subset LD matrices based on SNP IDs (assuming LD has row/col names)
ld1_sub <- EU_cov[common_snps, common_snps]
ld2_sub <- BB_cov[common_snps, common_snps]

summary_stat_sd_list = list(g1_sub, g2_sub)
names(summary_stat_sd_list) <- c("EUR", "AFR")
R_mat_list = list(ld1_sub, ld2_sub)
names(R_mat_list) <- c("EUR", "AFR")

mesusie_res <- meSuSie_core(R_mat_list, summary_stat_sd_list, L = 10)

MESuSiE_cs<-rep(0,length(mesusie_res$pip))
MESuSiE_cs[unlist(mesusie_res$cs$cs)]<-1
MESuSiE_res_dataframe<-data.frame(SNP = g1_sub$SNP,MESuSiE_PIP_Either = mesusie_res$pip,MESuSiE_PIP_WB =mesusie_res$pip_config[,1] ,MESuSiE_PIP_BB =mesusie_res$pip_config[,2] ,MESuSiE_PIP_Shared =mesusie_res$pip_config[,3],MESuSiE_cs = MESuSiE_cs )

mfd_result_all <- mfd_result_all %>%
  left_join(MESuSiE_res_dataframe, by = "SNP")%>% replace(is.na(.), 0)

##############################################
#
#		Run SuSiEx
#
##############################################

chr = g1_sub$CHR[current_block_id]
bp_start = min(g1_sub$POS)
bp_end = max(g1_sub$POS)

output_file <- paste0(data_dir, "susiex_input/eur_loci_", current_block_id, "_summ")
if(!dir.exists(dirname(output_file))) {
  dir.create(dirname(output_file), recursive = TRUE)
}
fwrite(g1_sub, file = output_file, quote = FALSE, row.names = FALSE, col.names = TRUE,sep = "\t")

output_file <- paste0(data_dir, "susiex_input/afr_loci_", current_block_id, "_summ")
fwrite(g2_sub,file = paste0(data_dir,"susiex_input/afr_loci_",current_block_id,"_summ"), quote = FALSE, row.names=FALSE, col.names = T,sep = "\t")



# Construct the sst_file argument
sst_file <- paste0(data_dir,"susiex_input/eur_loci_", current_block_id, "_summ,",data_dir, "susiex_input/afr_loci_",current_block_id, "_summ")


susiex_ref_eur <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/risk_loci_ld_eur_",trait)
susiex_ref_afr <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/risk_loci_ld_afr_",trait)


# Construct the ref_file argument
ref_file <- paste0(susiex_ref_eur,"/loci_", current_block_id, 
                   ",",susiex_ref_afr,"/loci_", current_block_id)

susiex_out<- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/results_",trait,"/result/susiex_result/")
if(!dir.exists(dirname(susiex_out))) {
  dir.create(dirname(susiex_out), recursive = TRUE)
}

# Construct the ld_file argument
ld_file <- paste0(susiex_out, "/eur_loci_", current_block_id, 
                  ".ld,",susiex_out,"/afr_loci_", current_block_id, ".ld")

# Define the out_name variable
out_name <- paste0("SuSiEx_",current_block_id,"_output_cs95")

n_gwas_chr = paste(median(g1_sub$N),median(g2_sub$N),sep = ",")

# Construct the command
command <- paste(
  "~/gwas_software/SuSiEx/bin/SuSiEx",
  paste0("--sst_file=", sst_file),
  paste0("--ref_file=", ref_file),
  paste0("--ld_file=", ld_file),
  paste0("--n_gwas=", n_gwas_chr),
  paste0("--out_dir=", susiex_out),
  paste0("--out_name=", out_name),
  paste0("--chr=", chr),
  paste0("--bp=", bp_start, ",", bp_end),
  "--chr_col=2,2",
  "--bp_col=3,3",
  "--snp_col=1,1",
  "--a1_col=5,5",
  "--a2_col=6,6",
  "--eff_col=9,9",
  "--se_col=10,10",
  "--pval_col=4,4",
  "--plink=/home/chen4422/gwas_software/plink2/plink",
  "--keep-ambig=True",
  "--mult-step=True",
  "--maf=0.01",
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
system(command)


susiex_res_name<-paste0(susiex_out,out_name,".snp")
susiex_cred_name<-paste0(susiex_out,out_name,".cs")

if (file.exists(susiex_res_name) && file.exists(susiex_cred_name)) {
  susiex_res = read.table(susiex_res_name, header = T)
  susiex_cs = read.table(susiex_cred_name,header = T)
  #susiex_cs_num = max(susiex_cs$CS_ID)
  susiex_cred = which(mfd_result_all$POS %in% susiex_cs$BP)
  
  # df: data.frame with per-CS PIP columns for all SNPs
  susiex_pip_cols <- grep("^PIP(\\(CS[0-9]+\\)|\\.CS[0-9]+\\.)$", names(susiex_res), value = TRUE)
  
  # fast matrix calc with stability (row-wise)
  M <- as.matrix(susiex_res[, susiex_pip_cols, drop = FALSE])
  susiex_res$pip <- 1 - apply(1 - M, 1, prod)  # per row/SNP overall pip
  susiex_res_to_merge = susiex_res %>% select(SNP,pip) %>% mutate(SuSiEx_cs = 0) %>% dplyr::rename(SuSiEx_PIP_Either = "pip")
  
  mfd_result_all <- mfd_result_all %>%
    left_join(susiex_res_to_merge, by = "SNP")
  
  mfd_result_all <- mfd_result_all%>% replace(is.na(.), 0)
  mfd_result_all$SuSiEx_cs[susiex_cred] <- 1
  
}else {
  mfd_result_all$SuSiEx_PIP_Either = 0
  mfd_result_all$SuSiEx_cs  = 0
}
##############################################
#
#		Run Paintor
#
##############################################
system(paste0("mkdir -p ",paste0(data_dir,"paintor_all/")))
paintor_dir<-paste0(data_dir,"paintor_all/")

input_file_name<-paste0(paintor_dir,"loci_",current_block_id,".input")
write.table(paste0("loci_",current_block_id),input_file_name,col.names = F,row.names = F,quote=F)
sumstats = cbind(g1_sub[,c(1,8)],g2_sub[,c(8)])
colnames(sumstats) = c("SNP","zscore_1","zscore_2")
fwrite(sumstats,paste0(paintor_dir,"loci_",current_block_id),col.names = T,row.names = F,quote=F,sep=" ")


ld_EU_paintor_name<-paste0(paintor_dir,"loci_",current_block_id,".LD1")
ld_BB_paintor_name<-paste0(paintor_dir,"loci_",current_block_id,".LD2")
fwrite(ld1_sub,ld_EU_paintor_name,col.names = F,row.names = F,quote=F,sep=" ")
fwrite(ld2_sub,ld_BB_paintor_name,col.names = F,row.names = F,quote=F,sep=" ")


paintor_suffix<-paste0("loci_",current_block_id)
input_dir<-paintor_dir
setwd(input_dir)
annotation_file<-matrix(rep(1,nrow(g1_sub)),ncol=1)
colnames(annotation_file)<-"coding"
annotation_paintor_name<-paste0(paintor_dir,"loci_",current_block_id,".annotations")
write.table(annotation_file,annotation_paintor_name,col.names = T,row.names = F,quote=F)
system(paste0("~/gwas_software/PAINTOR_V3.0/PAINTOR -input ",input_file_name," -Zhead  zscore_1,zscore_2 -LDname LD1,LD2 -in ",input_dir," -out ",result_dir," -mcmc "," -annotations coding -Gname ", paintor_suffix," -Lname ",paste0(paintor_suffix,"_bayes_factor")," -RESname ",paste0("mcmc.paintor")))


find_index_until_threshold <- function(vec, threshold) {
  # Sort the input vector in decreasing order and obtain the sorted indices
  sorted_indices <- order(vec, decreasing = TRUE)
  
  # Calculate the cumulative sum of the sorted vector
  cumsum_vec <- cumsum(vec[sorted_indices])
  
  # Find the index where the cumulative sum reaches or exceeds the threshold
  index <- which(cumsum_vec >= threshold)[1]
  
  # Return the original indices of the elements up to the found index
  sorted_indices[1:index]
}


paintor_res<-read.table(paste0(result_dir,paintor_suffix,".mcmc.paintor"),header=T)

if(sum(paintor_res$Posterior_Prob)==0){
  paintor_cs<-rep(0,nrow(paintor_res))
}else{
  paintor_cs_index<-find_index_until_threshold(paintor_res$Posterior_Prob/sum(paintor_res$Posterior_Prob),0.95)
  paintor_cs<-rep(0,nrow(paintor_res))
  paintor_cs[paintor_cs_index]<-1
}

paintor_res_dataframe<-data.frame(SNP = paintor_res$SNP,Paintor_PIP_Either = paintor_res$Posterior_Prob,Paintor_cs = paintor_cs)

mfd_result_all <- mfd_result_all %>%
  left_join(paintor_res_dataframe, by = "SNP")

mfd_result_all <- mfd_result_all%>% replace(is.na(.), 0)

#################Run XMAP########################

library(XMAP)

load(paste0(data_dir,trait,"_xmap_omega.Rdata"))
z_eur <- g1_sub$Z
z_afr <- g2_sub$Z




xmap <- XMAP(simplify2array(list(ld1_sub, ld2_sub)), cbind(z_eur, z_afr), n=c(median(g1_sub$N), median(g2_sub$N)),
             K = 10, Omega = OmegaHat, Sig_E = c(c2,c1), tol = 1e-6,           #c1 is for afr, c2 is for eur
             maxIter = 300, estimate_residual_variance = F, estimate_prior_variance = T,
             estimate_background_variance = F)

#xmap <- XMAP(simplify2array(list(ld1_sub, ld2_sub)), cbind(z_eur, z_afr), n=c(median(g1_sub$N), median(g2_sub$N)),
#             K = 10, tol = 1e-6,           #c1 is for afr, c2 is for eur
#             maxIter = 300, estimate_residual_variance = F, estimate_prior_variance = T,
#             estimate_background_variance = F)


cs1 <- get_CS(xmap, Xcorr = ld1_sub, coverage = 0.95, min_abs_corr = 0.5)
cs2 <- get_CS(xmap, Xcorr = ld2_sub, coverage = 0.95, min_abs_corr = 0.5)
cs_xmap <- unique(unlist(cs1$cs[intersect(names(cs1$cs), names(cs2$cs))],use.names = FALSE))
#cs_xmap
pip_xmap <- get_pip(xmap$gamma)
#plot_CS(pip_xmap, cs_xmap, main = "XMAP")

xmap_res_dataframe<-data.frame(SNP = g1_sub$SNP,XMAP_PIP_Either = pip_xmap, XMAP_cs = 0)
xmap_res_dataframe$XMAP_cs[cs_xmap] <- 1

mfd_result_all <- mfd_result_all %>%
  left_join(xmap_res_dataframe, by = "SNP") %>% replace(is.na(.), 0)


#################Run CARMA-X########################

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
# EUR ancestry
g <- 1
z.list[[g]] <- as.matrix(g1_sub$`Z`)
lambda.list[[g]] <- lambda.power

# AFR ancestry
g <- 2
z.list[[g]] <- as.matrix(g2_sub$`Z`)
lambda.list[[g]] <- lambda.power


ld.list[[1]] <- as.matrix(ld1_sub)
ld.list[[2]] <- as.matrix(ld2_sub)
ld.list[[3]] <- as.matrix(bdiag(ld1_sub, ld2_sub))


# Define the prior distribution
prior.name <- 'Spike-slab'

# Set a random seed
set.seed(123)
# Run CARMA-X

carmax_results <- CARMAX(z.list, ld.list, lambda.list = lambda.list, effect.size.prior = prior.name, LD.estimation = FALSE, outlier.switch = FALSE)

carmax_cs <- tryCatch({
  cs_list_1 <- carmax_results[[1]][["Credible set"]][[2]]
  cs_list_2 <- carmax_results[[2]][["Credible set"]][[2]]
  cs_1 <- if (is.list(cs_list_1) && length(cs_list_1) > 0) sort(unique(unlist(cs_list_1))) else integer(0)
  cs_2 <- if (is.list(cs_list_2) && length(cs_list_2) > 0) sort(unique(unlist(cs_list_2))) else integer(0)
  sort(intersect(cs_1, cs_2))
}, error = function(e) integer(0))
carmax_pip_eur <- carmax_results[[1]]$PIPs
carmax_pip_afr <- carmax_results[[2]]$PIPs
carmax_pip_either <- pmax(carmax_pip_eur, carmax_pip_afr, na.rm = TRUE)
carmax_pip_shared <- pmin(carmax_pip_eur, carmax_pip_afr, na.rm = TRUE)
carmax_res_dataframe<-data.frame(SNP = g1_sub$SNP,CARMAX_PIP_Either = carmax_pip_either, CARMAX_PIP_WB =carmax_pip_eur ,CARMAX_PIP_BB = carmax_pip_afr ,CARMAX_PIP_Shared =carmax_pip_shared, CARMAX_cs = 0)
carmax_res_dataframe$CARMAX_cs[carmax_cs] <- 1
mfd_result_all <- mfd_result_all %>%
  left_join(carmax_res_dataframe, by = "SNP") %>% replace(is.na(.), 0)


mf_output_name<-paste0(result_dir,"MF_result_LOCI_",current_block_id,".RData")
save(mfd_result_all,file =mf_output_name )
  

