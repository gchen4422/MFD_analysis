args <- as.numeric(commandArgs(TRUE))
if(length(args) < 2) stop("Please provide Block ID and Trait Index arguments.")
current_block_id <- args[1]
trait_idx <- args[2]
setwd("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision")

# Preload GSL dependencies
#gsl_lib <- "/apps/spack/negishi/apps/gsl/2.7.1-gcc-12.2.0-oxu7x5x/lib"
#dyn.load(file.path(gsl_lib, "libgslcblas.so.0.0.0"))
#dyn.load(file.path(gsl_lib, "libgsl.so.27.0.0"))




library(dplyr)
library(data.table)
library(MESuSiE)
library(susieR)
library(MFD)

traits_list <- c("BMI", "DBP", "SBP")
trait <- traits_list[trait_idx] # Select trait based on argument


# Directories
wrk_dir <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/results_", trait, "/")
data_dir <- paste0(wrk_dir, "summary_data/")
result_dir <- paste0(wrk_dir, "result/")
out_dir <- paste0(wrk_dir, "out/")

# Read Candidate Regions
candidate_region_file <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_risk_loci")
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
eas_sum_file <- paste0(data_dir, "eas_loci_", current_block_id, "_summ")

eur_ld_file <- paste0(data_dir, "eur_loci_", current_block_id, ".ld")
afr_ld_file <- paste0(data_dir, "afr_loci_", current_block_id, ".ld")
eas_ld_file <- paste0(data_dir, "eas_loci_", current_block_id, ".ld")


# Check if files exist
if(!file.exists(eur_sum_file) | !file.exists(afr_sum_file) | !file.exists(eas_sum_file) ) {
  stop("Summary statistics file missing")
}

# Load Summary Stats
summary_stat_eur <- fread(eur_sum_file)
colnames(summary_stat_eur)[which(colnames(summary_stat_eur) == "BETA")] <- "Beta"
colnames(summary_stat_eur)[which(colnames(summary_stat_eur) == "SE")] <- "Se"

summary_stat_afr <- fread(afr_sum_file)
colnames(summary_stat_afr)[which(colnames(summary_stat_afr) == "BETA")] <- "Beta"
colnames(summary_stat_afr)[which(colnames(summary_stat_afr) == "SE")] <- "Se"

summary_stat_eas <- fread(eas_sum_file)
setnames(summary_stat_eas, old = "BETA", new = "Beta", skip_absent = TRUE)
setnames(summary_stat_eas, old = "SE",   new = "Se",   skip_absent = TRUE)
summary_stat_eas[, INFO := 1L]
afr_cols <- names(summary_stat_afr)
missing_cols <- setdiff(afr_cols, names(summary_stat_eas))
if (length(missing_cols) > 0) {
  stop("Missing columns in EAS file: ", paste(missing_cols, collapse = ", "))
}
summary_stat_eas <- summary_stat_eas[, ..afr_cols]


# Load LD Matrices
EU_cov <- as.matrix(fread(eur_ld_file))
BB_cov <- as.matrix(fread(afr_ld_file))
EA_cov <- as.matrix(fread(eas_ld_file))

# Ensure Matrix Dimensions match Summary Stats (Basic error checking)
# SNP order is established by the upstream preparation step
rownames(EU_cov) <- summary_stat_eur$SNP
colnames(EU_cov) <- summary_stat_eur$SNP
rownames(BB_cov) <- summary_stat_afr$SNP
colnames(BB_cov) <- summary_stat_afr$SNP
rownames(EA_cov) <- summary_stat_eas$SNP
colnames(EA_cov) <- summary_stat_eas$SNP


##############################################
#
#		Run MF Decision tree (3-ancestry: EUR, AFR, EAS)
#
##############################################

gwas_list <- list(EUR = summary_stat_eur, AFR = summary_stat_afr, EAS = summary_stat_eas)
ld_list   <- list(EUR = EU_cov, AFR = BB_cov, EAS = EA_cov)

#mfd_result_fm = run_mf_decision_fm(gwas_list, ld_list, pop_names = c("EUR", "AFR", "EAS"))

mfd_result_full = run_mf_decision(gwas_list, ld_list, pop_names = c("EUR", "AFR", "EAS"))

mfd_result_ori = mfd_result_full
mfd_result = mfd_result_full$results
##############################################
#
#		Run SuSiE
#
##############################################
susie_EU<-susie_rss(summary_stat_eur$Z,EU_cov, n = median(summary_stat_eur$N))
susie_BB<-susie_rss(summary_stat_afr$Z,BB_cov, n = median(summary_stat_afr$N))
susie_EA<-susie_rss(summary_stat_eas$Z,EA_cov, n = median(summary_stat_eas$N))


susie_BB_cs<-rep(0,length(susie_BB$pip))
susie_BB_cs[unlist(susie_BB$sets$cs)]<-1
susie_EU_cs<-rep(0,length(susie_EU$pip))
susie_EU_cs[unlist(susie_EU$sets$cs)]<-1
susie_EA_cs<-rep(0,length(susie_EA$pip))
susie_EA_cs[unlist(susie_EA$sets$cs)]<-1

SuSiE_EU_res_dataframe<-data.frame(SNP = summary_stat_eur$SNP,SuSiE_PIP_EU = susie_EU$pip,susie_EU_cs = susie_EU_cs )
SuSiE_BB_res_dataframe<-data.frame(SNP = summary_stat_afr$SNP,SuSiE_PIP_BB = susie_BB$pip,susie_BB_cs = susie_BB_cs )
SuSiE_EA_res_dataframe<-data.frame(SNP = summary_stat_eas$SNP,SuSiE_PIP_EA = susie_EA$pip,susie_EA_cs = susie_EA_cs )


mfd_result_all <- mfd_result %>%
  left_join(SuSiE_EU_res_dataframe, by = "SNP")%>%
  left_join(SuSiE_BB_res_dataframe, by = "SNP")%>%
  left_join(SuSiE_EA_res_dataframe, by = "SNP")

mfd_result_all <- mfd_result_all%>% replace(is.na(.), 0)
mfd_result_all <- mfd_result_all%>%mutate(SuSiE_cs = ifelse(susie_EU_cs+susie_BB_cs+susie_EA_cs==0,0,1),SuSiE_PIP_Either = pmax(SuSiE_PIP_EU,SuSiE_PIP_BB,SuSiE_PIP_EA),SuSiE_PIP_Shared = pmin(SuSiE_PIP_EU,SuSiE_PIP_BB,SuSiE_PIP_EA))




##############################################
#
#		Run MESuSiE
#
##############################################

# MESuSiE requires variants present in both datasets
common_snps_two <- intersect(summary_stat_eur$SNP, summary_stat_afr$SNP)
common_snps <- intersect(common_snps_two, summary_stat_eas$SNP)

g1_sub <- summary_stat_eur[SNP %in% common_snps][order(match(SNP, common_snps))]
g2_sub <- summary_stat_afr[SNP %in% common_snps][order(match(SNP, common_snps))]
g3_sub <- summary_stat_eas[SNP %in% common_snps][order(match(SNP, common_snps))]

# Harmonize EAS alleles to EUR: flip, or drop if inconsistent
same_al <- (g3_sub$ALT == g1_sub$ALT) & (g3_sub$REF == g1_sub$REF)
flip_eas <- (g3_sub$ALT == g1_sub$REF) & (g3_sub$REF == g1_sub$ALT)
drop_eas <- !same_al & !flip_eas

if (any(drop_eas)) {
  cat(sprintf("  Dropping %d SNPs with inconsistent alleles\n", sum(drop_eas)))
  keep <- !drop_eas
  common_snps <- common_snps[keep]
  g1_sub <- g1_sub[keep]
  g2_sub <- g2_sub[keep]
  g3_sub <- g3_sub[keep]
  flip_eas <- flip_eas[keep]
}

if (any(flip_eas)) {
  g3_sub[flip_eas, Z    := -Z]
  g3_sub[flip_eas, Beta := -Beta]
  g3_sub[flip_eas, MAF  := 1 - MAF]
  old_alt <- g3_sub$ALT[flip_eas]
  g3_sub[flip_eas, ALT := REF]
  g3_sub[flip_eas, REF := old_alt]
}
cat(sprintf("  EAS allele harmonization: %d same, %d flipped, %d dropped\n",
            sum(same_al), sum(flip_eas), sum(drop_eas)))

# Subset LD matrices based on SNP IDs (assuming LD has row/col names)
ld1_sub <- EU_cov[common_snps, common_snps]
ld2_sub <- BB_cov[common_snps, common_snps]
ld3_sub <- EA_cov[common_snps, common_snps]

# Flip EAS LD to match harmonized allele coding
flip_idx <- which(flip_eas)
if (length(flip_idx) > 0) {
  ld3_sub[flip_idx, ] <- -ld3_sub[flip_idx, ]
  ld3_sub[, flip_idx] <- -ld3_sub[, flip_idx]
}

summary_stat_sd_list = list(g1_sub, g2_sub,g3_sub)
names(summary_stat_sd_list) <- c("EUR", "AFR","EAS")
R_mat_list = list(ld1_sub, ld2_sub,ld3_sub)
names(R_mat_list) <- c("EUR", "AFR","EAS")

mesusie_res <- meSuSie_core(R_mat_list, summary_stat_sd_list, L = 10)

MESuSiE_cs<-rep(0,length(mesusie_res$pip))
MESuSiE_cs[unlist(mesusie_res$cs$cs)]<-1
MESuSiE_res_dataframe<-data.frame(SNP = g1_sub$SNP,MESuSiE_PIP_Either = mesusie_res$pip,MESuSiE_PIP_EUR =mesusie_res$pip_config[,1] ,MESuSiE_PIP_AFR =mesusie_res$pip_config[,2] ,MESuSiE_PIP_EAS =mesusie_res$pip_config[,3] ,MESuSiE_PIP_Shared =mesusie_res$pip_config[,7],MESuSiE_cs = MESuSiE_cs )

mfd_result_all <- mfd_result_all %>%
  left_join(MESuSiE_res_dataframe, by = "SNP")%>% replace(is.na(.), 0)

##############################################
#
#		Run SuSiEx
#
##############################################

chr = g1_sub$CHR[1]
bp_start = min(g1_sub$POS)
bp_end = max(g1_sub$POS)

output_file <- paste0(data_dir, "susiex_input/eur_loci_", current_block_id, "_summ")
if(!dir.exists(dirname(output_file))) {
  dir.create(dirname(output_file), recursive = TRUE)
}
fwrite(g1_sub, file = output_file, quote = FALSE, row.names = FALSE, col.names = TRUE,sep = "\t")

output_file <- paste0(data_dir, "susiex_input/afr_loci_", current_block_id, "_summ")
fwrite(g2_sub,file = paste0(data_dir,"susiex_input/afr_loci_",current_block_id,"_summ"), quote = FALSE, row.names=FALSE, col.names = T,sep = "\t")
output_file <- paste0(data_dir, "susiex_input/eas_loci_", current_block_id, "_summ")
fwrite(g3_sub,file = paste0(data_dir,"susiex_input/eas_loci_",current_block_id,"_summ"), quote = FALSE, row.names=FALSE, col.names = T,sep = "\t")



# Construct the sst_file argument
sst_file <- paste0(data_dir,"susiex_input/eur_loci_", current_block_id, "_summ,",data_dir, "susiex_input/afr_loci_",current_block_id, "_summ,",data_dir,"susiex_input/eas_loci_",current_block_id, "_summ")


susiex_ref_eur <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/risk_loci_ld_eur_",trait)
susiex_ref_afr <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/risk_loci_ld_afr_",trait)
susiex_ref_eas <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/risk_loci_ld_eas_",trait)


# Construct the ref_file argument
ref_file <- paste0(susiex_ref_eur,"/loci_", current_block_id,
                   ",",susiex_ref_afr,"/loci_", current_block_id,
                   ",",susiex_ref_eas,"/loci_", current_block_id)

susiex_out<- paste0(result_dir, "susiex_result/")
if(!dir.exists(susiex_out)) {
  dir.create(susiex_out, recursive = TRUE)
}

# Construct the ld_file argument
ld_file <- paste0(susiex_out, "/eur_loci_", current_block_id,
                  ".ld,",susiex_out,"/afr_loci_", current_block_id, ".ld,",susiex_out,"/eas_loci_", current_block_id, ".ld")

# Define the out_name variable
out_name <- paste0("SuSiEx_",current_block_id,"_output_cs95")

n_gwas_chr = paste(median(g1_sub$N),median(g2_sub$N),median(g3_sub$N),sep = ",")

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
  "--chr_col=2,2,2",
  "--bp_col=3,3,3",
  "--snp_col=1,1,1",
  "--a1_col=5,5,5",
  "--a2_col=6,6,6",
  "--eff_col=9,9,9",
  "--se_col=10,10,10",
  "--pval_col=4,4,4",
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
sumstats = cbind(g1_sub[,c(1,8)],g2_sub[,c(8)],g3_sub[,c(8)])
colnames(sumstats) = c("SNP","zscore_1","zscore_2","zscore_3")
fwrite(sumstats,paste0(paintor_dir,"loci_",current_block_id),col.names = T,row.names = F,quote=F,sep=" ")


ld_EU_paintor_name<-paste0(paintor_dir,"loci_",current_block_id,".LD1")
ld_BB_paintor_name<-paste0(paintor_dir,"loci_",current_block_id,".LD2")
ld_EA_paintor_name<-paste0(paintor_dir,"loci_",current_block_id,".LD3")
fwrite(ld1_sub,ld_EU_paintor_name,col.names = F,row.names = F,quote=F,sep=" ")
fwrite(ld2_sub,ld_BB_paintor_name,col.names = F,row.names = F,quote=F,sep=" ")
fwrite(ld3_sub,ld_EA_paintor_name,col.names = F,row.names = F,quote=F,sep=" ")


paintor_suffix<-paste0("loci_",current_block_id)
input_dir<-paintor_dir
setwd(input_dir)
annotation_file<-matrix(rep(1,nrow(g1_sub)),ncol=1)
colnames(annotation_file)<-"coding"
annotation_paintor_name<-paste0(paintor_dir,"loci_",current_block_id,".annotations")
write.table(annotation_file,annotation_paintor_name,col.names = T,row.names = F,quote=F)
system(paste0("~/gwas_software/PAINTOR_V3.0/PAINTOR -input ",input_file_name," -Zhead  zscore_1,zscore_2,zscore_3 -LDname LD1,LD2,LD3 -in ",input_dir," -out ",result_dir," -mcmc "," -annotations coding -Gname ", paintor_suffix," -Lname ",paste0(paintor_suffix,"_bayes_factor")," -RESname ",paste0("mcmc.paintor")))


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

# Load pre-computed Omega matrix and intercepts from LD-score regression
# OmegaHat is ordered AFR(1), EUR(2), EAS(3); c1=AFR, c2=EUR, c3=EAS
load(paste0(data_dir, trait, "_xmap_omega.Rdata"))

# Reorder to match data ordering: EUR(1), AFR(2), EAS(3)
perm <- c(2, 1, 3)
OmegaHat <- OmegaHat[perm, perm]

# Pairwise LD-score regression can yield a non-PD 3x3 Omega; project to nearest PD
OmegaHat <- as.matrix(Matrix::nearPD(OmegaHat, corr = FALSE, keepDiag = FALSE)$mat)

z_eur <- g1_sub$Z
z_afr <- g2_sub$Z
z_eas <- g3_sub$Z

#xmap <- XMAP(simplify2array(list(ld1_sub, ld2_sub, ld3_sub)),
#             cbind(z_eur, z_afr, z_eas),
#             n=c(median(g1_sub$N), median(g2_sub$N), median(g3_sub$N)),
#             K = 10, Omega = OmegaHat, Sig_E = c(c2, c1, c3), tol = 1e-6,
#             maxIter = 300, estimate_residual_variance = F, estimate_prior_variance = T,
#             estimate_background_variance = F)

xmap <- XMAP(simplify2array(list(ld1_sub, ld2_sub, ld3_sub)),
             cbind(z_eur, z_afr, z_eas),
             n=c(median(g1_sub$N), median(g2_sub$N), median(g3_sub$N)),
             K = 10, tol = 1e-6,
             maxIter = 300, estimate_residual_variance = F, estimate_prior_variance = T,
             estimate_background_variance = F)

cs1 <- get_CS(xmap, Xcorr = ld1_sub, coverage = 0.95, min_abs_corr = 0.5)
cs2 <- get_CS(xmap, Xcorr = ld2_sub, coverage = 0.95, min_abs_corr = 0.5)
cs3 <- get_CS(xmap, Xcorr = ld3_sub, coverage = 0.95, min_abs_corr = 0.5)
# Joint CS: signals that pass purity in at least 2 of K ancestries
# (strict all-K intersection is too conservative for K=3 due to LD diversity)
cs_list <- list(cs1$cs, cs2$cs, cs3$cs)
cs_names <- lapply(cs_list, function(x) if (is.null(x)) character(0) else names(x))
all_cs_names <- unique(unlist(cs_names))
cs_name_counts <- table(unlist(cs_names))
cs_shared_names <- names(cs_name_counts[cs_name_counts >= 2])
# Get SNP indices from the first ancestry that has each shared CS
cs_xmap <- unique(unlist(lapply(cs_shared_names, function(nm) {
  for (cs in cs_list) { if (!is.null(cs) && nm %in% names(cs)) return(cs[[nm]]) }
  return(NULL)
}), use.names = FALSE))
pip_xmap <- get_pip(xmap$gamma)

xmap_res_dataframe<-data.frame(SNP = g1_sub$SNP,XMAP_PIP_Either = pip_xmap, XMAP_cs = 0)
xmap_res_dataframe$XMAP_cs[cs_xmap] <- 1

mfd_result_all <- mfd_result_all %>%
  left_join(xmap_res_dataframe, by = "SNP") %>% replace(is.na(.), 0)


#################CARMA-X########################
# CARMA-X is excluded because the benchmarked implementation supports two ancestries.


mf_output_name<-paste0(result_dir,"MF_result_LOCI_",current_block_id,".RData")
save(mfd_result_all, mfd_result_ori, file =mf_output_name )
