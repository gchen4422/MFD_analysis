args<-as.numeric(commandArgs(TRUE))
LD_BLOCK = args[1]
h2_num = args[2]
num_causal = args[3]

library(data.table)

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_admix/causal_num_",num_causal,"/")
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

# Admix data is RDS (no .bim), so use dummy alleles (A/T) consistently
summary_stat_1_allele = data.frame("A1" = rep("A", nrow(zfile)), "A2" = rep("T", nrow(zfile)))
summary_stat_2_allele = summary_stat_1_allele

summary_stat_1_susiex = cbind(zfile[,c(1,2,3)],summary_stat_1_allele,summary_stat_1[,c(2:5)])
summary_stat_2_susiex = cbind(zfile[,c(1,2,3)],summary_stat_2_allele,summary_stat_2[,c(2:5)])

summary_stat_1_susiex$P = 2*pnorm(abs(summary_stat_1_susiex$Z),lower.tail = F)
summary_stat_2_susiex$P = 2*pnorm(abs(summary_stat_2_susiex$Z),lower.tail = F)

write.table(summary_stat_1_susiex,file = paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_eur"),col.names = T,row.names = F,quote=F, sep ="\t")
write.table(summary_stat_2_susiex,file = paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_afr"),col.names = T,row.names = F,quote=F, sep ="\t")



N_eur <- zfile$N_1[1]
N_afr <- zfile$N_2[1]

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
MESuSiE_res<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10,residual_variance=NULL,prior_weights=NULL,ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =T,max_iter =100)
MESuSiE.time<-Sys.time()-start.time
MESuSiE.time<-as.numeric(MESuSiE.time,units = "mins")
MESuSiE_name<-paste0(result_dir,"MESuSiE_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".RData")

#################Run SuSiEx########################
# Write SuSiEx-format LD files from the LD matrices (once per region)
susiex_ld_eur_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_admix/susiex_ld_eur"
susiex_ld_afr_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_admix/susiex_ld_afr"
system(paste0("mkdir -p ", susiex_ld_eur_dir))
system(paste0("mkdir -p ", susiex_ld_afr_dir))

ld_prefix_eur <- file.path(susiex_ld_eur_dir, paste0("region_", LD_BLOCK, ".ld"))
ld_prefix_afr <- file.path(susiex_ld_afr_dir, paste0("region_", LD_BLOCK, ".ld"))

# Only create if not already present (shared across causal/h2 combos)
if (!file.exists(paste0(ld_prefix_eur, ".ld.bin"))) {
  n_snp <- nrow(EU_cov)
  # Write LD binary (float32, p*p matrix)
  con <- file(paste0(ld_prefix_eur, ".ld.bin"), "wb")
  writeBin(as.vector(EU_cov), con, size = 4)
  close(con)
  con <- file(paste0(ld_prefix_afr, ".ld.bin"), "wb")
  writeBin(as.vector(BB_cov), con, size = 4)
  close(con)
  # Write _ref.bim (CHR SNP CM POS A1 A2) with dummy alleles
  snp_parts <- strsplit(as.character(zfile$RSID), ":")
  bim_df <- data.frame(
    CHR = zfile$CHR,
    SNP = zfile$RSID,
    CM  = 0,
    POS = zfile$POS,
    A1  = "A",
    A2  = "T"
  )
  write.table(bim_df, paste0(ld_prefix_eur, "_ref.bim"), col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(bim_df, paste0(ld_prefix_afr, "_ref.bim"), col.names = F, row.names = F, quote = F, sep = "\t")
  # Compute real MAFs from genotype data
  geno_dir_eur <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_admix/real_data_eur_100regions"
  geno_dir_afr <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_admix/real_data_afr_100regions"
  EUR_raw <- readRDS(file.path(geno_dir_eur, paste0("region", LD_BLOCK, ".rds")))
  AFR_raw <- readRDS(file.path(geno_dir_afr, paste0("region", LD_BLOCK, ".rds")))
  EU_geno <- t(EUR_raw)
  BB_geno <- t(AFR_raw)
  colnames(EU_geno) <- rownames(EUR_raw)
  colnames(BB_geno) <- rownames(AFR_raw)
  EU_geno <- EU_geno[, as.character(zfile$RSID), drop = FALSE]
  BB_geno <- BB_geno[, as.character(zfile$RSID), drop = FALSE]
  freq_eur <- colMeans(EU_geno) / 2
  freq_afr <- colMeans(BB_geno) / 2
  maf_eur <- pmin(freq_eur, 1 - freq_eur)
  maf_afr <- pmin(freq_afr, 1 - freq_afr)

  # Write _maf.txt for MultiSuSiE (one MAF per line, no header)
  write.table(maf_eur, paste0(susiex_ld_eur_dir, "/region_", LD_BLOCK, "_maf.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)
  write.table(maf_afr, paste0(susiex_ld_afr_dir, "/region_", LD_BLOCK, "_maf.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)

  # Write _frq.frq with real MAFs (for SuSiEx)
  frq_df_eur <- data.frame(
    CHR = zfile$CHR,
    SNP = zfile$RSID,
    A1  = "A",
    A2  = "T",
    MAF = maf_eur,
    NCHROBS = 2 * N_eur
  )
  frq_df_afr <- data.frame(
    CHR = zfile$CHR,
    SNP = zfile$RSID,
    A1  = "A",
    A2  = "T",
    MAF = maf_afr,
    NCHROBS = 2 * N_afr
  )
  frq_header <- " CHR         SNP   A1   A2          MAF  NCHROBS"
  writeLines(frq_header, paste0(ld_prefix_eur, "_frq.frq"))
  write.table(frq_df_eur, paste0(ld_prefix_eur, "_frq.frq"), col.names = F, row.names = F, quote = F, sep = "\t", append = TRUE)
  writeLines(frq_header, paste0(ld_prefix_afr, "_frq.frq"))
  write.table(frq_df_afr, paste0(ld_prefix_afr, "_frq.frq"), col.names = F, row.names = F, quote = F, sep = "\t", append = TRUE)
}

loci_info = fread("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_admix/simulation_locus_info.txt")

chr = loci_info$CHR[LD_BLOCK]
bp_start = loci_info$MinBP[LD_BLOCK]
bp_end = loci_info$MaxBP[LD_BLOCK]

sst_file <- paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_eur,",
                   data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_afr")

# No PLINK ref files needed — SuSiEx uses pre-computed LD directly
# ref_file still needed for SNP matching; use the _ref.bim prefix
ref_file <- paste0(ld_prefix_eur, "_ref,", ld_prefix_afr, "_ref")

ld_file <- paste0(ld_prefix_eur, ",", ld_prefix_afr)

out_name <- paste0("SuSiEx_CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_output_cs95")

command <- paste(
  "~/gwas_software/SuSiEx/bin/SuSiEx",
  paste0("--sst_file=", sst_file),
  paste0("--ref_file=", ref_file),
  paste0("--ld_file=", ld_file),
  paste0("--n_gwas=", N_eur, ",", N_afr),
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
susie_EU<-susie_rss(zfile$zscore_1,EU_cov,check_prior=F, n = N_eur)
susie_BB<-susie_rss(zfile$zscore_2,BB_cov,check_prior=F, n = N_afr)

susie.time<-Sys.time()-start.time
susie.time<-as.numeric(susie.time,units = "mins")

#kriging_rss(zfile$zscore_1,EU_cov, n = 300000)
#kriging_rss(zfile$zscore_1,EU_cov, n = 300000)


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
# Create PAINTOR .input file (not generated by script 3 for admix)
write.table(paintor_suffix, input_file_name, col.names = F, row.names = F, quote = F)
# Create PAINTOR .annotations file (single "coding" column of all 1s, same as shared_50)
annotation_file <- matrix(rep(1, nrow(zfile)), ncol = 1)
colnames(annotation_file) <- "coding"
write.table(annotation_file, paste0(paintor_dir, paintor_suffix, ".annotations"), col.names = T, row.names = F, quote = F)
start.time<-Sys.time()
system(paste0("~/gwas_software/PAINTOR_V3.0/PAINTOR -input ",input_file_name," -Zhead  zscore_1,zscore_2 -LDname LD1,LD2 -in ",input_dir," -out ",result_dir," -mcmc "," -annotations coding -Gname ",paintor_suffix," -Lname ",paste0(paintor_suffix,"_bayes_factor")," -RESname ",paste0("mcmc.paintor")))
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

xmap <-   XMAP(simplify2array(list(EU_cov, BB_cov)), cbind(z_eur, z_afr), n = c(N_eur, N_afr),
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
# EUR ancestry
g <- 1
z.list[[g]] <- as.matrix(summary_stat_1$`Z`)
lambda.list[[g]] <- lambda.power

# AFR ancestry
g <- 2
z.list[[g]] <- as.matrix(summary_stat_2$`Z`)
lambda.list[[g]] <- lambda.power


ld.list[[1]] <- as.matrix(EU_cov)
ld.list[[2]] <- as.matrix(BB_cov)


# Define the prior distribution
prior.name <- 'Spike-slab'

# Set a random seed
set.seed(123)
# Run CARMA-X
start.time<-Sys.time()

carmax_results <- CARMAX(z.list, ld.list, lambda.list = lambda.list, effect.size.prior = prior.name, LD.estimation = TRUE, outlier.switch = FALSE)

carmax.time<-Sys.time()-start.time
carmax.time<-as.numeric(carmax.time,units = "mins")




save(susie_EU,susie_BB,MESuSiE_res,xmap,carmax_results, file = MESuSiE_name)
#################Run PIPSORT########################
# library(dplyr)
#
# z_file_eur = summary_stat_1_susiex %>% dplyr::select(RSID,Z)
# z_file_afr = summary_stat_2_susiex %>% dplyr::select(RSID,Z)
#
#
# write.table(z_file_eur,file = paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_z_eur"),col.names = F,row.names = F,quote=F, sep ="\t")
# write.table(z_file_afr,file = paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_z_afr"),col.names = F,row.names = F,quote=F, sep ="\t")
#
# line1 <- paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_z_eur")
# line2 <- paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_z_afr")
# newfile <- data.frame(z_file = c(line1, line2), stringsAsFactors = FALSE)
# out_path <- paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_zfiles.txt")
# write.table(newfile, file = out_path, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
#
#
# line1_ld <- ld_EU_paintor_name
# line2_ld <- ld_BB_paintor_name
# newfile_ld <- data.frame(z_file = c(line1_ld, line2_ld), stringsAsFactors = FALSE)
# out_path_ld <- paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_ldfiles.txt")
# write.table(newfile_ld, file = out_path_ld, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
#
#
#
#
# snp_map = summary_stat_1_susiex %>% dplyr::select(RSID) %>% mutate(SNP1 = 0:(nrow(summary_stat_1_susiex)-1), SNP2 =0:(nrow(summary_stat_1_susiex)-1) )
# snp_map_name = paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_snp_map")
# write.table(snp_map,file =snp_map_name ,col.names = F,row.names = F,quote=F, sep =",")
#
#
# pipsort_outname = paste0(result_dir,paintor_suffix,"_pipsort")
#
#
# #limit the maximum number of causal variants at 1 if either study has more than 1700 variants and at 2 if either study has more than 325 variants
#
# n_max <- max(nrow(summary_stat_1), nrow(summary_stat_2))
# causal_number_pipsort <- if (n_max > 1700) 1 else if (n_max > 325) 2 else causal_number_pipsort
#
#
# cmd <- sprintf(
#   '/home/chen4422/gwas_software/pipsort/PIPSORT -c %d -l %s -z %s -m %s -n 300000,300000 -p 0.25 -o %s',
#   causal_number_pipsort, shQuote(out_path_ld), shQuote(out_path), shQuote(snp_map_name), shQuote(pipsort_outname)
# )
#
#
#
# start.time<-Sys.time()
#
# system(cmd)
#
# pipsort.time<-Sys.time()-start.time
# pipsort.time<-as.numeric(pipsort.time,units = "mins")
#
#
#
#
#
#
# py <- "/apps/spack/negishi/apps/anaconda/2022.10-py39-gcc-8.5.0-sjvibry/bin/python"
# script <- "/home/chen4422/gwas_software/pipsort/utils/get_global_pips.py"
# args <- c(
#   paste0(pipsort_outname, "_study0_post.txt"),
#   paste0(pipsort_outname, "_study1_post.txt"),
#   paste0(pipsort_outname, "_shared_pips.txt"),
#   paste0(pipsort_outname, "_global_pips.txt")
# )
#
# system2(py, args = c(script, args))
#
# script <- "/home/chen4422/gwas_software/pipsort/utils/get_not_shared_pips.py"
# system2(
#   command = py,
#   args = c(
#     script,
#     paste0(pipsort_outname, "_shared_pips.txt"),
#     paste0(pipsort_outname, "_global_pips.txt"),
#     paste0(pipsort_outname, "_not_shared_pips.txt")
#   )
# )
##############################################
#
#		Save results/time used
#
##############################################
time_out_name<-paste0(result_dir,"time_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".txt")
time_vec<-c(MESuSiE.time,susie.time,Paintor.time,xmap.time,susiex.time,carmax.time)
write.table(time_vec,time_out_name,col.names = F,row.names = F,quote=F)
