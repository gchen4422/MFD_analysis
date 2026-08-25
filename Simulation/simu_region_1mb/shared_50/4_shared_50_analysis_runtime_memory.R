args <- as.numeric(commandArgs(TRUE))
LD_BLOCK = args[1]
h2_num = args[2]
num_causal = args[3]

library(data.table)
library(peakRAM) # Required for memory tracking

wrk_dir <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_", num_causal, "/")
system(paste0("mkdir -p ", wrk_dir))
data_dir <- paste0(wrk_dir, "summary_data/")
system(paste0("mkdir -p ", data_dir))
result_dir <- paste0(wrk_dir, "result/")
system(paste0("mkdir -p ", result_dir))
paintor_dir <- data_dir

##############################################
#
#		Format MESuSiE Input
#
##############################################
zfile <- read.table(paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num), header = T)
summary_stat_1 = data.frame("SNP" = zfile$RSID, "Beta" = zfile$zscore_1 / sqrt(zfile$N_1), "Se" = 1 / sqrt(zfile$N_1), "Z" = zfile$zscore_1, "N" = zfile$N_1)
summary_stat_2 = data.frame("SNP" = zfile$RSID, "Beta" = zfile$zscore_2 / sqrt(zfile$N_2), "Se" = 1 / sqrt(zfile$N_2), "Z" = zfile$zscore_2, "N" = zfile$N_2)

summary_stat_1_bim = read.table(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur/loci_", LD_BLOCK, ".bim"))
summary_stat_2_bim = read.table(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_", LD_BLOCK, ".bim"))
summary_stat_2_bim$V2 = summary_stat_1_bim$V2

summary_stat_1_allele = data.frame("A1" = summary_stat_1_bim$V5, "A2" = summary_stat_1_bim$V6)
summary_stat_2_allele = data.frame("A1" = summary_stat_2_bim$V5, "A2" = summary_stat_2_bim$V6)

summary_stat_1_susiex = cbind(zfile[, c(1, 2, 3)], summary_stat_1_allele, summary_stat_1[, c(2:5)])
summary_stat_2_susiex = cbind(zfile[, c(1, 2, 3)], summary_stat_2_allele, summary_stat_2[, c(2:5)])

summary_stat_1_susiex$P = 2 * pnorm(abs(summary_stat_1_susiex$Z), lower.tail = F)
summary_stat_2_susiex$P = 2 * pnorm(abs(summary_stat_2_susiex$Z), lower.tail = F)

write.table(summary_stat_1_susiex, file = paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_eur"), col.names = T, row.names = F, quote = F, sep = "\t")
write.table(summary_stat_2_susiex, file = paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_afr"), col.names = T, row.names = F, quote = F, sep = "\t")
write.table(summary_stat_2_bim, file = paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_", LD_BLOCK, ".bim"), col.names = F, row.names = F, quote = F, sep = "\t")

summary_stat_sd_list = list("EUR" = summary_stat_1, "AFR" = summary_stat_2)

EU_cov <- as.matrix(fread(paste0(paintor_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, ".LD1")))
colnames(EU_cov) <- zfile$RSID
BB_cov <- as.matrix(fread(paste0(paintor_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, ".LD2")))
colnames(BB_cov) <- zfile$RSID
R_mat_list = list("EUR" = EU_cov, "AFR" = BB_cov)

##############################################
#
#		Run MESuSiE (with peakRAM)
#
##############################################

library(MESuSiE)
bayes_fac = 3
ancestry_weight = c(bayes_fac / (bayes_fac * 2 + 1), bayes_fac / (bayes_fac * 2 + 1), 1 / (bayes_fac * 2 + 1))

# Wrap execution in peakRAM
mesusie_perf <- peakRAM(
  MESuSiE_res <- meSuSie_core(R_mat_list, summary_stat_sd_list, L = 10, residual_variance = NULL, prior_weights = NULL, ancestry_weight = ancestry_weight, optim_method = "optim", estimate_residual_variance = T, max_iter = 100)
)
#MESuSiE_name <- paste0(result_dir, "MESuSiE_CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, ".RData")

#################Run SuSiEx (with peakRAM)########################

loci_info = fread("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/simulation_locus_info.txt")

chr = loci_info$CHR[LD_BLOCK]
bp_start = loci_info$MinBP[LD_BLOCK]
bp_end = loci_info$MaxBP[LD_BLOCK]

sst_file <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_", num_causal, 
                   "/summary_data/CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_eur,/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_", 
                   num_causal, "/summary_data/CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_afr")

ref_file <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur/loci_", LD_BLOCK, 
                   ",/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_", LD_BLOCK)

ld_file <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/susiex_ld_eur/loci_", LD_BLOCK, 
                  ".ld,/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/susiex_ld_afr/loci_", LD_BLOCK, ".ld")

out_name <- paste0("SuSiEx_CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_output_cs95")

# Define a temporary log file for memory usage
susiex_log <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/susiex_result/",out_name,"_log.txt")

# Construct the command with /usr/bin/time -v
# -v: verbose (gives memory info)
# -o: outputs stats to a file so R can read it

command <- paste(
  "/usr/bin/time -v -o", susiex_log,
  "~/gwas_software/SuSiEx/bin/SuSiEx",
  paste0("--sst_file=", sst_file),
  paste0("--ref_file=", ref_file),
  paste0("--ld_file=", ld_file),
  "--n_gwas=300000,300000",
  "--out_dir=./susiex_result",
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


# Run system command
start.time<-Sys.time()

system(command)

susiex.time<-Sys.time()-start.time
susiex.time<-as.numeric(susiex.time,units = "mins")

# --- Parse the Memory Log ---
# Record maximum resident set size and elapsed wall-clock time
log_content <- readLines(susiex_log)
max_rss_line <- grep("Maximum resident set size", log_content, value = TRUE)
elapsed_line <- grep("Elapsed \\(wall clock\\) time", log_content, value = TRUE)

# Extract numbers
# Memory is in kbytes, convert to MiB (1024 kbytes = 1 MiB)
susiex_ram_kb <- as.numeric(gsub("[^0-9]", "", max_rss_line))
susiex_ram_mib <- susiex_ram_kb / 1024

# Use R timing for elapsed time and the process log for peak memory.

# Create the manual performance row for SuSiEx
susiex_perf <- data.frame(
  Function_Call = "SuSiEx",
  Elapsed_Time_sec = susiex.time*60, # Elapsed time from the R timing block
  Total_RAM_Used_MiB = NA,
  Peak_RAM_Used_MiB = susiex_ram_mib
)

##############################################
#
#		Run SuSiE (with peakRAM)
#
##############################################
library(susieR)

# Wrap execution in peakRAM
susie_perf <- peakRAM({
  susie_EU <- susie_rss(zfile$zscore_1, EU_cov, check_prior = F, n = 300000)
  susie_BB <- susie_rss(zfile$zscore_2, BB_cov, check_prior = F, n = 300000)
})

##############################################
#
#		Run Paintor (with peakRAM)
#
##############################################
paintor_out_name <- paste0(paintor_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num)
input_file_name <- paste0(paintor_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, ".input")
ld_EU_paintor_name <- paste0(paintor_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, ".LD1")
ld_BB_paintor_name <- paste0(paintor_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, ".LD2")
paintor_suffix <- paste0("CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num)
input_dir <- paintor_dir


paintor_log <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/paintor_runtime_memory/", paintor_suffix, "_log.txt")


# Wrap execution in peakRAM

paintor_cmd <- paste0(
  "/usr/bin/time -v -o ", paintor_log, " ",
  "~/gwas_software/PAINTOR_V3.0/PAINTOR -input ", input_file_name, 
  " -Zhead zscore_1,zscore_2 -LDname LD1,LD2 -in ", input_dir, 
  " -out ", result_dir, " -mcmc ", 
  " -annotations coding -Gname ", paintor_suffix, 
  " -Lname ", paste0(paintor_suffix, "_bayes_factor"), 
  " -RESname ", paste0("mcmc.paintor")
)
start.time<-Sys.time()
system(paintor_cmd)
Paintor.time<-Sys.time()-start.time
Paintor.time<-as.numeric(Paintor.time,units = "mins")

log_content_p <- readLines(paintor_log)
max_rss_line_p <- grep("Maximum resident set size", log_content_p, value = TRUE)
paintor_ram_kb <- as.numeric(gsub("[^0-9]", "", max_rss_line_p))
paintor_ram_mib <- paintor_ram_kb / 1024

paintor_perf <- data.frame(
  Function_Call = "Paintor",
  Elapsed_Time_sec = Paintor.time*60,
  Total_RAM_Used_MiB = NA,
  Peak_RAM_Used_MiB = paintor_ram_mib
)

#################Run XMAP (with peakRAM)########################

library(XMAP)
z_afr <- zfile$zscore_2
z_eur <- zfile$zscore_1

# Wrap execution in peakRAM
xmap_perf <- peakRAM(
  xmap <- XMAP(simplify2array(list(EU_cov, BB_cov)), cbind(z_eur, z_afr), n = c(300000, 300000),
               K = 10, tol = 1e-6,
               maxIter = 300, estimate_residual_variance = F, estimate_prior_variance = T,
               estimate_background_variance = F)
)

#################Run CARMA-X (with peakRAM)########################
gsl_lib <- "/apps/spack/negishi/apps/gsl/2.7.1-gcc-12.2.0-oxu7x5x/lib"
dyn.load(file.path(gsl_lib, "libgslcblas.so.0.0.0"))
dyn.load(file.path(gsl_lib, "libgsl.so.27.0.0"))

library(CARMAX)
library(Matrix)

lambda.power <- 1
z.list <- list()
ld.list <- list()
full.w.list <- list()
lambda.list <- c()
label.list <- list()
input.prior.list <- list()

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

# Wrap execution in peakRAM
carmax_perf <- peakRAM(
  carmax_results <- CARMAX(z.list, ld.list, lambda.list = lambda.list, effect.size.prior = prior.name, LD.estimation = FALSE, outlier.switch = FALSE)
)

# Save result objects
#save(susie_EU, susie_BB, MESuSiE_res, xmap, carmax_results, file = MESuSiE_name)

##############################################
#
#		Save Performance Results
#
##############################################

# Combine all performance dataframes
# peakRAM returns columns: Function_Call, Elapsed_Time_sec, Total_RAM_Used_MiB, Peak_RAM_Used_MiB
final_perf_df <- rbind(mesusie_perf, susie_perf, paintor_perf, xmap_perf, susiex_perf, carmax_perf)

# Add Method Name for clarity
final_perf_df$Method <- c("MESuSiE", "SuSiE", "Paintor", "XMAP", "SuSiEx", "CARMAX")

# Add Time in Minutes
final_perf_df$Time_Min <- final_perf_df$Elapsed_Time_sec / 60

# Reorder columns for readability
final_perf_df <- final_perf_df[, c("Method", "Time_Min", "Elapsed_Time_sec", "Peak_RAM_Used_MiB")]

# Print to console for log
print(final_perf_df)

# Save to file
time_out_name<-paste0(result_dir,"Runtime_memory_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".txt")
write.table(final_perf_df,time_out_name,col.names = T,row.names = F,quote=F)
