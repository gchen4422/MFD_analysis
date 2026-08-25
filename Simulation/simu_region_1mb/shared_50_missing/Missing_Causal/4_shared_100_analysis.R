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
 
if(External_index ==1 ){
	num_snp_missing = round(nrow(zfile)*0.3)
}else if(External_index ==2){
	num_snp_missing = round(nrow(zfile)*0.3+nrow(zfile)*(1-0.3)*0.3)
}


non_causal_SNP = zfile%>%filter(Signal==0)%>%pull(RSID)
EU_missing_SNP = c(sample(non_causal_SNP,num_snp_missing-1),zfile%>%filter(Signal!=0)%>%pull(RSID))
BB_missing_SNP = EU_missing_SNP
zfile<-zfile%>%mutate(EU_missing = ifelse(RSID%in%EU_missing_SNP,1,0),BB_missing = ifelse(RSID%in%BB_missing_SNP,1,0))
 

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
R_mat_list=list("EUR" = EU_cov_subset,"AFR" = BB_cov_subset)

summary_stat_1 = data.frame("SNP" = zfile_subset$RSID, "Beta"=zfile_subset$zscore_1/sqrt(num_EU),"Se"=1/sqrt(num_EU), "Z" =zfile_subset$zscore_1 ,  "N" =num_EU )
summary_stat_2 = data.frame("SNP" = zfile_subset$RSID, "Beta"=zfile_subset$zscore_2/sqrt(num_BB),"Se"=1/sqrt(num_BB), "Z" =zfile_subset$zscore_2 ,  "N" =num_BB )
summary_stat_sd_list = list("EUR" = summary_stat_1,"AFR"=summary_stat_2 )  



summary_stat_1_bim = read.table(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur/loci_",LD_BLOCK,".bim")) %>% filter(V2 %in% zfile_subset$RSID)
summary_stat_2_bim = read.table(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_",LD_BLOCK,".bim")) %>% filter(V2 %in% zfile_subset$RSID)
summary_stat_2_bim$V2 = summary_stat_1_bim$V2


summary_stat_1_allele = data.frame("A1" = summary_stat_1_bim$V5, "A2" = summary_stat_1_bim$V6) 
summary_stat_2_allele = data.frame("A1" = summary_stat_2_bim$V5, "A2" = summary_stat_2_bim$V6)

summary_stat_1_susiex = cbind(zfile_subset[,c(1,2,3)],summary_stat_1_allele,summary_stat_1[,c(2:5)])
summary_stat_2_susiex = cbind(zfile_subset[,c(1,2,3)],summary_stat_2_allele,summary_stat_2[,c(2:5)])

summary_stat_1_susiex$P = 2*pnorm(abs(summary_stat_1_susiex$Z),lower.tail = F)
summary_stat_2_susiex$P = 2*pnorm(abs(summary_stat_2_susiex$Z),lower.tail = F)


write.table(summary_stat_1_susiex,file = paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_eur"),col.names = T,row.names = F,quote=F, sep ="\t")
write.table(summary_stat_2_susiex,file = paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_afr"),col.names = T,row.names = F,quote=F, sep ="\t")
write.table(summary_stat_2_bim,file = paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_missing_causal_",LD_BLOCK,".bim"),col.names = F,row.names = F,quote=F, sep ="\t")



##############################################
#
#		Run MESuSiE Input
#
##############################################


library(MESuSiE)
bayes_fac=3
ancestry_weight=c(bayes_fac/(bayes_fac*2+1),bayes_fac/(bayes_fac*2+1),1/(bayes_fac*2+1))

start.time<-Sys.time()
MESuSiE_res<-meSuSie_core(R_mat_list,summary_stat_sd_list,L=10,residual_variance=NULL,prior_weights=NULL,ancestry_weight=ancestry_weight,optim_method ="optim",estimate_residual_variance =T,max_iter =100)
MESuSiE_cs<-rep(0,length(MESuSiE_res$pip))
MESuSiE_cs[unlist(MESuSiE_res$cs$cs)]<-1
MESuSiE_res_dataframe<-data.frame(RSID = zfile_subset$RSID,MESuSiE_PIP_Either = MESuSiE_res$pip,MESuSiE_PIP_WB =MESuSiE_res$pip_config[,1] ,MESuSiE_PIP_BB =MESuSiE_res$pip_config[,2] ,MESuSiE_PIP_Shared =MESuSiE_res$pip_config[,3],MESuSiE_cs = MESuSiE_cs )

zfile_MESuSiE <- zfile %>%
  left_join(MESuSiE_res_dataframe, by = "RSID")


#################Run SuSiEx########################

#loci_info = fread("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/simulation_locus_info.txt")

chr = median(zfile_subset$CHR)
bp_start = min(as.numeric(zfile_subset$POS))
bp_end = max(as.numeric(zfile_subset$POS))


# Construct the sst_file argument
sst_file <- paste0(data_dir,
                   "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_eur,",data_dir, 
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
susiex.time<-Sys.time()-start.time
susiex.time<-as.numeric(susiex.time,units = "mins")


##############################################
#
#		Run SuSiE
#
##############################################



library(susieR)
if(causal_index==1){
  susie_EU_zfile<-zfile%>%filter(EU_missing==0)
  susie_BB_zfile<-zfile%>%filter(BB_missing==0)
  missing_EU_index = which(zfile$EU_missing!=0)
  missing_BB_index = which(zfile$BB_missing!=0)
  susie_EU_cov<-EU_cov[-missing_EU_index,-missing_EU_index]
  susie_BB_cov<-BB_cov[-missing_BB_index,-missing_BB_index]
}else if(causal_index==2){
  susie_EU_zfile<-zfile%>%filter(EU_missing==0)
  susie_BB_zfile<-zfile
  missing_EU_index = which(zfile$EU_missing!=0)
  susie_EU_cov<-EU_cov[-missing_EU_index,-missing_EU_index]
  susie_BB_cov<-BB_cov
}

susie_EU<-susie_rss(susie_EU_zfile$zscore_1,susie_EU_cov,check_prior=F)
susie_BB<-susie_rss(susie_BB_zfile$zscore_2,susie_BB_cov,check_prior=F)

susie_BB_cs<-rep(0,length(susie_BB$pip))
susie_BB_cs[unlist(susie_BB$sets$cs)]<-1
susie_EU_cs<-rep(0,length(susie_EU$pip))
susie_EU_cs[unlist(susie_EU$sets$cs)]<-1

SuSiE_EU_res_dataframe<-data.frame(RSID = susie_EU_zfile$RSID,SuSiE_PIP_EU = susie_EU$pip,susie_EU_cs = susie_EU_cs )
SuSiE_BB_res_dataframe<-data.frame(RSID = susie_BB_zfile$RSID,SuSiE_PIP_BB = susie_BB$pip,susie_BB_cs = susie_BB_cs )


zfile_MESuSiE_SuSiE <- zfile_MESuSiE %>%
  left_join(SuSiE_EU_res_dataframe, by = "RSID")%>%
  left_join(SuSiE_BB_res_dataframe, by = "RSID")

zfile_MESuSiE_SuSiE <- zfile_MESuSiE_SuSiE%>% replace(is.na(.), 0)
zfile_MESuSiE_SuSiE <- zfile_MESuSiE_SuSiE%>%mutate(SuSiE_cs = ifelse(susie_EU_cs+susie_BB_cs==0,0,1),SuSiE_PIP_Either = pmax(SuSiE_PIP_EU,SuSiE_PIP_BB),SuSiE_PIP_Shared = pmin(SuSiE_PIP_EU,SuSiE_PIP_BB))


##############################################
#
#		Run Paintor
#
##############################################

ld_EU_paintor_name<-paste0(paintor_dir,"Missing_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")
ld_BB_paintor_name<-paste0(paintor_dir,"Missing_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")
fwrite(EU_cov_subset,ld_EU_paintor_name,sep=" ",col.names=F)
fwrite(BB_cov_subset,ld_BB_paintor_name,sep=" ",col.names=F)

 z_file_out<-paste0(paintor_dir,"Missing_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num)
 write.table(zfile_subset,z_file_out,col.names = T,row.names = F,quote=F,sep=" ")  
  
	annotation_file<-matrix(rep(1,nrow(zfile_subset)),ncol=1)
    colnames(annotation_file)<-"coding"
    annotation_paintor_name<-paste0(paintor_dir,"Missing_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".annotations")
    write.table(annotation_file,annotation_paintor_name,col.names = T,row.names = F,quote=F)
    
    input_file_name<-paste0(paintor_dir,"Missing_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".input")
    write.table(paste0("Missing_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),input_file_name,col.names = F,row.names = F,quote=F)


paintor_out_name<-paste0(paintor_dir,"Missing_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num)

paintor_suffix<-paste0("Missing_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num)
input_dir<-paintor_dir
start.time<-Sys.time()
system(paste0("~/gwas_software/PAINTOR_V3.0/PAINTOR -input ",input_file_name," -Zhead  zscore_1,zscore_2 -LDname LD1,LD2 -in ",input_dir," -out ",result_dir," -mcmc "," -annotations coding -Gname ",paintor_suffix," -Lname ",paste0(paintor_suffix,"_bayes_factor")," -RESname ",paste0("mcmc.paintor")))
time.Paintor<-Sys.time()-start.time
time.Paintor<-as.numeric(time.Paintor,units = "mins")
print(paste0("Paintor time is ",time.Paintor))



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


paintor_res_dataframe<-data.frame(RSID = paintor_res$RSID,Paintor_PIP_Either = paintor_res$Posterior_Prob,Paintor_cs = paintor_cs)

zfile_MESuSiE_SuSiE_Paintor <- zfile_MESuSiE_SuSiE %>%
  left_join(paintor_res_dataframe, by = "RSID")

zfile_MESuSiE_SuSiE_Paintor <- zfile_MESuSiE_SuSiE_Paintor%>% replace(is.na(.), 0)
MESuSiE_name<-paste0(result_dir,"MESuSiE_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".RData")
save(zfile_MESuSiE_SuSiE_Paintor,file =MESuSiE_name )

#################Run XMAP########################

library(XMAP)
z_eur <- zfile_subset$zscore_1
z_afr <- zfile_subset$zscore_2

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


ld.list[[1]] <- as.matrix(EU_cov_subset)
ld.list[[2]] <- as.matrix(BB_cov_subset)
ld.list[[3]] <- as.matrix(bdiag(EU_cov_subset, BB_cov_subset))


# Define the prior distribution
prior.name <- 'Spike-slab'

# Set a random seed
set.seed(123)
# Run CARMA-X
start.time<-Sys.time()

carmax_results <- CARMAX(z.list, ld.list, lambda.list = lambda.list, effect.size.prior = prior.name, LD.estimation = FALSE, outlier.switch = FALSE)

carmax.time<-Sys.time()-start.time
carmax.time<-as.numeric(carmax.time,units = "mins")



add_result_name<-paste0(result_dir,"MESuSiE_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_additional.RData")

save(xmap,carmax_results, file = add_result_name)


#system(paste0("rm ",ld_EU_paintor_name))
#system(paste0("rm ",ld_BB_paintor_name))
#system(paste0("rm ",ld_EU_name))
#system(paste0("rm ",ld_BB_name))
