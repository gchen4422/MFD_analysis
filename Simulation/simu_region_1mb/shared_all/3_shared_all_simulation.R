# Loading necessary libraries
library(data.table)
library(snpStats)
library(mvtnorm)
library(dplyr)

# Setting constants and initializing parameters
#args <- as.numeric(commandArgs(TRUE))
#LD_BLOCK <- args[1]
num_EU <- 300000
num_BB <- 300000

geno_dir_eur <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur")
LD_dir_eur <- geno_dir_eur

geno_dir_afr <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr")
LD_dir_afr <- geno_dir_afr




# Function to process genotype data
process_geno <- function(geno_data) {
  apply(geno_data, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x <- scale(x - mean(x))
    return(x)
  })
}

#Loci 76 only contains 1 snp
for (LD_BLOCK in 1:100){
  # Seed initialization
  set.seed(1128)
  seed_generator <- sample.int(1000, 100)[LD_BLOCK]
  set.seed(seed_generator)
  
  
  
  # Reading and processing EU data
  outfile_EU <- paste0(geno_dir_eur, "/loci_", LD_BLOCK)
  EU_plink <- read.plink(paste0(outfile_EU, ".bed"))
  #order_CHR_POS <- EU_plink$map %>% mutate(CHR_POS = paste0(chromosome, "_", position)) %>% pull(CHR_POS)
  EU_plink_geno <- as(EU_plink$genotypes, "numeric")
  EU_plink_geno <- process_geno(EU_plink_geno)
  
  ld_EU_name <- paste0(LD_dir_eur, "/loci_", LD_BLOCK, ".ld")
  EU_cov <- as.matrix(fread(ld_EU_name))
  #ld_EU_name_1kg <- paste0(LD_dir_eur, "/loci_1kg_", LD_BLOCK, ".ld")
  #EU_cov_1kg = as.matrix(fread(ld_EU_name_1kg))
  colnames(EU_cov) = colnames(EU_plink_geno)
  
  
  # Reading and processing BB data
  outfile_BB <- paste0(geno_dir_afr, "/loci_", LD_BLOCK)
  BB_plink <- read.plink(paste0(outfile_BB, ".bed"))
  BB_plink$map$CHR_POS = paste0(BB_plink$map$chromosome,":",BB_plink$map$position)
  BB_plink_geno <- as(BB_plink$genotypes, "numeric")
  colnames(BB_plink_geno) = BB_plink$map$CHR_POS
  BB_plink_geno <- process_geno(BB_plink_geno)
  
  ld_BB_name <- paste0(LD_dir_afr, "/loci_", LD_BLOCK, ".ld")
  BB_cov<-as.matrix(fread(ld_BB_name))
  #ld_BB_name_1kg <- paste0(LD_dir_afr, "/loci_1kg_", LD_BLOCK, ".ld")
  #BB_cov_1kg<-as.matrix(fread(ld_BB_name_1kg))
  
  colnames(BB_cov) = BB_plink$map$CHR_POS
  
  
  # Check SNP position mismatch 
  if(any(colnames(EU_cov)!=EU_plink$map$snp.name)|any(colnames(BB_cov)!=BB_plink$map$CHR_POS)){
    cat("SNP mismatch of geno and LD reference")
  }
  # Checking allele match
  if (any(EU_plink$map$allele.1 != BB_plink$map$allele.1)) {
    cat("Allele mismatch of the two ancestries")
  }
  # Initializing true prior probabilities
  true_prior_prob <- rep(1 / nrow(BB_cov), nrow(BB_cov))

  
  for(num_causal in 1:3){
    
    wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/causal_num_",num_causal,"/")
    system(paste0("mkdir -p ",wrk_dir))
    data_dir<-paste0(wrk_dir,"summary_data/")
    result_dir<-paste0(wrk_dir,"result/")
    out_dir<-paste0(wrk_dir,"out/")
    system(paste0("mkdir -p ",data_dir))
    system(paste0("mkdir -p ",result_dir))
    system(paste0("mkdir -p ",out_dir))
    
    ##sample number of ancestry-specific causal SNP, and then sample from SNP list
    ##Better way can be sample causal SNP list first, and sample 0,1 to define as EU/BB specific
    
    num_selected_SNP<-c(1,3,5)[num_causal]
    
    shared_causal_SNP_list<-BB_plink$map$CHR_POS[sample(seq(1,nrow(BB_cov)),num_selected_SNP,replace=F,true_prior_prob)]
    EU_plink_causal<-matrix(EU_plink_geno[,colnames(EU_plink_geno)%in%shared_causal_SNP_list],ncol=num_selected_SNP)
    BB_plink_causal<-matrix(BB_plink_geno[,colnames(BB_plink_geno)%in%shared_causal_SNP_list],ncol=num_selected_SNP)
    
    
    beta_shared<-rmvnorm(num_selected_SNP,mean=rep(0,2),sigma=matrix(c(1,0.8,0.8,1),ncol=2,nrow=2))
    
    for(h2_num in 1:2){
      
      h2= c(1e-4,2e-4)[h2_num] 
      h2_BB<-h2*(num_selected_SNP)
      beta_BB_shared<-c(beta_shared[,2])*sqrt(h2_BB/drop(var(BB_plink_causal%*%c(beta_shared[,2]))))         
      beta_BB_all<-rep(0,nrow(BB_cov))
      beta_BB_all[which(BB_plink$map$CHR_POS%in%c(shared_causal_SNP_list))]<-beta_BB_shared
      beta_BB_marginal<-as.matrix(BB_cov)%*%beta_BB_all
      
      h2_EU<-h2*(num_selected_SNP)
      beta_EU_shared<-c(beta_shared[,1])*sqrt(h2_EU/drop(var(EU_plink_causal%*%c(beta_shared[,1]))))
      beta_EU_all<-rep(0,ncol(EU_cov))
      beta_EU_all[which(EU_plink$map$snp.name%in%c(shared_causal_SNP_list))]<-beta_EU_shared
      beta_EU_marginal<-as.matrix(EU_cov)%*%beta_EU_all
      
      y_null_EU<-rnorm(nrow(EU_plink_geno),0,sqrt(1-h2_EU))
      y_null_EU<-y_null_EU-mean(y_null_EU)
      y_null_BB<-rnorm(nrow(BB_plink_geno),0,sqrt(1-h2_BB))
      y_null_BB<-y_null_BB-mean(y_null_BB)
      err_beta_EU<-t(EU_plink_geno)%*%y_null_EU/nrow(EU_plink_geno)
      err_beta_BB<-t(BB_plink_geno)%*%y_null_BB/nrow(BB_plink_geno)
      
      err_beta_EU_scale<-sqrt(nrow(EU_plink_geno)/num_EU)*err_beta_EU
      err_beta_BB_scale<-sqrt(nrow(BB_plink_geno)/num_BB)*err_beta_BB
      
      z_EU<-(beta_EU_marginal+err_beta_EU_scale)*sqrt(num_EU)
      z_BB<-(beta_BB_marginal+err_beta_BB_scale)*sqrt(num_BB)
      
      signal<-rep(0,length(EU_plink$map$snp.name))
      signal[which(EU_plink$map$snp.name%in%shared_causal_SNP_list)]<-3
      
      z_file<-data.frame("CHR" = as.numeric(EU_plink$map$chr[1]),"POS" =EU_plink$map$position ,"RSID" =EU_plink$map$snp.name ,"zscore_1"=z_EU,"zscore_2"=z_BB,"Signal"=signal,"N_1" = num_EU,"N_2" = num_BB)
      
      z_file_out<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num)
      write.table(z_file,z_file_out,col.names = T,row.names = F,quote=F,sep=" ")  
      
      ld_EU_paintor_name<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")
      ld_BB_paintor_name<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")
      fwrite(EU_cov,ld_EU_paintor_name,sep=" ",col.names=F)
      fwrite(BB_cov,ld_BB_paintor_name,sep=" ",col.names=F)
      
      annotation_file<-matrix(rep(1,nrow(z_file)),ncol=1)
      colnames(annotation_file)<-"coding"
      annotation_paintor_name<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".annotations")
      write.table(annotation_file,annotation_paintor_name,col.names = T,row.names = F,quote=F)
      
      input_file_name<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".input")
      write.table(paste0("CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),input_file_name,col.names = F,row.names = F,quote=F)
      
    }
  }
}

