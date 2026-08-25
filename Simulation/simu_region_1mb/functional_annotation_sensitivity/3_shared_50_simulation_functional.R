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

annot_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/functional_annotation_sensitivity/annot/")


# Function to process genotype data
process_geno <- function(geno_data) {
  apply(geno_data, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x <- scale(x - mean(x))
    return(x)
  })
}



for (LD_BLOCK in 1:100){



  # Seed initialization
  set.seed(1128)
  seed_LOCIrator <- sample.int(1000, 100)[LD_BLOCK]
  set.seed(seed_LOCIrator)




  # Reading and processing EU data
  outfile_EU <- paste0(geno_dir_eur, "/loci_", LD_BLOCK)
  EU_plink <- read.plink(paste0(outfile_EU, ".bed"))
  #order_CHR_POS <- EU_plink$map %>% mutate(CHR_POS = paste0(chromosome, "_", position)) %>% pull(CHR_POS)
  EU_plink_geno <- as(EU_plink$genotypes, "numeric")
  EU_plink_geno <- process_geno(EU_plink_geno)

  ld_EU_name <- paste0(LD_dir_eur, "/loci_", LD_BLOCK, ".ld")
  EU_cov <- as.matrix(fread(ld_EU_name))

  # Reading and processing BB data
  outfile_BB <- paste0(geno_dir_afr, "/loci_", LD_BLOCK)
  BB_plink <- read.plink(paste0(outfile_BB, ".bed"))
  BB_plink$map$CHR_POS = paste0(BB_plink$map$chromosome,":",BB_plink$map$position)
  BB_plink_geno <- as(BB_plink$genotypes, "numeric")
  colnames(BB_plink_geno) = BB_plink$map$CHR_POS
  BB_plink_geno <- process_geno(BB_plink_geno)

  ld_BB_name <- paste0(LD_dir_afr, "/loci_", LD_BLOCK, ".ld")
  BB_cov<-as.matrix(fread(ld_BB_name))
  colnames(BB_cov) = BB_plink$map$CHR_POS

  # Check SNP position mismatch
  if(any(colnames(EU_cov)!=EU_plink$map$snp.name)|any(colnames(BB_cov)!=BB_plink$map$snp.name)){
    cat("SNP mismatch of geno and LD reference")
  }
  # Checking allele match
  if (any(EU_plink$map$allele.1 != BB_plink$map$allele.1)) {
    cat("Allele mismatch of the two ancestries")
  }



  # Add functional annotation information to simulation
  annotation_file_name<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_annot/loci_",LD_BLOCK,".annot")
  annot_file<-fread(annotation_file_name)
  evo2_file<-fread(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_annot/CAUSAL_1_LOCI_",LD_BLOCK,"_h2_1_eur_scored.csv"))

  # --- Check evo2 alignment with genotype data ---
  # 1. Row count
  stopifnot("evo2 row count != number of SNPs" = nrow(evo2_file) == nrow(BB_plink$map))

  # 2. SNP order (RSID should match snp.name which is CHR:POS)
  stopifnot("evo2 RSID order != BB_plink snp.name order" = all(evo2_file$RSID == BB_plink$map$snp.name))
  stopifnot("evo2 RSID order != EU_plink snp.name order" = all(evo2_file$RSID == EU_plink$map$snp.name))

  # 3. Allele match with BB_plink (allele.1 = minor/alt, allele.2 = major/ref in snpStats)
  bb_a1_match <- all(evo2_file$A1 == BB_plink$map$allele.1 & evo2_file$A2 == BB_plink$map$allele.2)
  bb_a1_flip  <- all(evo2_file$A1 == BB_plink$map$allele.2 & evo2_file$A2 == BB_plink$map$allele.1)
  cat("evo2 vs BB_plink: A1/A2 direct match =", bb_a1_match, ", flipped match =", bb_a1_flip, "\n")

  # 4. Allele match with EU_plink
  eu_a1_match <- all(evo2_file$A1 == EU_plink$map$allele.1 & evo2_file$A2 == EU_plink$map$allele.2)
  eu_a1_flip  <- all(evo2_file$A1 == EU_plink$map$allele.2 & evo2_file$A2 == EU_plink$map$allele.1)
  cat("evo2 vs EU_plink: A1/A2 direct match =", eu_a1_match, ", flipped match =", eu_a1_flip, "\n")

  # 5. Per-SNP breakdown if not a clean match
  if (!bb_a1_match & !bb_a1_flip) {
    mismatch_idx <- which(!(evo2_file$A1 == BB_plink$map$allele.1 & evo2_file$A2 == BB_plink$map$allele.2) &
                            !(evo2_file$A1 == BB_plink$map$allele.2 & evo2_file$A2 == BB_plink$map$allele.1))
    cat("BB_plink: ", length(mismatch_idx), "SNPs with allele mismatch out of", nrow(evo2_file), "\n")
  }
  if (!eu_a1_match & !eu_a1_flip) {
    mismatch_idx <- which(!(evo2_file$A1 == EU_plink$map$allele.1 & evo2_file$A2 == EU_plink$map$allele.2) &
                            !(evo2_file$A1 == EU_plink$map$allele.2 & evo2_file$A2 == EU_plink$map$allele.1))
    cat("EU_plink: ", length(mismatch_idx), "SNPs with allele mismatch out of", nrow(evo2_file), "\n")
  }

  # Check annot_file alignment with genotype data
  # annot uses rsIDs, plink uses CHR:POS — compare by position instead
  stopifnot("annot row count != BB_plink" = nrow(annot_file) == nrow(BB_plink$map))
  stopifnot("annot row count != EU_plink" = nrow(annot_file) == nrow(EU_plink$map))

  # Position order check
  stopifnot("annot BP order != BB_plink position" = all(annot_file$BP == BB_plink$map$position))
  stopifnot("annot BP order != EU_plink position" = all(annot_file$BP == EU_plink$map$position))

  # Also confirm annot aligns with evo2 file by position
  stopifnot("annot BP order != evo2 POS" = all(annot_file$BP == evo2_file$POS))

  cat("All aligned — safe to do annot_file_subset$evo2_score = evo2_file$evo2_delta_score\n")

  #$annot_file<-annot_file%>%filter(SNP%in%EU_plink$map$snp.name)%>%arrange(match(SNP,EU_plink$map$snp.name))
  annot_file_subset<-annot_file%>%select(c("H3K27ac_Hnisz_common","DHS_Trynka_common","H3K4me3_Trynka_common","Conserved_LindbladToh_common","non_synonymous_common"))
  annot_file_subset$evo2_score  = evo2_file$evo2_delta_score
  annot_file_subset$evo2_score = rank(-annot_file_subset$evo2_score) / nrow(annot_file_subset)
  enrichment_beta<-rep(0,6)
  true_prior_prob<-exp(as.matrix(as.matrix(annot_file_subset)%*%enrichment_beta))/sum(exp(as.matrix(annot_file_subset)%*%enrichment_beta))

  write.table(annot_file_subset,paste0(annot_dir,"Region_",LD_BLOCK,".annot"),quote=F,row.names=F,col.names=T,sep=" ")

  for(num_causal in 1:3){


    wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/functional_annotation_sensitivity/causal_num_",num_causal,"/")
    system(paste0("mkdir -p ",wrk_dir))
    data_dir<-paste0(wrk_dir,"summary_data/w0_startover/")
    result_dir<-paste0(wrk_dir,"result/")
    out_dir<-paste0(wrk_dir,"out/")
    system(paste0("mkdir -p ",data_dir))
    system(paste0("mkdir -p ",result_dir))
    system(paste0("mkdir -p ",out_dir))

    ##sample number of ancestry-specific causal SNP, and then sample from SNP list
    ## 50/50 configuration matching baseline: half shared, half unique
    num_shared_causal_SNP<-c(1,3,5)[num_causal]
    num_unique_causal_SNP<-c(1,3,5)[num_causal]
    num_pop_1_unique<-sample(seq(0,num_unique_causal_SNP),1)
    num_pop_2_unique<-num_unique_causal_SNP-num_pop_1_unique

    num_selected_SNP<-num_shared_causal_SNP+num_unique_causal_SNP

    causal_SNP_list<-BB_plink$map$snp.name[sample(seq(1,nrow(BB_cov)),num_selected_SNP,replace=F,true_prior_prob)]
    shared_causal_SNP_list<-sample(causal_SNP_list,num_shared_causal_SNP)
    pop_1_unique<-sample(setdiff(causal_SNP_list,shared_causal_SNP_list),num_pop_1_unique)
    pop_2_unique<-setdiff(setdiff(causal_SNP_list,shared_causal_SNP_list),pop_1_unique)

    # Draw raw betas ONCE before h2 loop (both h2 conditions share same effect directions)
    beta_shared<-rmvnorm(num_shared_causal_SNP,mean=rep(0,2),sigma=matrix(c(1,0.8,0.8,1),ncol=2,nrow=2))
    beta_pop_1<-rnorm(length(pop_1_unique),0,1)
    beta_pop_2<-rnorm(length(pop_2_unique),0,1)

    print(LD_BLOCK)
    for(h2_num in 1:2){

      h2= c(1e-4,2e-4)[h2_num]

      # Build full beta vectors and apply joint scaling
      num_EUR_causal<-num_shared_causal_SNP+num_pop_1_unique
      num_AFR_causal<-num_shared_causal_SNP+num_pop_2_unique

      h2_EUR<-h2*num_EUR_causal
      h2_AFR<-h2*num_AFR_causal

      # EUR betas: shared + pop_1_unique
      beta_EUR_raw<-rep(0,nrow(EU_cov))
      beta_EUR_raw[which(EU_plink$map$snp.name%in%shared_causal_SNP_list)]<-beta_shared[,1]
      if(num_pop_1_unique>0){
        beta_EUR_raw[which(EU_plink$map$snp.name%in%pop_1_unique)]<-beta_pop_1
      }
      if(num_EUR_causal>0){
        beta_EUR_all<-beta_EUR_raw*sqrt(h2_EUR/drop(var(EU_plink_geno%*%beta_EUR_raw)))
      }else{
        beta_EUR_all<-beta_EUR_raw
      }

      # AFR betas: shared + pop_2_unique
      beta_AFR_raw<-rep(0,nrow(BB_cov))
      beta_AFR_raw[which(BB_plink$map$snp.name%in%shared_causal_SNP_list)]<-beta_shared[,2]
      if(num_pop_2_unique>0){
        beta_AFR_raw[which(BB_plink$map$snp.name%in%pop_2_unique)]<-beta_pop_2
      }
      if(num_AFR_causal>0){
        beta_AFR_all<-beta_AFR_raw*sqrt(h2_AFR/drop(var(BB_plink_geno%*%beta_AFR_raw)))
      }else{
        beta_AFR_all<-beta_AFR_raw
      }

      # Signal vector
      signal<-rep(0,nrow(EU_cov))
      signal[which(EU_plink$map$snp.name%in%shared_causal_SNP_list)]<-3
      signal[which(EU_plink$map$snp.name%in%pop_1_unique)]<-1
      signal[which(EU_plink$map$snp.name%in%pop_2_unique)]<-2

      beta_EUR_marginal<-as.matrix(EU_cov)%*%beta_EUR_all
      beta_AFR_marginal<-as.matrix(BB_cov)%*%beta_AFR_all

      h2_causal_EUR<-drop(var(EU_plink_geno%*%beta_EUR_all))
      h2_causal_AFR<-drop(var(BB_plink_geno%*%beta_AFR_all))
      y_null_EU<-rnorm(nrow(EU_plink_geno),0,sqrt(1-h2_causal_EUR))
      y_null_EU<-y_null_EU-mean(y_null_EU)
      y_null_BB<-rnorm(nrow(BB_plink_geno),0,sqrt(1-h2_causal_AFR))
      y_null_BB<-y_null_BB-mean(y_null_BB)
      err_beta_EU<-t(EU_plink_geno)%*%y_null_EU/nrow(EU_plink_geno)
      err_beta_BB<-t(BB_plink_geno)%*%y_null_BB/nrow(BB_plink_geno)

      err_beta_EU_scale<-sqrt(nrow(EU_plink_geno)/num_EU)*err_beta_EU
      err_beta_BB_scale<-sqrt(nrow(BB_plink_geno)/num_BB)*err_beta_BB

      z_EU<-(beta_EUR_marginal+err_beta_EU_scale)*sqrt(num_EU)
      z_BB<-(beta_AFR_marginal+err_beta_BB_scale)*sqrt(num_BB)

      z_file<-data.frame("CHR" = as.numeric(EU_plink$map$chr[1]),"POS" =EU_plink$map$position ,"RSID" =EU_plink$map$snp.name ,"zscore_1"=z_EU,"zscore_2"=z_BB,"Signal"=signal,"N_1" = num_EU,"N_2" = num_BB)

      z_file_out<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num)
      write.table(z_file,z_file_out,col.names = T,row.names = F,quote=F,sep=" ")

      #ld_EU_paintor_name<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")
      #ld_BB_paintor_name<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")
      #fwrite(EU_cov,ld_EU_paintor_name,sep=" ",col.names=F)
      #fwrite(BB_cov,ld_BB_paintor_name,sep=" ",col.names=F)
    }
  }

}
