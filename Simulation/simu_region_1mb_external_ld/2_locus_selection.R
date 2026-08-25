library(data.table)
library(dplyr)
library(snpStats)


find_common_snps <- function(sumstat, ref) {
  
  # Identifying common SNPs between sumstat and ref datasets
  common_SNP <- intersect(sumstat$CHR_POS, ref$CHR_POS)
  
  # Subsetting data frames based on common SNPs and arranging by order of common_SNP
  sumstat <- sumstat %>% 
    filter(CHR_POS %in% common_SNP) %>% 
    arrange(match(CHR_POS, common_SNP))
  
  ref <- ref %>% 
    filter(CHR_POS %in% common_SNP) %>% 
    arrange(match(CHR_POS, common_SNP))
  
  matched_pos<-which((sumstat$REF==ref$REF&sumstat$ALT==ref$ALT)|(sumstat$REF==ref$ALT&sumstat$ALT==ref$REF))
  
  common_SNP<-sumstat[matched_pos,]%>%pull(CHR_POS)
  
  return(common_SNP)
}


# Read locus info and filter by region
selected.100.region<-read.table("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/simulation_locus_info.txt",header=T)
#candidate_region <- filter(selected.100.region, sim_region == region)
candidate_region = selected.100.region

eur_ref = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/eursel50k_nofilter_nodup.bim")
afr_ref = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/changedafr_nodup_flipped.bim")
ref_1kg = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build37/1kg_all_biallelic_unrelated.bim")


#Find snps present in both EUR and AFR ancestry
afr_ref = afr_ref %>% mutate(CHR_POS = paste0(V1,":",V4),REF = V6,ALT = V5) 
eur_ref = eur_ref %>% mutate(CHR_POS = V2,REF = V6,ALT = V5)






#Find multi-allelic snps in African population
multi_allelic = afr_ref$CHR_POS[which(duplicated(afr_ref$CHR_POS) == T)]
eur_biallelic = eur_ref %>% filter(!(V2 %in% multi_allelic))
afr_biallelic = afr_ref %>% filter(!(CHR_POS %in% multi_allelic))

common_SNP_bim = find_common_snps(eur_biallelic, afr_biallelic)

eur_biallelic <- eur_biallelic %>%
  filter(CHR_POS %in% common_SNP_bim)

afr_biallelic <- afr_biallelic %>%
  filter(CHR_POS %in% common_SNP_bim)


eur_ref = eur_biallelic
afr_ref = afr_biallelic



#identify snps to be flipped
flip_indices <- which(eur_ref$V5 == afr_ref$V6 & eur_ref$V6 == afr_ref$V5)

afr_ref$flip = 0
afr_ref$flip[flip_indices] = 1


#hrc_snp_list = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/HRC_snp_list")

ref_1kg = ref_1kg %>% mutate(CHR_POS = paste0(V1,":",V4),REF = V6,ALT = V5) %>% filter(CHR_POS %in% afr_ref$CHR_POS)
multi_allelic = ref_1kg$CHR_POS[which(duplicated(ref_1kg$CHR_POS) == T)]

ref_1kg = ref_1kg %>% filter(!(CHR_POS %in% multi_allelic))


eur_ref<- eur_ref %>%
  filter(CHR_POS %in% ref_1kg$CHR_POS)

afr_ref <- afr_ref %>%
  filter(CHR_POS %in% ref_1kg$CHR_POS)




# ref_1kg_sample = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_all_biallelic_unrelated.fam")
# ref_1kg_sample_toexclude = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/deg2_phase3.king.cutoff.out.id")
# ref_1kg_sample_ans_eur = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/all_phase3.psam") %>% filter(SuperPop == "EUR")
# ref_1kg_sample_ans_afr = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/all_phase3.psam") %>% filter(SuperPop == "AFR")
# ref_1kg_sample_eur = ref_1kg_sample %>% filter(V2 %in% ref_1kg_sample_ans_eur$`#IID`)
# ref_1kg_sample_afr = ref_1kg_sample %>% filter(V2 %in% ref_1kg_sample_ans_afr$`#IID`)
# 
# 
# write.table(c("#IID",ref_1kg_sample_eur$V2), file = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_all_biallelic_unrelated_eur_sample", quote = FALSE, row.names=FALSE, col.names = F)
# write.table(c("#IID",ref_1kg_sample_afr$V2), file = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_all_biallelic_unrelated_afr_sample", quote = FALSE, row.names=FALSE, col.names = F)

# write.table(ref_1kg$V2, file = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_all_biallelic_unrelated_snplist", quote = FALSE, row.names=FALSE, col.names = F)
# write.table(ref_1kg_afr$V2[which(duplicated(ref_1kg_afr$V2) == T)], file = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/snp_to_exclude", quote = FALSE, row.names=FALSE, col.names = F)
#


ref_1kg_eur = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build37/1kg_eur_biallelic_unrelated_filtered.bim")
ref_1kg_afr = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build37/1kg_afr_biallelic_unrelated_filtered.bim")



ref_1kg_eur = ref_1kg_eur %>% mutate(CHR_POS = paste0(V1,":",V4),REF = V6,ALT = V5)
multi_allelic = ref_1kg_eur$CHR_POS[which(duplicated(ref_1kg_eur$V2) == T)]
ref_1kg_eur = ref_1kg_eur %>% filter(!(CHR_POS %in% multi_allelic))

ref_1kg_afr = ref_1kg_afr %>% mutate(CHR_POS = paste0(V1,":",V4),REF = V6,ALT = V5)
multi_allelic = ref_1kg_afr$CHR_POS[which(duplicated(ref_1kg_afr$V2) == T)]
ref_1kg_afr = ref_1kg_afr %>% filter(!(CHR_POS %in% multi_allelic))

eur_ref = eur_ref %>% filter(CHR_POS %in% ref_1kg_eur$CHR_POS)
afr_ref = afr_ref %>% filter(CHR_POS %in% ref_1kg_afr$CHR_POS)



#identify snps to be flipped
flip_indices <- which(eur_ref$V5 == ref_1kg_eur$V6 & eur_ref$V6 == ref_1kg_eur$V5)

ref_1kg_eur$flip = 0
ref_1kg_eur$flip[flip_indices] = 1


# write.table(cbind(ref_1kg_eur$V2,eur_ref$REF), file = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_ref_snp_update.txt", quote = FALSE, row.names=FALSE, col.names = F)
# #The below equivalent to flip allele in summary statistics when running fine-mapping
# command = paste("~/gwas_software/plink2.0/plink2 --bfile /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_eur_biallelic_unrelated_filtered --ref-allele /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_ref_snp_update.txt 2 1 --make-bed --out /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_eur_biallelic_unrelated_filtered_flipped")
# t1 <- system(command, intern = TRUE)
# command = paste("~/gwas_software/plink2.0/plink2 --bfile /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_afr_biallelic_unrelated_filtered --ref-allele /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_ref_snp_update.txt 2 1 --make-bed --out /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_afr_biallelic_unrelated_filtered_flipped")
# t1 <- system(command, intern = TRUE)


ref_1kg_eur = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build37/1kg_eur_biallelic_unrelated_filtered_flipped_clean.bim")
ref_1kg_afr = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build37/1kg_afr_biallelic_unrelated_filtered_flipped_clean.bim")



ref_1kg_eur = ref_1kg_eur %>% mutate(CHR_POS = paste0(V1,":",V4),REF = V6,ALT = V5)
ref_1kg_afr = ref_1kg_afr %>% mutate(CHR_POS = paste0(V1,":",V4),REF = V6,ALT = V5)


eur_ref = eur_ref %>% filter(CHR_POS %in% ref_1kg_eur$CHR_POS)
afr_ref = afr_ref %>% filter(CHR_POS %in% ref_1kg_afr$CHR_POS)




#identify snps to be flipped
flip_indices <- which(eur_ref$V5 == ref_1kg_eur$V6 & eur_ref$V6 == ref_1kg_eur$V5)
flip_indices <- which(afr_ref$V5 == ref_1kg_afr$V6 & afr_ref$V6 == ref_1kg_afr$V5)

ref_1kg_eur$flip = 0
ref_1kg_eur$flip[flip_indices] = 1


# write.table(cbind(ref_1kg_eur$V2,eur_ref$REF), file = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_ref_snp_update_second.txt", quote = FALSE, row.names=FALSE, col.names = F)
# #The below equivalent to flip allele in summary statistics when running fine-mapping
# command = paste("~/gwas_software/plink2.0/plink2 --bfile /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_eur_biallelic_unrelated_filtered_flipped_clean --ref-allele /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_ref_snp_update_second.txt 2 1 --make-bed --out /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_eur_biallelic_unrelated_filtered_flipped_clean_final")
# t1 <- system(command, intern = TRUE)
# command = paste("~/gwas_software/plink2.0/plink2 --bfile /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_afr_biallelic_unrelated_filtered_flipped_clean --ref-allele /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_ref_snp_update_second.txt 2 1 --make-bed --out /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_remove_relatedness/1kg_afr_biallelic_unrelated_filtered_flipped_clean_final")
# t1 <- system(command, intern = TRUE)

ref_1kg_eur = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build37/1kg_eur_biallelic_unrelated_filtered_flipped_clean_final.bim")
ref_1kg_afr = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build37/1kg_afr_biallelic_unrelated_filtered_flipped_clean_final.bim")

# 
# flip_indices <- which(eur_ref$V5 == test_eur$V6 & eur_ref$V6 == test_eur$V5)
# flip_indices <- which(afr_ref$V5 == test_afr$V6 & eur_ref$V6 == test_afr$V5)


ref_1kg_eur = ref_1kg_eur %>% mutate(CHR_POS = paste0(V1,":",V4),REF = V6,ALT = V5)
ref_1kg_afr = ref_1kg_afr %>% mutate(CHR_POS = paste0(V1,":",V4),REF = V6,ALT = V5)

eur_ref = eur_ref %>% filter(CHR_POS %in% ref_1kg_eur$CHR_POS)
afr_ref = afr_ref %>% filter(CHR_POS %in% ref_1kg_afr$CHR_POS)

ref_1kg_eur = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build37/1kg_eur_biallelic_unrelated_filtered_flipped_clean_final_maf_filtered.bim")
ref_1kg_afr = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build37/1kg_afr_biallelic_unrelated_filtered_flipped_clean_final_maf_filtered.bim")


ref_1kg_eur = ref_1kg_eur %>% mutate(CHR_POS = paste0(V1,":",V4),REF = V6,ALT = V5)
ref_1kg_afr = ref_1kg_afr %>% mutate(CHR_POS = paste0(V1,":",V4),REF = V6,ALT = V5)

common_SNP_bim = find_common_snps(ref_1kg_eur, ref_1kg_afr)

ref_1kg_eur = ref_1kg_eur %>% filter(CHR_POS %in% common_SNP_bim)
ref_1kg_afr = ref_1kg_afr %>% filter(CHR_POS %in% common_SNP_bim)


eur_ref = eur_ref %>% filter(CHR_POS %in% ref_1kg_eur$CHR_POS)
afr_ref = afr_ref %>% filter(CHR_POS %in% ref_1kg_afr$CHR_POS)


#write.table(cbind(eur_ref$V2,eur_ref$REF), file = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/ukb_ref_snp_update_afr.txt", quote = FALSE, row.names=FALSE, col.names = F)
# #The below equivalent to flip allele in summary statistics when running fine-mapping
#command = paste("~/gwas_software/plink2.0/plink2 --bfile /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/changedafr_nodup --ref-allele /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/ukb_ref_snp_update_afr.txt 2 1 --make-bed --out /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/changedafr_nodup_flipped")
#t1 <- system(command, intern = TRUE)



# Function to preprocess data by replacing missing values with column means 
# and scaling columns to have zero mean and unit variance
preprocess_data <- function(data_matrix) {
  # Replace missing values with column means and scale columns
  data_matrix <- apply(data_matrix, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x <- scale(x)
    return(x)
  })
  
  return(data_matrix)
}

# Function to compute covariance matrix and write it to a file
compute_cov <- function(data_matrix, file_name) {
  # Compute the covariance matrix
  cov_matrix <- cov2cor(crossprod(data_matrix))
  
  # Ensure the covariance matrix is symmetric
  cov_matrix <- (cov_matrix + t(cov_matrix)) / 2
  
  # Write the covariance matrix to a file
  fwrite(cov_matrix, file_name,sep =" ")
}



chr = 0
for(i in 1:nrow(candidate_region)){
  
  start = candidate_region[i,]$MinBP
  end = candidate_region[i,]$MaxBP
  loci_first = eur_ref[(eur_ref$V4 >= start) & (eur_ref$V4 <= end)]
  loci = loci_first[loci_first$V1 == candidate_region[i,]$CHR,]
  loci_snp = loci$V2
  snplist = paste("loci",i, sep = "_")
  write.table(loci_snp,file = paste0("risk_loci_ld_eur/",snplist), quote = FALSE, row.names=FALSE, col.names = F)
  ldname = paste("./risk_loci_ld_eur/loci", i,sep = "_")
  command = paste("~/gwas_software/plink2/plink --bfile /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/eursel50k_nofilter_nodup --extract",paste0("risk_loci_ld_eur/",snplist)," --allow-no-sex --make-bed --out",ldname)
  t1 <- system(command, intern = TRUE)
  
  #Check LD and summary stats allele consistency
  #info = fread(paste0("risk_loci_ld_eur/loci_",i,".bim"))
  eur_plink <- read.plink(paste0(ldname,".bed"))
  eur_plink_geno <- as(eur_plink$genotypes, "numeric")
  eur_mat = as.matrix(eur_plink_geno)
  eur_mat = preprocess_data(eur_mat)
  compute_cov(eur_mat, paste0(ldname,".ld"))
  
  
  start = candidate_region[i,]$MinBP
  end = candidate_region[i,]$MaxBP
  loci_first = ref_1kg_eur[(ref_1kg_eur$V4 >= start) & (ref_1kg_eur$V4 <= end)]
  loci = loci_first[loci_first$V1 == candidate_region[i,]$CHR,]
  loci_snp = loci$V2
  snplist = paste("loci_1kg",i, sep = "_")
  write.table(loci_snp,file = paste0("risk_loci_ld_eur/",snplist), quote = FALSE, row.names=FALSE, col.names = F)
  ldname = paste("./risk_loci_ld_eur/loci_1kg", i,sep = "_")
  command = paste("~/gwas_software/plink2.0/plink2 --bfile /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build37/1kg_eur_biallelic_unrelated_filtered_flipped_clean_final --extract",paste0("risk_loci_ld_eur/",snplist)," --allow-no-sex --make-bed --out",ldname)
  t1 <- system(command, intern = TRUE)
  # 
  # 
  # 
  # #Check LD and summary stats allele consistency
  #info = fread(paste0("risk_loci_ld_eur/loci_",i,".bim"))
  eur_plink <- read.plink(paste0(ldname,".bed"))
  eur_plink_geno <- as(eur_plink$genotypes, "numeric")
  eur_mat = as.matrix(eur_plink_geno)
  eur_mat = preprocess_data(eur_mat)
  # 
  compute_cov(eur_mat, paste0(ldname,".ld"))
  #save(eur_mat, file = paste0("risk_loci_ld/loci_",i,".Rdata"))
  
  loci_first = afr_ref[(afr_ref$V4 >= start) & (afr_ref$V4 <= end)]
  loci = loci_first[loci_first$V1 == candidate_region[i,]$CHR,]
  loci_snp = loci$V2
  snplist = paste("loci",i, sep = "_")
  write.table(loci_snp,file = paste0("risk_loci_ld_afr/",snplist), quote = FALSE, row.names=FALSE, col.names = F)
  ldname = paste("./risk_loci_ld_afr/loci", i,sep = "_")
  #command = paste("~/gwas_software/plink2.0/plink2 --bfile /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/changed_map_afr --extract",paste0("risk_loci_ld_afr/",snplist)," --allow-no-sex --make-bed --out",ldname)
  command = paste("~/gwas_software/plink2/plink --bfile /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/changedafr_nodup_flipped --keep-allele-order --extract",paste0("risk_loci_ld_afr/",snplist)," --allow-no-sex --make-bed --out",ldname)
  t1 <- system(command, intern = TRUE)
  
  
  #Check LD and summary stats allele consistency
  #info = fread(paste0("risk_loci_ld_afr/loci_",i,".bim"))
  afr_plink <- read.plink(paste0(ldname,".bed"))
  afr_plink_geno <- as(afr_plink$genotypes, "numeric")
  afr_mat = as.matrix(afr_plink_geno)
  afr_mat = preprocess_data(afr_mat)
  
  compute_cov(afr_mat, paste0(ldname,".ld"))
  
  start = candidate_region[i,]$MinBP
  end = candidate_region[i,]$MaxBP
  loci_first = ref_1kg_afr[(ref_1kg_afr$V4 >= start) & (ref_1kg_afr$V4 <= end)]
  loci = loci_first[loci_first$V1 == candidate_region[i,]$CHR,]
  loci_snp = loci$V2
  snplist = paste("loci_1kg",i, sep = "_")
  write.table(loci_snp,file = paste0("risk_loci_ld_afr/",snplist), quote = FALSE, row.names=FALSE, col.names = F)
  ldname = paste("./risk_loci_ld_afr/loci_1kg", i,sep = "_")
  command = paste("~/gwas_software/plink2.0/plink2 --bfile /scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build37/1kg_afr_biallelic_unrelated_filtered_flipped_clean_final --extract",paste0("risk_loci_ld_afr/",snplist)," --allow-no-sex --make-bed --out",ldname)
  t1 <- system(command, intern = TRUE)
  # 
  # 
  # #Check LD and summary stats allele consistency
  # #info = fread(paste0("risk_loci_ld_afr/loci_",i,".bim"))
  afr_plink <- read.plink(paste0(ldname,".bed"))
  afr_plink_geno <- as(afr_plink$genotypes, "numeric")
  afr_mat = as.matrix(afr_plink_geno)
  afr_mat = preprocess_data(afr_mat)
  # 
  compute_cov(afr_mat, paste0(ldname,".ld"))
  #save(afr_mat, file = paste0("risk_loci_ld/loci_",i,".Rdata"))
  

  # if(chr != loci$V1[1]){
  # annotation_file_name<-paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/baselineLF_v2.2.UKB/baselineLF2.2.UKB.",loci$V1[1],".annot.gz")
  # annotation_file<-fread(annotation_file_name)
  # }
  # chr = loci$V1[1]
  # loci_first = annotation_file[(annotation_file$BP >= start) & (annotation_file$BP <= end)]
  # loci = loci_first[loci_first$CHR == candidate_region[i,]$CHR,] %>% mutate(CHR_POS = paste0(CHR,":",BP))
  # snplist_extract = fread(paste0("risk_loci_ld_eur/loci_",i),header = F)
  # loci_sub = loci %>% filter(CHR_POS %in% snplist_extract$V1) %>% select(-CHR_POS)
  # #Remove potential duplication, i.e. snp in same location but with different rsid
  # loci_sub = loci_sub[!duplicated(loci_sub$BP), ]
  # snplist = paste("loci",i, sep = "_")
  # write.table(loci_sub,file = paste0("risk_loci_annot/",snplist,".annot"), quote = FALSE, row.names=FALSE, col.names = T)
  
  
}


