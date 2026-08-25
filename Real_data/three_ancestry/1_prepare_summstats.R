library(data.table)
library(tidyverse)
source("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/utility.R")

#all_sample_1kg = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/all_hg38.psam")
#eas_sample_1kg = all_sample_1kg %>% filter(SuperPop == "EAS")

#fwrite(
#    eas_sample_1kg,
#   file = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/eas_with_related.psam",sep = "\t",
#   quote = FALSE, row.names = FALSE, col.names = TRUE
#)


ref = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/1kg_eas_maf_unrelated_biallelic_nodup.bim")
colnames(ref) = c("CHR","SNP","CM","POS","ALT","REF")
ref <- GWAS_QC(ref, 0.01)  # Do maf qc, remove strand Ambiguous Variants, Remove Multi-allelic Variants, Remove Indel and MHC region SNPs


#test = fread("../eur_w_ld_chr/w_hm3.snplist")
#ref_hm3 = ref %>% filter(SNP %in% test$SNP) %>% select(-3)

#fwrite(ref_hm3,file = "../eur_w_ld_chr/hm3_1kg38_snplist",
#       quote = FALSE, row.names = FALSE, col.names = TRUE,sep = " ")


# fwrite(
#   ref[,c(2)],
#   file = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/1kg_snplist",
#   quote = FALSE, row.names = FALSE, col.names = FALSE
# )


#annot = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/baselineLF_v2.2.UKB/annot_combine/first_three_columns.annot")

# List all files
all_files <- list.files(path = "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/formated_summstats/", pattern = "_format_noMHC.csv.gz", full.names = TRUE)

# Identify eas and EUR files
eas_files <- all_files[grep("_EAS_", all_files)]
eur_files <- all_files[grep("_eur_", all_files)]
afr_files <- all_files[grep("_afr_", all_files)]


# Ensure corresponding AFR/EUR files are matched
for (eas_file in eas_files) {
  # Extract the trait name (assuming it's consistent across AFR and EUR files)
  #trait <- sub("*._afr_.*", "", afr_file)
  trait <- sub("^[^_]+_([^_]+)_.*", "\\1", basename(eas_file))

  # Find the matching EUR file
  eur_file <- eur_files[grep(trait, eur_files)]
  afr_file <- afr_files[grep(trait, afr_files)]

  if (length(eur_file) == 0) {
    message(paste("No EUR file found for EAS file:", afr_file))
    next
  }

  if (length(afr_file) == 0) {
    message(paste("No AFR file found for EAS file:", afr_file))
    next
  }

  # Read AFR summary stats
  afrdata <- fread(afr_file)
  colnames(afrdata)[which(colnames(afrdata) == "BP")] <- "POS"
  colnames(afrdata)[which(colnames(afrdata) == "A1")] <- "ALT"  # This is the effective allele from hail table allele2, see https://saigegit.github.io/SAIGE-doc/docs/single_step2.html
  colnames(afrdata)[which(colnames(afrdata) == "A2")] <- "REF"  # Reference allele
  colnames(afrdata)[which(colnames(afrdata) == "FRQ")] <- "MAF" # With respect to A1
  #afrdata$PVAL = 10^(-1*afrdata$PVAL) #For regenie output give p as -log10p

  # Perform QC and filtering
  afrdata_qced <- GWAS_QC(afrdata, 0.01)  # Do maf qc, remove strand Ambiguous Variants, Remove Multi-allelic Variants, Remove Indel and MHC region SNPs
  afrdata_qced_tidy <- afrdata_qced %>%
    filter(CHR %in% 1:22) %>%
    #filter(MISSINGRATE <= 0.01) %>%
    filter(CHR_POS %in% ref$CHR_POS) %>% arrange(CHR,POS) %>% select(-CHISQ)

  # Read EUR summary stats
  eurdata <- fread(eur_file)
  colnames(eurdata)[which(colnames(eurdata) == "BP")] <- "POS"
  colnames(eurdata)[which(colnames(eurdata) == "A1")] <- "ALT"
  colnames(eurdata)[which(colnames(eurdata) == "A2")] <- "REF"
  colnames(eurdata)[which(colnames(eurdata) == "FRQ")] <- "MAF"
  # lambda <- median(eurdata$CHISQ) / qchisq(0.5, df = 1)                #Perform genomic control for eur HDL and TG summary stats
  # eurdata$CHISQ_gc <- eurdata$CHISQ / lambda
  # eurdata$PVAL <- pchisq(eurdata$CHISQ_gc, df = 1, lower.tail = FALSE)
  #eurdata$PVAL = 10^(-1*eurdata$PVAL) #For regenie output give p as -log10p



  # Perform QC and filtering
  eurdata_qced <- GWAS_QC(eurdata, 0.01)
  eurdata_qced_tidy <- eurdata_qced %>%
    filter(CHR %in% 1:22) %>%
    #filter(MISSINGRATE <= 0.01) %>%
    filter(CHR_POS %in% ref$CHR_POS)  %>% arrange(CHR,POS) %>% select(-CHISQ)


  # Read eas summary stats
  easdata <- fread(eas_file)
  colnames(easdata)[which(colnames(easdata) == "BP")] <- "POS"
  colnames(easdata)[which(colnames(easdata) == "A1")] <- "ALT"  # This is the effective allele from hail table allele2, see https://saigegit.github.io/SAIGE-doc/docs/single_step2.html
  colnames(easdata)[which(colnames(easdata) == "A2")] <- "REF"  # Reference allele
  colnames(easdata)[which(colnames(easdata) == "FRQ")] <- "MAF" # With respect to A1
  #easdata$PVAL = 10^(-1*easdata$PVAL) #For regenie output give p as -log10p

  # Perform QC and filtering
  easdata_qced <- GWAS_QC(easdata, 0.01)  # Do maf qc, remove strand Ambiguous Variants, Remove Multi-allelic Variants, Remove Indel and MHC region SNPs
  easdata_qced_tidy <- easdata_qced %>%
    filter(CHR %in% 1:22) %>%
    #filter(MISSINGRATE <= 0.01) %>%
    filter(CHR_POS %in% ref$CHR_POS) %>% arrange(CHR,POS) #%>% select(-CHISQ)


  # Find common SNPs
  #common_snp <- find_common_snps(eurdata_qced_tidy, afrdata_qced_tidy,easdata_qced_tidy)

  # Filter for common SNPs
  #eurdata_qced_tidy_tofinemap <- eurdata_qced_tidy %>% filter(CHR_POS %in% common_snp)
  #afrdata_qced_tidy_tofinemap <- afrdata_qced_tidy %>% filter(CHR_POS %in% common_snp)
  #easdata_qced_tidy_tofinemap <- easdata_qced_tidy %>% filter(CHR_POS %in% common_snp)

  # Perform allele flipping if necessary
  #afrdata_qced_tidy_tofinemap <- allele_flip(eurdata_qced_tidy_tofinemap, afrdata_qced_tidy_tofinemap)
  #afrdata_qced_tidy_tofinemap$SNP = ref$SNP[match(afrdata_qced_tidy_tofinemap$CHR_POS,ref$CHR_POS)]

  eurdata_qced_tidy$SNP = ref$SNP[match(eurdata_qced_tidy$CHR_POS,ref$CHR_POS)]
  afrdata_qced_tidy$SNP = ref$SNP[match(afrdata_qced_tidy$CHR_POS,ref$CHR_POS)]
  easdata_qced_tidy$SNP = ref$SNP[match(easdata_qced_tidy$CHR_POS,ref$CHR_POS)]


  #eurdata_qced_tidy_tofinemap = eurdata_qced_tidy_tofinemap %>% filter(SNP %in% annot$SNP)
  #afrdata_qced_tidy_tofinemap = afrdata_qced_tidy_tofinemap %>% filter(SNP %in% annot$SNP)


  common_snp_ref_eur <- find_common_snps(eurdata_qced_tidy, ref)
  common_snp_ref_afr <- find_common_snps(afrdata_qced_tidy, ref)
  common_snp_ref_eas <- find_common_snps(easdata_qced_tidy, ref)

  eurdata_qced_tidy = eurdata_qced_tidy %>% filter(CHR_POS %in% common_snp_ref_eur) %>% filter(SNP != "rs1417105319")  #This SNP is a duplication in 1kg
  afrdata_qced_tidy = afrdata_qced_tidy %>% filter(CHR_POS %in% common_snp_ref_afr) %>% filter(SNP != "rs1417105319")  #This SNP is a duplication in 1kg
  easdata_qced_tidy = easdata_qced_tidy %>% filter(CHR_POS %in% common_snp_ref_eas) %>% filter(SNP != "rs1417105319")  #This SNP is a duplication in 1kg

  # Write results
  fwrite(
    eurdata_qced_tidy,
    file = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_EUR_GRCh38_combined_format_noMHC_to_finemap.csv"),
    quote = FALSE, row.names = FALSE, col.names = TRUE
  )


  fwrite(
    eurdata_qced_tidy[,c(1,6)],
    file = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_EUR_GRCh38_combined_format_noMHC_to_finemap_ref_flip"),
    quote = FALSE, row.names = FALSE, col.names = FALSE,sep = " "
  )

  fwrite(
    eurdata_qced_tidy[,c(1)],
    file = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_EUR_GRCh38_combined_format_noMHC_to_finemap_snplist"),
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )

  fwrite(
    afrdata_qced_tidy,
    file = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_AFR_GRCh38_combined_format_noMHC_to_finemap.csv"),
    quote = FALSE, row.names = FALSE, col.names = TRUE
  )

  fwrite(
    afrdata_qced_tidy[,c(1,6)],
    file = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_AFR_GRCh38_combined_format_noMHC_to_finemap_ref_flip"),
    quote = FALSE, row.names = FALSE, col.names = FALSE,sep = " "
  )

  fwrite(
    afrdata_qced_tidy[,c(1)],
    file = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_AFR_GRCh38_combined_format_noMHC_to_finemap_snplist"),
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )

  fwrite(
    easdata_qced_tidy,
    file = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_EAS_GRCh38_combined_format_noMHC_to_finemap.csv"),
    quote = FALSE, row.names = FALSE, col.names = TRUE
  )

  fwrite(
    easdata_qced_tidy[,c(1,6)],
    file = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_EAS_GRCh38_combined_format_noMHC_to_finemap_ref_flip"),
    quote = FALSE, row.names = FALSE, col.names = FALSE,sep = " "
  )

  fwrite(
    easdata_qced_tidy[,c(1)],
    file = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_EAS_GRCh38_combined_format_noMHC_to_finemap_snplist"),
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )

  message(paste("Processed files for trait:", trait))
}



for (eas_file in eas_files) {
  # Extract the trait name (assuming it's consistent across eas and EUR files)
  # For eas
  trait <- sub("^[^_]+_([^_]+)_.*", "\\1", basename(eas_file))
  bfile = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/1kg_eas_maf_unrelated_biallelic_nodup"
  outfile = paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/test_MF_ukb_tpmi/1kg_eas_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped")
  snplist = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_EAS_GRCh38_combined_format_noMHC_to_finemap_snplist")
  ref_flip = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_EAS_GRCh38_combined_format_noMHC_to_finemap_ref_flip")
  command = paste0("~/gwas_software/plink2.0/plink2 --bfile ",bfile," --extract ",snplist ," --ref-allele ",ref_flip," 2 1 --make-bed --out ",outfile)
  t1 <- system(command, intern = TRUE)
  message(paste("Processed files for trait:", trait))

  # Extract the trait name (assuming it's consistent across eas and EUR files)
  # For EUR
  #bfile = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/1kg_eur_maf_unrelated_biallelic_nodup_extracted"
  bfile = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/1kg_eur_maf_unrelated_biallelic_nodup_extracted"
  snplist = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/", trait, "_EUR_GRCh38_combined_format_noMHC_to_finemap_snplist")
  ref_flip = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/", trait, "_EUR_GRCh38_combined_format_noMHC_to_finemap_ref_flip")
  outfile = paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/test_MF/1kg_eur_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped")
  command = paste0("~/gwas_software/plink2.0/plink2 --bfile ",bfile," --extract ",snplist ," --ref-allele ",ref_flip," 2 1 --make-bed --out ",outfile)
  t1 <- system(command, intern = TRUE)


  # Extract the trait name (assuming it's consistent across AFR and EUR files)
  # For AFR
  bfile = "/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/1kg_afr_maf_unrelated_biallelic_nodup_extracted"
  outfile = paste0("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/test_MF/1kg_afr_maf_unrelated_biallelic_nodup_extracted_",trait,"_flipped")
  snplist = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/", trait, "_AFR_GRCh38_combined_format_noMHC_to_finemap_snplist")
  ref_flip = paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/", trait, "_AFR_GRCh38_combined_format_noMHC_to_finemap_ref_flip")
  command = paste0("~/gwas_software/plink2.0/plink2 --bfile ",bfile," --extract ",snplist ," --ref-allele ",ref_flip," 2 1 --make-bed --out ",outfile)
  t1 <- system(command, intern = TRUE)
  message(paste("Processed files for trait:", trait))
}


bim_test_eur = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/test_MF/1kg_eur_maf_unrelated_biallelic_nodup_extracted_BMI_flipped.bim")
all(eurdata_qced_tidy$SNP == bim_test_eur$V2)
all(eurdata_qced_tidy$ALT == bim_test_eur$V5)
all(eurdata_qced_tidy$REF == bim_test_eur$V6)


bim_test_afr = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/test_MF/1kg_afr_maf_unrelated_biallelic_nodup_extracted_BMI_flipped.bim")
all(afrdata_qced_tidy$SNP == bim_test_afr$V2)
all(afrdata_qced_tidy$ALT == bim_test_afr$V5)
all(afrdata_qced_tidy$REF == bim_test_afr$V6)


bim_test_eas = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/lipids_v8/test_MF_ukb_tpmi/1kg_eas_maf_unrelated_biallelic_nodup_extracted_BMI_flipped.bim")
all(easdata_qced_tidy$SNP == bim_test_eas$V2)
all(easdata_qced_tidy$ALT == bim_test_eas$V5)
all(easdata_qced_tidy$REF == bim_test_eas$V6)
