library(data.table)
library(dplyr)
library(rtracklayer)
library(IRanges)
library(GenomicRanges)
library(tidyr)
wrk_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/Annotation/"
system(paste0("mkdir -p ",wrk_dir))



traits_list <- c("BMI", "DBP", "SBP")
for (trait in traits_list) {

  candidate_region <- fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_risk_loci"))
  data_dir <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/results_", trait, "/summary_data/")
  out_dir  <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/results_", trait, "/out/")

  summary_stat_list_eur <- list()
  summary_stat_list_afr <- list()
  summary_stat_list_eas <- list()
  summary_stat_list <- list()

  for (i in 1:nrow(candidate_region)) {

    summary_stat_eur <- fread(paste0(data_dir, "eur_loci_", i, "_summ"))
    summary_stat_eur$PHENONAME <- trait
    summary_stat_afr <- fread(paste0(data_dir, "afr_loci_", i, "_summ"))
    summary_stat_afr$PHENONAME <- trait
    summary_stat_eas <- fread(paste0(data_dir, "eas_loci_", i, "_summ"))
    summary_stat_eas$PHENONAME <- trait

    # Harmonize EAS alleles to EUR before joining
    overlap_eur <- intersect(summary_stat_eur$SNP, summary_stat_eas$SNP)
    if (length(overlap_eur) > 0) {
      eur_al <- summary_stat_eur[match(overlap_eur, SNP), .(SNP, ALT, REF)]
      eas_al <- summary_stat_eas[match(overlap_eur, SNP), .(SNP, ALT, REF)]
      flip <- (eas_al$ALT == eur_al$REF) & (eas_al$REF == eur_al$ALT)
      flip_snps <- overlap_eur[flip]
      if (length(flip_snps) > 0) {
        idx <- summary_stat_eas$SNP %in% flip_snps
        summary_stat_eas[idx, BETA := -BETA]
        summary_stat_eas[idx, Z    := -Z]
        summary_stat_eas[idx, MAF  := 1 - MAF]
        old_alt <- summary_stat_eas$ALT[idx]
        summary_stat_eas[idx, ALT := REF]
        summary_stat_eas[idx, REF := old_alt]
      }
    }

    # Harmonize remaining EAS alleles to AFR (for SNPs absent in EUR)
    eas_remaining <- setdiff(summary_stat_eas$SNP, summary_stat_eur$SNP)
    overlap_afr <- intersect(summary_stat_afr$SNP, eas_remaining)
    if (length(overlap_afr) > 0) {
      afr_al <- summary_stat_afr[match(overlap_afr, SNP), .(SNP, ALT, REF)]
      eas_al2 <- summary_stat_eas[match(overlap_afr, SNP), .(SNP, ALT, REF)]
      flip2 <- (eas_al2$ALT == afr_al$REF) & (eas_al2$REF == afr_al$ALT)
      flip_snps2 <- overlap_afr[flip2]
      if (length(flip_snps2) > 0) {
        idx2 <- summary_stat_eas$SNP %in% flip_snps2
        summary_stat_eas[idx2, BETA := -BETA]
        summary_stat_eas[idx2, Z    := -Z]
        summary_stat_eas[idx2, MAF  := 1 - MAF]
        old_alt2 <- summary_stat_eas$ALT[idx2]
        summary_stat_eas[idx2, ALT := REF]
        summary_stat_eas[idx2, REF := old_alt2]
      }
    }

    # Now ALT/REF are aligned — join on all keys including ALT, REF
    summary_stat_combined <- full_join(
      summary_stat_eur,
      summary_stat_afr,
      by = c("CHR", "POS", "CHR_POS", "SNP", "ALT", "REF", "PHENONAME"),
      suffix = c("_eur", "_afr")
    ) %>%
      full_join(
        summary_stat_eas,
        by = c("CHR", "POS", "CHR_POS", "SNP", "ALT", "REF", "PHENONAME")
      ) %>%
      dplyr::rename(MAF_eas = MAF, BETA_eas = BETA, Z_eas = Z, N_eas = N) %>%
      mutate(Region = i) %>%
      select(
        SNP, CHR, POS, CHR_POS, PHENONAME, Region,
        ALT, REF,
        MAF_eur, MAF_afr, MAF_eas,
        BETA_eur, BETA_afr, BETA_eas,
        Z_eur, Z_afr, Z_eas,
        N_eur, N_afr, N_eas
      )

    summary_stat_eur_save <- summary_stat_eur %>% mutate(Region = i)
    summary_stat_afr_save <- summary_stat_afr %>% mutate(Region = i)
    summary_stat_eas_save <- summary_stat_eas %>% mutate(Region = i)

    summary_stat_list_eur[[i]] <- summary_stat_eur_save
    summary_stat_list_afr[[i]] <- summary_stat_afr_save
    summary_stat_list_eas[[i]] <- summary_stat_eas_save
    summary_stat_list[[i]] <- summary_stat_combined
  }

  summary_stat_all_eur <- rbindlist(summary_stat_list_eur)
  summary_stat_all_afr <- rbindlist(summary_stat_list_afr)
  summary_stat_all_eas <- rbindlist(summary_stat_list_eas)
  summary_stat_all     <- rbindlist(summary_stat_list)

  fwrite(summary_stat_all_eur, paste0(out_dir, trait, "_REF_EUR"), col.names = TRUE, row.names = FALSE, quote = FALSE)
  fwrite(summary_stat_all_afr, paste0(out_dir, trait, "_REF_AFR"), col.names = TRUE, row.names = FALSE, quote = FALSE)
  fwrite(summary_stat_all_eas, paste0(out_dir, trait, "_REF_EAS"), col.names = TRUE, row.names = FALSE, quote = FALSE)
  fwrite(summary_stat_all,     paste0(out_dir, trait, "_REF_ALL"), col.names = TRUE, row.names = FALSE, quote = FALSE)
}

summary_stat_list_all_loci <- list()

for (trait in traits_list){


  #trait <- gsub(".*/(.*?)_.*", "\\1", afr_file)
  candidate_region <- fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait, "_risk_loci"))
  out_dir  <- paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/results_", trait, "/out/")

  summary_stat_trait = fread(paste0(out_dir,trait,"_REF_ALL")) %>%mutate(ID = CHR_POS)
  summary_stat_list_all_loci[[trait]] <- summary_stat_trait



}

summary_stat_all_combined <- bind_rows(summary_stat_list_all_loci)
fwrite(summary_stat_all_combined,paste0(wrk_dir,"EU_BB_EA_all_combined_updated"),col.names = T,row.names = F,quote=F)
#fwrite(summary_stat_all_combined,paste0(wrk_dir,"EU_all_combined_bmi"),col.names = T,row.names = F,quote=F)





###############For annotation#####################
all_dat<-fread(paste0(wrk_dir,"EU_BB_EA_all_combined_updated"))
#all_dat<-fread(paste0(wrk_dir,"EU_all_combined_bp"))
all_dat<-all_dat%>%distinct(CHR_POS,.keep_all = TRUE)%>%mutate(ID = CHR_POS)

vcf_file <- all_dat %>% select(CHR, POS, ID, REF, ALT, CHR_POS) %>% arrange(CHR, POS)
annotation_out_name<-paste0(wrk_dir,"EU_BB_EA_all_updated.vcf")
#annotation_out_name<-paste0(wrk_dir,"EU_AA_bp.vcf")
write.table(vcf_file,annotation_out_name,col.names=T,row.names=F,quote=F)
system(paste0("sed -i '1s/^/#/' ",annotation_out_name))
system(paste0("sed -i '1i###fileformat=VCFv4.0 \' ",annotation_out_name))
most_severe_out<-paste0(wrk_dir,"Most_Severe")
#vep_cmd<-paste0("/net/fantasia/home/borang/software/ensembl-vep/vep -i ",annotation_out_name," --cache --offline --most_severe --force_overwrite -o ",most_severe_out)
#system(vep_cmd)


loftee_out<-paste0(wrk_dir,"Loftee")
#vep_loftee_cmd<-paste0("/net/fantasia/home/borang/software/ensembl-vep/vep -i ",annotation_out_name," --cache --offline --plugin LoF,loftee_path:/net/fantasia/home/borang/.vep/Plugins/loftee,human_ancestor_fa:/net/fantasia/home/borang/.vep/Plugins/human_ancestor.fa.gz,phylocsf_data:/net/fantasia/home/borang/.vep/Plugins/phylocsf_gerp.sql,gerp_file:/net/fantasia/home/borang/.vep/Plugins/loftee/GERP_scores.final.sorted.txt.gz -dir_plugins /net/fantasia/home/borang/.vep/Plugins/loftee --force_overwrite -o ",loftee_out)
#system(vep_loftee_cmd)
loftee_out_HC<-paste0(wrk_dir,"Loftee_HC")
loftee_out_LC<-paste0(wrk_dir,"Loftee_LC")
hc_cmd<-paste0("cat ",loftee_out,"|grep LoF=HC>",loftee_out_HC)
system(hc_cmd)
lc_cmd<-paste0("cat ",loftee_out,"|grep LoF=LC>",loftee_out_LC)
system(lc_cmd)
hc_annotation<-read.table(loftee_out_HC,header=F)
lc_annotation<-read.table(loftee_out_LC,header=F)
hc_snp<-unique(hc_annotation$V1)
lc_snp<-unique(lc_annotation$V1)
most_severe_annotation<-fread(most_severe_out,skip=42)
colnames(most_severe_annotation)[1]<-"variant"

missense_snp<-c(hc_snp,lc_snp,most_severe_annotation%>%filter(Consequence=="missense_variant")%>%pull(variant))
synonymous_snp<-most_severe_annotation%>%filter(!(variant%in%missense_snp),Consequence=="synonymous_variant")%>%pull(variant)
utr_5_snp<-most_severe_annotation%>%filter(!(variant%in%c(missense_snp,synonymous_snp)),Consequence=="5_prime_UTR_variant")%>%pull(variant)
utr_3_snp<-most_severe_annotation%>%filter(!(variant%in%c(missense_snp,synonymous_snp,utr_5_snp)),Consequence=="3_prime_UTR_variant")%>%pull(variant)



Promotor_CRE_all<-c()
for(chr in 1:22){
  Promotor_CRE_all<-rbind(Promotor_CRE_all,fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/Annotation/LDSC_base_annotation/baselineLD.",chr,".annot.gz"))) #This annotation is on build 38
}
Promotor_CRE_all<-Promotor_CRE_all%>%mutate(CHR_POS=paste0(CHR,"_",BP))%>%filter(CHR_POS%in%most_severe_annotation$variant)
promotor_snp<-Promotor_CRE_all%>%filter(!(CHR_POS%in%c(missense_snp,synonymous_snp,utr_5_snp,utr_3_snp)),Promoter_UCSC.extend.500==1)%>%pull(CHR_POS)
cre_snp<-Promotor_CRE_all%>%filter(!(CHR_POS%in%c(missense_snp,synonymous_snp,utr_5_snp,utr_3_snp,promotor_snp)),H3K27ac_Hnisz.extend.500==1,DHS_Trynka.extend.500==1)%>%pull(CHR_POS)




annotation_data<-data.frame("SNP"=most_severe_annotation$variant)
annotation_data<-annotation_data%>%mutate(missense = ifelse(SNP%in%missense_snp,1,0),synonymous = ifelse(SNP%in%synonymous_snp,1,0),utr_5 = ifelse(SNP%in%utr_5_snp,1,0),utr_3 = ifelse(SNP%in%utr_3_snp,1,0),promotor = ifelse(SNP%in%promotor_snp,1,0),CRE = ifelse(SNP%in%cre_snp,1,0))

apply(annotation_data[,2:ncol(annotation_data)],2,table)


annotation_data<-annotation_data%>%separate(SNP, into = c("CHR", "POS"), sep = "_", convert = TRUE,remove = FALSE)

write.table(annotation_data,paste0(wrk_dir,"/Summarized_Annotation_all_updated.txt"),col.names=T,row.names=F,quote=F,sep=" ")







wrk_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/Annotation/"

format_eqtl<-function(dat,CHR_POS){
  variant_gh38<-do.call(rbind,lapply(dat$variant_id,function(x)unlist(strsplit(x,"_"))))
  df <- data.frame(cbind(variant_gh38[,1], variant_gh38[,2], variant_gh38[,2]))
  colnames(df) <- c('chr', 'start', 'end')
  df$chr<-gsub("chr","",df$chr)
  df$chr_pos<-paste0(df$chr,"_",df$start)
  gwas_eqtl_intersect<-intersect(df$chr_pos,CHR_POS)
  return(gwas_eqtl_intersect)
}

annotation_data = fread(paste0(wrk_dir,"/Summarized_Annotation_all_updated.txt"))

############################

# For brain eqtl

############################


# List of brain tissue names (without file extension)
tissues <- c(
  "Brain_Amygdala",
  "Brain_Cortex",
  "Brain_Putamen_basal_ganglia",
  "Brain_Anterior_cingulate_cortex_BA24",
  "Brain_Frontal_Cortex_BA9",
  "Brain_Spinal_cord_cervical_c-1",
  "Brain_Caudate_basal_ganglia",
  "Brain_Hippocampus",
  "Brain_Substantia_nigra",
  "Brain_Cerebellar_Hemisphere",
  "Brain_Hypothalamus",
  "Brain_Cerebellum",
  "Brain_Nucleus_accumbens_basal_ganglia"
)

# Base directory where the files are located
base_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MESuSiE_inf_all_by_all/combined/formatted/formatted_sumstats/eQTL/GTEx_Analysis_v8_eQTL_independent/"

# Create an empty list to store formatted eQTL results for each tissue
formatted_eqtls <- list()

# Loop over each tissue file
for (tissue in tissues) {
  # Construct the full file path (assuming the naming pattern is consistent)
  file_path <- file.path(base_dir, paste0(tissue, ".v8.independent_eqtls.txt.gz"))

  # Read the eQTL file
  eqtl_data <- fread(file_path)

  # Format the eQTL data using the SNP annotation data
  formatted_eqtls[[tissue]] <- format_eqtl(eqtl_data, annotation_data$SNP)

  # Optionally, print a message to track progress
  message("Processed: ", tissue)
}

# Now, 'formatted_eqtls' is a named list where each element corresponds to the formatted eQTL results for a brain tissue.

all_brain_eqtl_snps <- unique(unlist(lapply(formatted_eqtls, function(x) x)))

annotation_data<-annotation_data%>%mutate(brain_ind_eQTL = ifelse(SNP%in%all_brain_eqtl_snps,1,0))
table(annotation_data$brain_ind_eQTL)


#annotation_data = annotation_data %>% select(-Human_Promoter_Villar_ExAC.flanking.500_common,-brain_eQTL,-Human_Promoter_Villar_ExAC_common)
############################

# For heart eqtl

############################

# List of brain tissue names (without file extension)
tissues <- c(
  "Heart_Atrial_Appendage",
  "Heart_Left_Ventricle"
)

# Base directory where the files are located
base_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MESuSiE_inf_all_by_all/combined/formatted/formatted_sumstats/eQTL/GTEx_Analysis_v8_eQTL_independent/"

# Create an empty list to store formatted eQTL results for each tissue
formatted_eqtls <- list()

# Loop over each tissue file
for (tissue in tissues) {
  # Construct the full file path (assuming the naming pattern is consistent)
  file_path <- file.path(base_dir, paste0(tissue, ".v8.independent_eqtls.txt.gz"))

  # Read the eQTL file
  eqtl_data <- fread(file_path)

  # Format the eQTL data using the SNP annotation data
  formatted_eqtls[[tissue]] <- format_eqtl(eqtl_data, annotation_data$SNP)

  # Optionally, print a message to track progress
  message("Processed: ", tissue)
}

# Now, 'formatted_eqtls' is a named list where each element corresponds to the formatted eQTL results for a brain tissue.

all_heart_eqtl_snps <- unique(unlist(lapply(formatted_eqtls, function(x) x)))




annotation_data<-annotation_data%>%mutate(heart_ind_eQTL = ifelse(SNP%in%all_heart_eqtl_snps,1,0))
table(annotation_data$heart_ind_eQTL)
which(is.na(annotation_data))


############################

# For artery eqtl

############################


# List of brain tissue names (without file extension)
tissues <- c(
  "Artery_Aorta",
  "Artery_Coronary"
  #"Artery_Tibial"
)

# Base directory where the files are located
base_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MESuSiE_inf_all_by_all/combined/formatted/formatted_sumstats/eQTL/GTEx_Analysis_v8_eQTL_independent/"

# Create an empty list to store formatted eQTL results for each tissue
formatted_eqtls <- list()

# Loop over each tissue file
for (tissue in tissues) {
  # Construct the full file path (assuming the naming pattern is consistent)
  file_path <- file.path(base_dir, paste0(tissue, ".v8.independent_eqtls.txt.gz"))

  # Read the eQTL file
  eqtl_data <- fread(file_path)

  # Format the eQTL data using the SNP annotation data
  formatted_eqtls[[tissue]] <- format_eqtl(eqtl_data, annotation_data$SNP)

  # Optionally, print a message to track progress
  message("Processed: ", tissue)
}

# Now, 'formatted_eqtls' is a named list where each element corresponds to the formatted eQTL results for a brain tissue.

all_artery_eqtl_snps <- unique(unlist(lapply(formatted_eqtls, function(x) x)))




annotation_data<-annotation_data%>%mutate(artery_ind_eQTL = ifelse(SNP%in%all_artery_eqtl_snps,1,0))
table(annotation_data$artery_ind_eQTL)






##Update for rsid

ref = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/1kg_build38/flipped_finemap_allbyall/1kg_afr_maf_unrelated_biallelic_nodup_extracted.bim")
colnames(ref) = c("CHR","SNP","CM","POS","ALT","REF")
ref$CHR_POS = paste0(ref$CHR,"_",ref$POS)
annotation_data$RSID = ref$SNP[match(annotation_data$SNP,ref$CHR_POS)]


#write.table(annotation_data,paste0(wrk_dir,"/Summarized_Annotation_bmi.txt"),col.names=T,row.names=F,quote=F,sep=" ")

write.table(annotation_data,paste0(wrk_dir,"/Summarized_Annotation_all_updated"),col.names=T,row.names=F,quote=F,sep=" ")




#############################################################################################################################
annotation_data = fread(paste0(wrk_dir,"/Summarized_Annotation_all_updated"))
snplist_annot = fread("/scratch/negishi/chen4422/UKBB_pheSVD-main/ukb_geno/in_sample_ld/baselineLF_v2.2.UKB/annot_combine/combined_baselineLF2.2.UKB.annot") %>% filter(SNP %in% annotation_data$RSID)
snplist_annot_sub = snplist_annot %>% select(SNP,non_synonymous_common) %>% filter(non_synonymous_common == 1)
annotation_data = annotation_data %>% mutate(non_synonymous_common = ifelse(RSID %in% snplist_annot_sub$SNP,1,0))
apply(annotation_data[,c(4:12,14)],2,table)

write.table(annotation_data,paste0(wrk_dir,"/Summarized_Annotation_all_updated.txt"),col.names=T,row.names=F,quote=F,sep=" ")
