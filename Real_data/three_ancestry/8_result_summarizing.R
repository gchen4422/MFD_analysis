library(ggpubr)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggpmisc)
library(VennDiagram)
library(gridExtra)




##################################################################################################
#
#		Step 1: Get the all result to be analyzed
#
##################################################################################################
res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summary_res/"
system(paste0("mkdir -p ",res_dir))
# List all files
traits_list <- c("BMI", "DBP", "SBP")



annotation_file<-fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/Annotation/Summarized_Annotation_all_updated.txt")
table(annotation_file$brain_ind_eQTL)
all_trait_locus_summary<-c()
region_id = 0
for (trait in traits_list){
  trait_name <- trait
  candidate_region<-fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/summstats_to_finemap/", trait_name, "_risk_loci"))
  colnames(candidate_region)[1] = "CHR"
  colnames(candidate_region)[2] = "MinBP"
  colnames(candidate_region)[3] = "MaxBP"
  region_length_file<-candidate_region%>%summarise(min_count = min(width),max_count = max(width))
  region_length_file$Trait=trait_name
  cat(trait_name)

  result_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/results_",trait_name,"/"),"result/")
  sumstat_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/MF_benchmarking/MFD_revision/results_",trait_name,"/"),"summary_data/")

  #all_gene_region_trait<-seq(region_length_file%>%filter(Trait==trait_name)%>%pull(min_count),region_length_file%>%filter(Trait==trait_name)%>%pull(max_count))
  all_gene_region_trait = 1:nrow(candidate_region)
  ###Check all results are available
  locus_selected_paintor<-which(paste0("loci_",all_gene_region_trait,".mcmc.paintor")%in%list.files(path=result_dir))
  locus_selected_MESuSiE<-which(paste0("MF_result_LOCI_merged_",all_gene_region_trait,".RData")%in%list.files(path=result_dir))
  locus_selected<-intersect(locus_selected_paintor,locus_selected_MESuSiE)
  cat("All regions are availabe",length(all_gene_region_trait)==length(locus_selected))

  all_locus_summary<-c()
  for(region in locus_selected){
    region_id = region_id + 1
    load(paste0(result_dir,"MF_result_LOCI_merged_",region,".RData"))

    WB<-fread(paste0(sumstat_dir,"eur_loci_",region,"_summ"))
    WB$PHENONAME = trait
    BB<-fread(paste0(sumstat_dir,"afr_loci_",region,"_summ"))
    BB$PHENONAME = trait
    EA<-fread(paste0(sumstat_dir,"eas_loci_",region,"_summ"))
    EA$PHENONAME = trait

    # Harmonize EAS alleles to EUR before joining
    overlap_eur <- intersect(WB$SNP, EA$SNP)
    if (length(overlap_eur) > 0) {
      eur_al <- WB[match(overlap_eur, SNP), .(SNP, ALT, REF)]
      eas_al <- EA[match(overlap_eur, SNP), .(SNP, ALT, REF)]
      flip <- (eas_al$ALT == eur_al$REF) & (eas_al$REF == eur_al$ALT)
      flip_snps <- overlap_eur[flip]
      if (length(flip_snps) > 0) {
        idx <- EA$SNP %in% flip_snps
        EA[idx, BETA := -BETA]
        EA[idx, Z    := -Z]
        EA[idx, MAF  := 1 - MAF]
        old_alt <- EA$ALT[idx]
        EA[idx, ALT := REF]
        EA[idx, REF := old_alt]
      }
    }

    # Harmonize remaining EAS alleles to AFR (for SNPs absent in EUR)
    eas_remaining <- setdiff(EA$SNP, WB$SNP)
    overlap_afr <- intersect(BB$SNP, eas_remaining)
    if (length(overlap_afr) > 0) {
      afr_al <- BB[match(overlap_afr, SNP), .(SNP, ALT, REF)]
      eas_al2 <- EA[match(overlap_afr, SNP), .(SNP, ALT, REF)]
      flip2 <- (eas_al2$ALT == afr_al$REF) & (eas_al2$REF == afr_al$ALT)
      flip_snps2 <- overlap_afr[flip2]
      if (length(flip_snps2) > 0) {
        idx2 <- EA$SNP %in% flip_snps2
        EA[idx2, BETA := -BETA]
        EA[idx2, Z    := -Z]
        EA[idx2, MAF  := 1 - MAF]
        old_alt2 <- EA$ALT[idx2]
        EA[idx2, ALT := REF]
        EA[idx2, REF := old_alt2]
      }
    }

    summary_stat_combined <- full_join(
      WB,
      BB,
      by = c("CHR", "POS", "CHR_POS", "SNP", "ALT", "REF", "PHENONAME"),
      suffix = c("_eur", "_afr")
    ) %>%
      full_join(
        EA,
        by = c("CHR", "POS", "CHR_POS", "SNP", "ALT", "REF", "PHENONAME")
      ) %>%
      dplyr::rename(MAF_eas = MAF, BETA_eas = BETA, Z_eas = Z, N_eas = N) %>%
      mutate(Region = region) %>%
      select(c(
        "SNP", "CHR", "POS", "CHR_POS", "PHENONAME", "Region",
        "ALT", "REF", "MAF_eur", "MAF_afr", "MAF_eas",
        "BETA_eur", "BETA_afr", "BETA_eas",
        "Z_eur", "Z_afr", "Z_eas",
        "N_eur", "N_afr", "N_eas"
      )) %>% mutate(PHENONAME = trait) %>% arrange(CHR,POS)

    if(!all(summary_stat_combined$SNP == mfd_result_all$SNP)) {
      warning("Mismatch in Credible Set variants after merge! Check SNP IDs error.")
    }

    #For PAINTOR result without annotation

    locus_summary = summary_stat_combined %>% left_join(mfd_result_all %>% select(-CHR,-POS), by = "SNP")

    all_locus_summary<-rbind(all_locus_summary,locus_summary)
    cat(region)

  }
  ###Adding functional annotations to the data
  annotation_file_order<-annotation_file[match(all_locus_summary$CHR_POS,annotation_file$SNP),]
  all_locus_summary<-cbind(all_locus_summary,annotation_file_order)

  all_trait_locus_summary<-rbind(all_trait_locus_summary,all_locus_summary)

}


all_trait_locus_summary<-all_trait_locus_summary%>%dplyr::select(unique(colnames(.)))
all_trait_locus_summary <- all_trait_locus_summary %>% select(-RSID)
res_all<-all_trait_locus_summary

#res_name<-paste0(res_dir,"res_lipids_pca_bic.RData")
res_name<-paste0(res_dir,"res_all_mf.RData")
save(res_all,file = res_name)
