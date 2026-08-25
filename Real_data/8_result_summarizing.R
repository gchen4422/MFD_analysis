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
res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summary_res/"
system(paste0("mkdir -p ",res_dir))
# List all files
traits_list <- c("BMI", "DBP", "SBP")



annotation_file<-fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/Annotation/Summarized_Annotation_all_updated")
table(annotation_file$brain_ind_eQTL)
all_trait_locus_summary<-c()
region_id = 0
for (trait in traits_list){
  trait_name <- trait
  candidate_region<-fread(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/", trait_name, "_risk_loci"))
  colnames(candidate_region)[1] = "CHR"
  colnames(candidate_region)[2] = "MinBP"
  colnames(candidate_region)[3] = "MaxBP"
  region_length_file<-candidate_region%>%summarise(min_count = min(width),max_count = max(width))
  region_length_file$Trait=trait_name
  cat(trait_name)
  
  result_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/results_",trait_name,"/"),"result/")
  sumstat_dir<-paste0(paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/results_",trait_name,"/"),"summary_data/")
  
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
    
    summary_stat_combined <- full_join(
      WB, 
      BB, 
      by = c("CHR", "POS", "CHR_POS", "SNP", "ALT", "REF", "PHENONAME"), 
      suffix = c("_eur", "_afr")
    ) %>%
      mutate(Region = region) %>%
      select(c(
        "SNP", "CHR", "POS", "CHR_POS", "PHENONAME", "Region", 
        "ALT", "REF", "MAF_eur", "MAF_afr", "BETA_eur", "BETA_afr", 
        "Z_eur", "Z_afr", "N_eur", "N_afr"
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
res_all<-all_trait_locus_summary

#res_name<-paste0(res_dir,"res_lipids_pca_bic.RData")
res_name<-paste0(res_dir,"res_all_mf.RData")
save(res_all,file = res_name)











res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/summary_res/"


load(paste0(res_dir,"res_bmi_pca_bic_second_archive.RData"))
#annotation_data = fread(paste0(wrk_dir,"/Summarized_Annotation.txt"))
test = paste0(res_all$CHR,"_",res_all$POS)
all(test == annotation_data$SNP)
res_all$brain_ind_eQTL = annotation_data$brain_ind_eQTL
table(res_all$brain_ind_eQTL)
res_all$CADD = annotation_data$CADD
res_name<-paste0(res_dir,"res_bmi_pca_bic_second.RData")
save(res_all,file = res_name)
