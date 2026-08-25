library(ggpubr)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggpmisc)
library(VennDiagram)
library(gridExtra)
library(ggbreak)
library(DescTools)
library(coin)
library(susieR)
library(ggrepel)
library(stringr)

res_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summary_res/"
load(paste0(res_dir,"res_all_mf_updated.RData"))   

res_all %>% dplyr::group_by(PHENONAME) %>% dplyr::summarise(n_regions = dplyr::n_distinct(Region), .groups = "drop")
res_all$Region = ifelse(res_all$PHENONAME == "DBP", res_all$Region+ 204, res_all$Region)
res_all$Region = ifelse(res_all$PHENONAME == "SBP", res_all$Region+ 243, res_all$Region)
res_all$CS <- ifelse(res_all$CS != 0, 1, res_all$CS)
res_all$ASV <- ifelse(is.na(res_all$MAF_eur) | is.na(res_all$MAF_afr), 1, 0)


############################################################################################################
#
#   Stats to report AS-Vs
#
#############################################################################################################

res_all %>% filter(PIP_Either >0.5, ASV == 1) %>% nrow()
res_all %>% filter(MESuSiE_PIP_Either  >0.5, , ASV == 0) %>% nrow()
res_all %>% filter(SuSiE_PIP_Either    >0.5, ASV == 1) %>% nrow()
res_all %>% filter(SuSiE_merged_PIP_Either     >0.5, ASV == 1) %>% nrow()


res_all %>% filter(MESuSiE_cs ==1 ) %>% nrow()
res_all %>% filter(CS == 1, ASV == 1) %>% nrow()
res_all %>% filter(SuSiE_cs == 1, ASV == 1) %>% nrow()
res_all %>% filter(SuSiE_metal_cs == 1, ASV == 1) %>% nrow()

############################################################################################################
#
#   Stats to report
#
#############################################################################################################


# Count unique regions
unique_regions <- res_all %>% summarise(n_unique_regions = n_distinct(Region))
print(unique_regions)
# n_unique_regions
# 1              108
# Summarize the number of SNPs by region
snp_summary <- res_all %>%group_by(Region) %>% summarise(NSNP = n()) %>%summarise(summary(NSNP))
print(snp_summary)
# A tibble: 6 × 1
#`summary(NSNP)`
#<table[1d]>    
#1  1908.000      
#2  3277.000      
#3  3719.000      
#4  4169.976      
#5  4443.000      
#6 20469.000 




# Calculate position length by region and summarize
pos_length_summary <- res_all %>%group_by(Region) %>%summarise(max_pos = max(as.numeric(POS)),min_pos = min(as.numeric(POS)),pos_length = as.numeric(max_pos) - as.numeric(min_pos)) %>%summarise(summary(pos_length / (1024 * 1024)))
print(pos_length_summary)
# A tibble: 6 × 1
# `summary(pos_length/(1024 * 1024))`
# <table[1d]>                        
#1 0.9519548                          
#2 0.9594374                          
#3 1.0095758                          
#4 1.1934466                          
#5 1.1740952                          
#6 4.8323755  


# Count unique regions by PHENONAME
unique_regions_by_PHENONAME <- res_all %>%group_by(PHENONAME) %>%summarise(n_unique_regions = n_distinct(Region))
print(unique_regions_by_PHENONAME)
# A tibble: 3 × 2
#PHENONAME n_unique_regions
#<chr>                <int>
#1 BMI                    204
#2 DBP                     39
#3 SBP                     46



# 95% credible set size 
# 675, 673, 338 Regions detected by MESuSiE,SuSiE, and Paintor with a median 95% credible set size of 13, 16, 46.5. 
region_cs<-res_all%>%group_by(Region) %>% summarise(across(c("CS","MESuSiE_cs", "SuSiE_cs","Paintor_cs","SuSiEx_cs","XMAP_cs","CARMAX_cs","SuSiE_metal_cs"), ~ sum(.x, na.rm = TRUE)))
region_cs%>%summarise(sum(CS!=0),sum(MESuSiE_cs!=0),sum(SuSiE_cs!=0),sum(Paintor_cs!=0),sum(SuSiEx_cs!=0),sum(XMAP_cs !=0),sum(CARMAX_cs !=0),sum(SuSiE_metal_cs !=0))
res_all %>%
  group_by(Region) %>%
  summarise(
    across(c("CS","MESuSiE_cs", "SuSiE_cs","Paintor_cs","SuSiEx_cs","XMAP_cs","CARMAX_cs","SuSiE_metal_cs"), 
           ~ sum(.x, na.rm = TRUE), .names = "sum_{col}")) %>%
  summarise(
    MFD= median(sum_CS[sum_CS != 0], na.rm = TRUE),
    MESuSiE_cs_median = median(sum_MESuSiE_cs[sum_MESuSiE_cs != 0], na.rm = TRUE),
    SuSiE_cs_median = median(sum_SuSiE_cs[sum_SuSiE_cs != 0], na.rm = TRUE),
    Paintor_cs_median = median(sum_Paintor_cs[sum_Paintor_cs != 0], na.rm = TRUE),
    SuSiEx_cs_median = median(sum_SuSiEx_cs[sum_SuSiEx_cs != 0], na.rm = TRUE),
    XMAP_cs_median = median(sum_XMAP_cs[sum_XMAP_cs != 0], na.rm = TRUE),
    CARMAX_cs_median = median(sum_CARMAX_cs[sum_CARMAX_cs != 0], na.rm = TRUE),
    SuSiE_metal_cs_median = median(sum_SuSiE_metal_cs[sum_SuSiE_metal_cs != 0], na.rm = TRUE)
    #MESuSiE_res_mlk_annot_20_baseline_cs_median = median(sum_MESuSiE_res_mlk_annot_20_baseline_cs[sum_MESuSiE_res_mlk_annot_20_baseline_cs != 0], na.rm = TRUE),
    #MESuSiE_res_glmnet_annot_20_baseline_cs_median = median(sum_MESuSiE_res_glmnet_annot_20_baseline_cs[sum_MESuSiE_res_glmnet_annot_20_baseline_cs != 0], na.rm = TRUE)
  )
  
#MESuSiE identified set overlaps with MFD, SuSiE, MESuSiE credible set
res_all%>%group_by(PHENONAME)%>%summarise(overlap_prop = sum(CS==1&(SuSiE_cs==1|MESuSiE_cs==1))/sum(CS==1))
res_all%>%summarise(overlap_prop = sum(MESuSiE_cs==1&(SuSiE_cs==1|Paintor_cs==1))/sum(MESuSiE_cs==1))



# Median |Z-scores| of SNPs within credible set to quantile the strength of detected SNPs
MFD_median_z<-res_all%>%filter(CS==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr),na.rm=T)))%>%pull(zmax)
MESuSiE_median_z<-res_all%>%filter(MESuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr))))%>%pull(zmax)
SuSiE_median_z<-res_all%>%filter(SuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr),na.rm=T)))%>%pull(zmax)
Paintor_median_z<-res_all%>%filter(Paintor_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr))))%>%pull(zmax)
SuSiEx_median_z<-res_all%>%filter(SuSiEx_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr))))%>%pull(zmax)
XMAP_median_z<-res_all%>%filter(XMAP_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr))))%>%pull(zmax)
CARMAX_median_z<-res_all%>%filter(CARMAX_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr))))%>%pull(zmax)
SuSiE_metal_median_z<-res_all%>%filter(SuSiE_metal_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr),na.rm=T)))%>%pull(zmax)



#Paintor_20_median_z<-res_all%>%filter(Paintor_cs_sub_20==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB))))%>%pull(zmax)


cat(round(c(MFD_median_z,MESuSiE_median_z,SuSiE_median_z,Paintor_median_z,SuSiEx_median_z,XMAP_median_z,CARMAX_median_z,SuSiE_metal_median_z),2))
# Median 5.55 5.52 5.49 5.54 5.47 6.04 2.51 2.45
# eQTL enrichment by PHENONAMEf
res_all<-res_all%>%mutate(utr_comb = ifelse((utr_3+utr_5)>0,1,0))
#res_all<-res_all%>%mutate(utr_comb_flanking = ifelse((UTR_3_UCSC.flanking.500_common+UTR_5_UCSC.flanking.500_common)>0,1,0))

#ann_col_name<-c("missense", "synonymous", "utr_comb", "promotor", "CRE","brain_ind_eQTL","non_synonymous_common","Conserved_LindbladToh_common","GTEx_eQTL_MaxCPP_common","H3K4me3_Trynka_common","utr_comb_flanking","BLUEPRINT_H3K27acQTL_MaxCPP_common","BLUEPRINT_H3K4me1QTL_MaxCPP_common","BLUEPRINT_DNA_methylation_MaxCPP_common")
ann_col_name<-c( "non_synonymous_common","synonymous", "promotor","utr_comb","CRE","heart_ind_eQTL","artery_ind_eQTL","brain_ind_eQTL")

# Functions for calculating fold enrichment
calc_fold_enrichment <- function(df, cs_col, ann_col_name) {
  df %>%
    group_by(Region) %>%
    filter(sum(!!sym(cs_col)) != 0) %>%
    group_by(PHENONAME, !!sym(cs_col)) %>%
    summarise(across(ann_col_name, ~ sum(.x, na.rm = TRUE) / n())) %>%
    group_by(PHENONAME) %>%
    summarise(across(ann_col_name, ~ .x[!!sym(cs_col) == 1] / .x[!!sym(cs_col) == 0]))
}

MFD_PIP_ann <- calc_fold_enrichment(res_all, "CS", ann_col_name)
MESuSiE_PIP_ann <- calc_fold_enrichment(res_all, "MESuSiE_cs", ann_col_name)
SuSiE_PIP_ann <- calc_fold_enrichment(res_all, "SuSiE_cs", ann_col_name)
Paintor_PIP_ann <- calc_fold_enrichment(res_all, "Paintor_cs", ann_col_name)
SuSiEx_PIP_ann<-calc_fold_enrichment(res_all, "SuSiEx_cs", ann_col_name)
XMAP_PIP_ann<-calc_fold_enrichment(res_all, "XMAP_cs", ann_col_name)
CARMAX_PIP_ann<-calc_fold_enrichment(res_all, "CARMAX_cs", ann_col_name)
SuSiE_metal_PIP_ann<-calc_fold_enrichment(res_all, "SuSiE_metal_cs", ann_col_name)



# Combine results
PHENONAME_CS_enrichment <- bind_rows(
  MFD_PIP_ann %>% mutate(Method = "MFD"),
  MESuSiE_PIP_ann %>% mutate(Method = "MESuSiE"),
  SuSiE_PIP_ann %>% mutate(Method = "SuSiE"),
  Paintor_PIP_ann %>% mutate(Method = "Paintor"),
  SuSiEx_PIP_ann %>% mutate(Method = "SuSiEx"),
  XMAP_PIP_ann %>% mutate(Method = "XMAP"),
  CARMAX_PIP_ann %>% mutate(Method = "CARMAx"),
  SuSiE_metal_PIP_ann %>% mutate(Method = "SuSiE_Meta")
) %>% mutate(Method = factor(Method, levels = c("MFD","MESuSiE","SuSiE", "Paintor","SuSiEx","XMAP","CARMAx","SuSiE_Meta")))%>%
  dplyr::select(PHENONAME,heart_ind_eQTL,artery_ind_eQTL,brain_ind_eQTL ,Method )#%>%dplyr::rename(eQTL = brain_ind_eQTL)
# Pivot to long format
PHENONAME_CS_enrichment_long <- PHENONAME_CS_enrichment %>%
  pivot_longer(cols = -c(Method, PHENONAME), names_to = "Cat", values_to = "Prop") %>%
  mutate(Method = factor(Method, levels = c("MFD","MESuSiE","SuSiE", "Paintor","SuSiEx","XMAP","CARMAx","SuSiE_Meta")))

PHENONAME_CS_enrichment_long%>%group_by(Method)%>%summarise(min(Prop),max(Prop),mean(Prop))

# Method         `min(Prop)` `max(Prop)` `mean(Prop)`
# <fct>                <dbl>       <dbl>        <dbl>
# 1 MESuSiE              4.18        11.1          8.33
# 2 MESuSiE_mlk_10       4.15        20.7          9.47
# 3 MESuSiE_glm_10       4.57        16.7          9.45
# 4 MESuSiE_mlk          3.64        26.7         12.3 
# 5 MESuSiE_glm          3.06        22.0         10.3 
# 6 SuSiE                3.10         8.64         7.00
# 7 Paintor              0.958        3.83         2.32
# 8 Paintor_annot        1.66         2.92         2.22

# Enrichment of 95% credible set without by PHENONAME
calc_fold_enrichment_marginal<-function(df, cs_col, ann_col_name) {
  df %>%group_by(Region) %>%
    filter(sum(!!sym(cs_col)) != 0) %>%
    group_by( !!sym(cs_col)) %>%
    summarise(across(ann_col_name, ~ sum(.x, na.rm = TRUE) / n())) %>%
    summarise(across(ann_col_name, ~ .x[!!sym(cs_col) == 1] / .x[!!sym(cs_col) == 0]))
}


MFD_PIP_ann <- calc_fold_enrichment_marginal(res_all, "CS", ann_col_name)
MESuSiE_PIP_ann <- calc_fold_enrichment_marginal(res_all, "MESuSiE_cs", ann_col_name)
SuSiE_PIP_ann <- calc_fold_enrichment_marginal(res_all, "SuSiE_cs", ann_col_name)
Paintor_PIP_ann <- calc_fold_enrichment_marginal(res_all, "Paintor_cs", ann_col_name)
SuSiEx_PIP_ann<-calc_fold_enrichment_marginal(res_all, "SuSiEx_cs", ann_col_name)
XMAP_PIP_ann<-calc_fold_enrichment_marginal(res_all, "XMAP_cs", ann_col_name)
CARMAX_PIP_ann<-calc_fold_enrichment_marginal(res_all, "CARMAX_cs", ann_col_name)
SuSiE_metal_PIP_ann<-calc_fold_enrichment_marginal(res_all, "SuSiE_metal_cs", ann_col_name)



CS_enrichment <- bind_rows(
  MFD_PIP_ann %>% mutate(Method = "MFD"),
  MESuSiE_PIP_ann %>% mutate(Method = "MESuSiE"),
  SuSiE_PIP_ann %>% mutate(Method = "SuSiE"),
  Paintor_PIP_ann %>% mutate(Method = "Paintor"),
  SuSiEx_PIP_ann %>% mutate(Method = "SuSiEx"),
  XMAP_PIP_ann %>% mutate(Method = "XMAP"),
  CARMAX_PIP_ann %>% mutate(Method = "CARMAx"),
  SuSiE_metal_PIP_ann %>% mutate(Method = "SuSiE_Meta")
) %>% mutate(Method = factor(Method, levels = c("MFD","MESuSiE","SuSiE", "Paintor","SuSiEx","XMAP","CARMAx","SuSiE_Meta")))%>% 
  dplyr::rename(`Non-synonymous` = non_synonymous_common ,Synonymous = synonymous,UTR = utr_comb,Promotor = promotor,Brain_eQTL = brain_ind_eQTL,Heart_eQTL = heart_ind_eQTL,Artery_eQTL = artery_ind_eQTL)


# Pivot to long format
CS_enrichment_long <- CS_enrichment %>%
  pivot_longer(cols = -c(Method), names_to = "Cat", values_to = "Prop") %>%
  mutate(Method = factor(Method, levels = c("MFD","MESuSiE","SuSiE","SuSiE_Meta","Paintor","SuSiEx","XMAP","CARMAx"))) %>%
  mutate(Cat = factor(Cat, levels = c("Non-synonymous", "Synonymous", "UTR", "Promotor", "CRE","Brain_eQTL","Heart_eQTL","Artery_eQTL")))%>%
  mutate(Prop = round(Prop, 2))




# Enrichment of top 100 signal

res_all<-res_all%>%mutate(MFD_Signal = case_when(
  PIP_Ancestry_1>0.5~1,
  PIP_Ancestry_2>0.5~2,
  PIP_Shared>0.5~3,
  .default =0))
res_all<-res_all%>%mutate(SuSiE_Signal = case_when(
  SuSiE_PIP_EU>0.5&SuSiE_PIP_BB >0.5~1,
  SuSiE_PIP_EU>0.5&SuSiE_PIP_BB<0.5~2,
  SuSiE_PIP_EU<0.5&SuSiE_PIP_BB>0.5~3,
  .default =0))
res_all<-res_all%>%mutate(MESuSiE_Signal = case_when(
  MESuSiE_PIP_Shared>0.5~1,
  MESuSiE_PIP_WB>0.5~2,
  MESuSiE_PIP_BB>0.5~3,
  .default =0))
res_all<-res_all%>%mutate(Paintor_Signal = ifelse(Paintor_PIP_Either>0.5,1,0))
res_all<-res_all%>%mutate(SuSiEx_Signal = ifelse(SuSiEx_PIP_Either>0.5,1,0))
res_all<-res_all%>%mutate(XMAP_Signal = ifelse(XMAP_PIP_Either>0.5,1,0))
res_all<-res_all%>%mutate(CARMAX_Signal = case_when(
  CARMAX_PIP_Shared>0.5~1,
  CARMAX_PIP_WB>0.5~2,
  CARMAX_PIP_BB>0.5~3,
  .default =0))
res_all<-res_all%>%mutate(SuSiE_merged_Signal = ifelse(SuSiE_merged_PIP_Either>0.5,1,0))



top_N_signal = res_all %>%filter(MESuSiE_Signal != 0, MESuSiE_PIP_Either > 0.5) %>%nrow()
top_N_signal = res_all %>%filter(MESuSiE_PIP_Either_all_mlk > 0.5) %>%nrow()

bg_an<-res_all%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/(n()-top_N_signal)))
MESuSiE_Signal_ann<-res_all%>%filter(MESuSiE_Signal!=0)%>% arrange(desc(MESuSiE_PIP_Either))%>%top_n(n = top_N_signal, wt = MESuSiE_PIP_Either)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an
#MESuSiE_MLK_10_Signal_ann<-res_all%>%filter(MESuSiE_MLK_10_Signal!=0)%>% arrange(desc(MESuSiE_PIP_Either_10_sub_mlk))%>%top_n(n = top_N_signal, wt = MESuSiE_PIP_Either_10_sub_mlk)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an
#MESuSiE_GLM_10_Signal_ann<-res_all%>%filter(MESuSiE_GLM_10_Signal!=0)%>% arrange(desc(MESuSiE_PIP_Either_10_sub_glmnet))%>%top_n(n = top_N_signal, wt = MESuSiE_PIP_Either_10_sub_glmnet)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an
MESuSiE_MLK_Signal_ann<-res_all%>%filter(MESuSiE_MLK_Signal!=0)%>% arrange(desc(MESuSiE_PIP_Either_all_mlk))%>%top_n(n = top_N_signal, wt = MESuSiE_PIP_Either_all_mlk)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an
MESuSiE_GLM_Signal_ann<-res_all%>%filter(MESuSiE_GLM_Signal!=0)%>% arrange(desc(MESuSiE_PIP_Either_all_glmnet))%>%top_n(n = top_N_signal, wt = MESuSiE_PIP_Either_all_glmnet)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an
SuSiE_Signal_ann<-res_all%>%filter(SuSiE_Signal!=0)%>% arrange(desc(SuSiE_PIP))%>%top_n(n = top_N_signal, wt = SuSiE_PIP)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an
Paintor_Signal_ann<-res_all%>%filter(Paintor_Signal!=0) %>% arrange(desc(Paintor_PIP))%>%top_n(n = top_N_signal, wt = Paintor_PIP)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an
#Paintor_Signal_10_ann<-res_all%>%filter(Paintor_10_Signal!=0) %>% arrange(desc(Paintor_PIP_sub_baseline))%>%top_n(n = top_N_signal, wt = Paintor_PIP_sub_baseline)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an
Paintor_Signal_all_ann<-res_all%>%filter(Paintor_all_Signal!=0) %>% arrange(desc(Paintor_PIP_all))%>%top_n(n = top_N_signal, wt = Paintor_PIP_all)%>%summarise(across(ann_col_name,~ sum(.x, na.rm = TRUE)/n()))/bg_an






Signal_enrichment <- bind_rows(
  MESuSiE_Signal_ann %>% mutate(Method = "MESuSiE"),
  #MESuSiE_MLK_10_Signal_ann %>% mutate(Method = "MESuSiE_mlk_10"),
  #MESuSiE_GLM_10_Signal_ann %>% mutate(Method = "MESuSiE_glm_10"),
  MESuSiE_MLK_Signal_ann %>% mutate(Method = "MESuSiE_mlk"),
  MESuSiE_GLM_Signal_ann %>% mutate(Method = "MESuSiE_glm"),
  SuSiE_Signal_ann %>% mutate(Method = "SuSiE"),
  Paintor_Signal_ann %>% mutate(Method = "Paintor"),
  #Paintor_Signal_10_ann %>% mutate(Method = "Paintor_annot")
  Paintor_Signal_all_ann %>% mutate(Method = "Paintor_annot_all"),
) %>% mutate(Method = factor(Method, levels = c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE","Paintor","Paintor_annot_all")))%>% 
  dplyr::rename(Missense = missense ,Synonymous = synonymous,UTR = utr_comb,Promotor = promotor,eQTL = brain_ind_eQTL)

# Pivot to long format
Signal_enrichment_long <- Signal_enrichment %>%
  pivot_longer(cols = -c(Method), names_to = "Cat", values_to = "Prop") %>%
  mutate(Method = factor(Method, levels = c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE", "Paintor","Paintor_annot_all"))) %>%
  mutate(Cat = factor(Cat, levels = c("Missense", "Synonymous", "UTR", "Promotor", "CRE","eQTL")))%>%
  mutate(Prop = round(Prop, 2))




# Number of Signal with PIP > 0.5

##Number of signals
res_all%>%group_by(MESuSiE_Signal)%>%summarise(N_signl = n())
#res_all%>%group_by(MESuSiE_MLK_10_Signal)%>%summarise(N_signl = n())
#res_all%>%group_by(MESuSiE_GLM_10_Signal)%>%summarise(N_signl = n())
res_all%>%group_by(SuSiE_Signal)%>%summarise(N_signl = n())
res_all%>%group_by(MFD_Signal)%>%summarise(N_signl = n())
res_all%>%group_by(CARMAX_Signal)%>%summarise(N_signl = n())
#res_all%>%group_by(Paintor_Signal)%>%summarise(N_signl = n())
#res_all%>%group_by(Paintor_10_Signal)%>%summarise(N_signl = n())


res_all%>%group_by(PHENONAME,MESuSiE_Signal)%>%summarise(N_signl = n())
res_all%>%group_by(PHENONAME,SuSiE_Signal)%>%summarise(N_signl = n())
res_all%>%group_by(PHENONAME,MFD_Signal)%>%summarise(N_signl = n())
res_all%>%group_by(PHENONAME,CARMAX_Signal)%>%summarise(N_signl = n())
#res_all%>%group_by(PHENONAME,MESuSiE_GLM_Signal)%>%summarise(N_signl = n())
#res_all%>%group_by(PHENONAME,SuSiE_Signal)%>%summarise(N_signl = n())
#res_all%>%group_by(PHENONAME,Paintor_Signal)%>%summarise(N_signl = n())
#res_all%>%group_by(PHENONAME,Paintor_10_Signal)%>%summarise(N_signl = n())

# Proportion of shared signal
res_all%>%filter(MESuSiE_Signal!=0)%>%group_by(MESuSiE_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
#res_all%>%filter(MESuSiE_MLK_10_Signal!=0)%>%group_by(MESuSiE_MLK_10_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
#res_all%>%filter(MESuSiE_GLM_10_Signal!=0)%>%group_by(MESuSiE_GLM_10_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
res_all%>%filter(MESuSiE_MLK_Signal!=0)%>%group_by(MESuSiE_MLK_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
res_all%>%filter(MESuSiE_GLM_Signal!=0)%>%group_by(MESuSiE_GLM_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
res_all%>%filter(SuSiE_Signal!=0)%>%group_by(SuSiE_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))

res_all%>%filter(MESuSiE_Signal!=0)%>%group_by(PHENONAME,MESuSiE_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
res_all%>%filter(MESuSiE_MLK_10_Signal!=0)%>%group_by(PHENONAME,MESuSiE_MLK_10_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
res_all%>%filter(MESuSiE_GLM_10_Signal!=0)%>%group_by(PHENONAME,MESuSiE_GLM_10_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
res_all%>%filter(MESuSiE_MLK_Signal!=0)%>%group_by(PHENONAME,MESuSiE_MLK_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
res_all%>%filter(MESuSiE_GLM_Signal!=0)%>%group_by(PHENONAME,MESuSiE_GLM_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))
res_all%>%filter(SuSiE_Signal!=0)%>%group_by(PHENONAME,SuSiE_Signal)%>%summarise(n = n())%>%summarise(freq = n / sum(n))

###Proportion of shared signal lowering threshold for SuSiE
res_all%>%group_by(PHENONAME)%>%summarise(SuSiE_Shared_n = sum(SuSiE_Shared>0.1),SuSiE_WB_n = sum(susie_EU>0.1&susie_BB<0.1),SuSiE_BB_n = sum(susie_EU<0.1&susie_BB>0.1))%>%summarise(SuSiE_Shared_n/(SuSiE_Shared_n+SuSiE_WB_n+SuSiE_BB_n))
res_all%>%summarise(SuSiE_Shared_n = sum(SuSiE_Shared>0.1),SuSiE_WB_n = sum(susie_EU>0.1&susie_BB<0.1),SuSiE_BB_n = sum(susie_EU<0.1&susie_BB>0.1))%>%summarise(SuSiE_Shared_n/(SuSiE_Shared_n+SuSiE_WB_n+SuSiE_BB_n))




################Signal number by PHENONAME, not running for MESuSiE + Annotation#######################
Signal_number<- res_all%>%group_by(PHENONAME)%>%
  summarise(Paintor_Either_n = sum(Paintor_Signal!=0),
            SuSiE_Shared_n = sum(SuSiE_Signal==1),
            SuSiE_WB_n = sum(SuSiE_Signal==2),
            SuSiE_BB_n = sum(SuSiE_Signal==3),
            MESuSiE_Shared_n = sum(MESuSiE_Signal==1),
            MESuSiE_WB_n = sum(MESuSiE_Signal==2),
            MESuSiE_BB_n = sum(MESuSiE_Signal==3))
Signal_number<-Signal_number%>%
  pivot_longer(cols = -c(PHENONAME), names_to = "Cat", values_to = "Num")%>%
  separate(Cat, into = c("Method", "Signal"), sep = "_", extra = "merge") %>%
  mutate(
    Method = case_when(
      str_detect(Method, "MESuSiE") ~ "MESuSiE",
      str_detect(Method, "SuSiE") ~ "SuSiE",
      str_detect(Method, "Paintor") ~ "Paintor",
      TRUE ~ Method
    ),
    Signal = case_when(
      str_detect(Signal, "BB_n") ~ "AFR",
      str_detect(Signal, "WB_n") ~ "EUR",
      str_detect(Signal, "Shared_n") ~ "Shared",
      str_detect(Signal, "Either") ~ "Either",
      TRUE ~ Signal
    )
  )
data.frame(Signal_number)

# Correlation of Shared and ancestry specific signal
# Ancestry-specific signal MESuSiE unique/SuSiE unique/Both
cor_ancestry_MESuSiE <- res_all %>%
  filter(MESuSiE_Signal%in%c(2,3)) %>% filter(!(SuSiE_Signal%in%c(2,3)))%>%
  mutate(Method = "MESuSiE") %>%
  dplyr::select(Beta_WB, Beta_BB, PHENONAME, Method)
cor_ancestry_SuSiE <- res_all %>%
  filter(SuSiE_Signal%in%c(2,3)) %>%filter(!(MESuSiE_Signal%in%c(2,3)))%>%
  mutate(Method = "SuSiE") %>%
  dplyr::select(Beta_WB, Beta_BB, PHENONAME, Method)
cor_ancestry_both <- res_all %>%filter(MESuSiE_Signal%in%c(2,3),SuSiE_Signal%in%c(2,3))%>%
  mutate(Method = "Both") %>%
  dplyr::select(Beta_WB, Beta_BB, PHENONAME, Method)
cor_ancestry <- rbind(cor_ancestry_MESuSiE, cor_ancestry_SuSiE,cor_ancestry_both)
cor_ancestry<-cor_ancestry%>%mutate(Method = factor(Method,levels = c("MESuSiE","SuSiE","Both")))	
###Ancestry specific signal marginal correlations
cor_ancestry %>%
  group_by(Method, PHENONAME) %>%
  summarise(cor(Beta_WB, Beta_BB))
res_all %>%filter(MESuSiE_Signal%in%c(2,3))%>%
  summarise(cor(Beta_WB, Beta_BB))
cor_ancestry  %>%
  group_by(Method) %>%
  summarise(cor(Beta_WB, Beta_BB))
###Shared signal marginal correlations
res_all %>%filter(MESuSiE_Signal==1) %>%group_by(PHENONAME) %>%summarise(cor(Beta_WB, Beta_BB))
res_all %>%filter(MESuSiE_Signal==1) %>%summarise(cor(Beta_WB, Beta_BB))
#	Test Correlation Difference ancestry vs shared
library(boot)
function_cor <- function(data, i){
  d2 <- data[i,] 
  return(cor(d2$Beta_WB, d2$Beta_BB))
}
set.seed(1128)
ancestry_signal_data<-res_all%>%filter(MESuSiE_Signal%in%c(2,3))
bootstrap_correlation_ancestry <- boot(ancestry_signal_data,function_cor,R=1000)
shared_signal_data<-res_all%>%filter(MESuSiE_Signal==1)
bootstrap_correlation_shared <- boot(shared_signal_data,function_cor,R=1000)
res_all%>%filter(MESuSiE_Signal==1)%>%group_by(PHENONAME)%>%summarise(margin_cor = cor(Beta_WB,Beta_BB))
res_all%>%filter(MESuSiE_Signal==1)%>%summarise(margin_cor = cor(Beta_WB,Beta_BB))
cor_dif<-mean(bootstrap_correlation_shared$t-bootstrap_correlation_ancestry$t)
cor_dif_se<-sqrt(var(bootstrap_correlation_shared$t-bootstrap_correlation_ancestry$t))
cor_dif_P<-2*pnorm(-abs(cor_dif/cor_dif_se))

# MAF Diff
maf_shared<-data.frame(res_all)%>%filter(MESuSiE_Signal==1)%>%mutate(MAF_diff = MAF_WB-MAF_BB)%>%pull(MAF_diff)
maf_ancestry<-data.frame(res_all)%>%filter(MESuSiE_Signal%in%c(2,3))%>%mutate(MAF_diff = MAF_WB-MAF_BB)%>%pull(MAF_diff)
wilcox.test(maf_ancestry, maf_shared, alternative = "two.sided")	
# Wilcoxon rank sum test with continuity correction
# 
# data:  maf_ancestry and maf_shared
# W = 1079, p-value = 0.1958
# alternative hypothesis: true location shift is not equal to 0
median(maf_shared) # -0.0367385
median(maf_ancestry) # -0.06986  



maf_shared<-data.frame(res_all)%>%filter(MESuSiE_MLK_Signal == 1)%>%mutate(MAF_diff = MAF_WB-MAF_BB)%>%pull(MAF_diff)
maf_ancestry<-data.frame(res_all)%>%filter(MESuSiE_MLK_Signal %in%c(2,3))%>%mutate(MAF_diff = MAF_WB-MAF_BB)%>%pull(MAF_diff)
wilcox.test(maf_ancestry, maf_shared, alternative = "two.sided")	
# Wilcoxon rank sum test with continuity correction
# 
# data:  maf_ancestry and maf_shared
# W = 2339, p-value = 0.09448
# alternative hypothesis: true location shift is not equal to 0
median(maf_shared) # -0.0245965
median(maf_ancestry) # -0.047286  

# maf_shared<-data.frame(res_all)%>%filter(MESuSiE_GLM_Signal==1)%>%mutate(MAF_diff = MAF_WB-MAF_BB)%>%pull(MAF_diff)
# maf_ancestry<-data.frame(res_all)%>%filter(MESuSiE_GLM_Signal%in%c(2,3))%>%mutate(MAF_diff = MAF_WB-MAF_BB)%>%pull(MAF_diff)
# wilcox.test(maf_ancestry, maf_shared, alternative = "two.sided")	
# Wilcoxon rank sum test with continuity correction
# 
# data:  maf_ancestry and maf_shared
# W = 2338, p-value = 0.4879
# alternative hypothesis: true location shift is not equal to 0
# median(maf_shared) # -0.031028
# median(maf_ancestry) # -0.0326425




#####################################################################
#
#     Plot directory and theme
#
####################################################################
plot_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/Figure/"
system(paste0("mkdir -p ",plot_dir))
custom_theme <- function() {
  theme(
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5),  
    axis.title.x = element_text(size = 7, face="bold"),
    axis.title.y = element_text(size = 7, face="bold"),
    strip.text.x = element_text(size = 5),
    strip.text.y = element_text(size = 5),
    strip.background = element_blank(),
    legend.text = element_text(size=7),
    legend.title = element_text(size=7, face="bold"),
    plot.title = element_text(size=7, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line = element_line(color = "black")
  )
}
################################################
#
#		Set Size/Z-score/eQTL 
#
#
###############################################
################################################
#
#		Set SiZe Part
#
###############################################		
###Median set size by PHENONAME
all_sets_info<-data.frame(res_all%>%group_by(PHENONAME,Region) %>% summarise(across(c("CS","MESuSiE_cs", "SuSiE_cs","Paintor_cs","SuSiEx_cs","XMAP_cs","CARMAX_cs","SuSiE_metal_cs"), ~ sum(.x, na.rm = TRUE))))%>%filter(CS!= 0,MESuSiE_cs!=0,SuSiE_cs!=0,Paintor_cs!=0,SuSiEx_cs !=0,XMAP_cs!=0,CARMAX_cs!=0) #Median Set Size across all locus
all_sets_info$SuSiE_metal_cs = 1
all_sets_info_long<-all_sets_info%>%pivot_longer(!(PHENONAME|Region), names_to = "Method", values_to = "Count")
all_sets_info_long$Method<-factor(all_sets_info_long$Method,levels=c("CS","MESuSiE_cs", "SuSiE_cs","Paintor_cs","SuSiEx_cs","XMAP_cs","CARMAX_cs","SuSiE_metal_cs"))
levels(all_sets_info_long$Method)<-c("MFD","MESuSiE","SuSiE", "Paintor","SuSiEx","XMAP","CARMAx","SuSiE_Meta")
all_sets_info_long$Count = log(all_sets_info_long$Count,base = 2)+1


# all_sets_info_long_save <- all_sets_info_long %>%
#   mutate(causal = num_causal,
#          h2 = h2_num)


p_set = ggplot(data =all_sets_info_long,aes(x = PHENONAME, y=Count,fill=Method))+geom_boxplot(aes(x = PHENONAME,fill=Method),outlier.size = 0.1,fatten = 0.5,color = "darkgray")+scale_fill_manual(values=c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#B2D3A4","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e"),guide=FALSE)
p_set =p_set + theme_bw() + xlab("") +ylab("log2(Set Size) + 1")+coord_cartesian(ylim=c(0,12))
p_set= p_set+custom_theme()

################################################
#
#		Z-score Part
#
###############################################		
MFD_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(CS==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr),na.rm=T)))
MESuSiE_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(MESuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr))))%>%pull(zmax)
SuSiE_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(SuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr),na.rm=T)))%>%pull(zmax)
Paintor_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(Paintor_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr))))%>%pull(zmax)
SuSiEx_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(SuSiEx_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr))))%>%pull(zmax)
XMAP_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(XMAP_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr))))%>%pull(zmax)
CARMAX_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(CARMAX_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr))))%>%pull(zmax)
SuSiE_metal_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(SuSiE_metal_cs==1)%>%summarise(zmax = median(pmax(abs(Z_eur),abs(Z_afr),na.rm=T)))%>%pull(zmax)

set_size_z_info<-data.frame(cbind(MFD_cs_Z,MESuSiE_cs_Z,SuSiE_cs_Z,Paintor_cs_Z,SuSiEx_cs_Z,XMAP_cs_Z,CARMAX_cs_Z,SuSiE_metal_cs_Z))
colnames(set_size_z_info)<-c("PHENONAME",c("MFD","MESuSiE","SuSiE", "Paintor","SuSiEx","XMAP","CARMAx","SuSiE_Meta"))
set_size_z_info_long<-set_size_z_info %>%pivot_longer(!(PHENONAME), names_to = "Method", values_to = "Z")%>%mutate(Method = factor(Method, levels=c("MFD","MESuSiE","SuSiE", "Paintor","SuSiEx","XMAP","CARMAx","SuSiE_Meta")))

p_z = ggplot(data = set_size_z_info_long,aes(x = PHENONAME, y=Z,fill=Method))+geom_bar( stat = "identity",position="dodge")+scale_fill_manual(values=c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#B2D3A4","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e"))
p_z = p_z + geom_text(label = round(set_size_z_info_long$Z,2),position = position_dodge(width = 1),vjust=-0.5,size = 5*5/14)
p_z = p_z + theme_bw() + xlab("") +ylab("Median |Z|")+ ylim(0,max(round(set_size_z_info_long$Z,2)+1))
p_z = p_z +custom_theme()
################################################
#
#		eQTL enrichment 
#
#
###############################################	

PHENONAME_CS_enrichment_long_sub = PHENONAME_CS_enrichment_long %>% filter(Cat == "brain_ind_eQTL")
p_eQTL_brain <- ggplot(PHENONAME_CS_enrichment_long_sub, aes(x = PHENONAME, y = Prop, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +scale_fill_manual(values =c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#B2D3A4","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e")) +
  geom_text(,label = round(PHENONAME_CS_enrichment_long_sub$Prop,2),position = position_dodge(width = 1),vjust=-0.5,size = 5*5/14)+
  xlab("") + ylab("Brain eQTL Fold Enrichment") + ylim(0,max(round(PHENONAME_CS_enrichment_long_sub$Prop))+1)+
  theme_bw() + custom_theme()

PHENONAME_CS_enrichment_long_sub = PHENONAME_CS_enrichment_long %>% filter(Cat == "heart_ind_eQTL")
p_eQTL_heart <- ggplot(PHENONAME_CS_enrichment_long_sub, aes(x = PHENONAME, y = Prop, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +scale_fill_manual(values =c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#B2D3A4","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e")) +
  geom_text(,label = round(PHENONAME_CS_enrichment_long_sub$Prop,2),position = position_dodge(width = 1),vjust=-0.5,size = 5*5/14)+
  xlab("") + ylab("Heart eQTL Fold Enrichment") + ylim(0,max(round(PHENONAME_CS_enrichment_long_sub$Prop))+1)+
  theme_bw() + custom_theme()

PHENONAME_CS_enrichment_long_sub = PHENONAME_CS_enrichment_long %>% filter(Cat == "artery_ind_eQTL")
p_eQTL_artery <- ggplot(PHENONAME_CS_enrichment_long_sub, aes(x = PHENONAME, y = Prop, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +scale_fill_manual(values =c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#B2D3A4","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e")) +
  geom_text(,label = round(PHENONAME_CS_enrichment_long_sub$Prop,2),position = position_dodge(width = 1),vjust=-0.5,size = 5*5/14)+
  xlab("") + ylab("Artery eQTL Fold Enrichment") + ylim(0,max(round(PHENONAME_CS_enrichment_long_sub$Prop))+1)+
  theme_bw() + custom_theme()



p_out<-p_set/p_z/p_eQTL_brain/p_eQTL_heart/p_eQTL_artery+plot_annotation(tag_levels = 'a')+plot_layout(guides = "collect",heights = c(1.5,1,1,1,1))&theme(legend.position = 'bottom',plot.tag = element_text(size = 6,face="bold"))
p_out
ggsave(paste0(plot_dir,"set_size_eqtl_all_traits.pdf"),p_out,width=140,height = 280,dpi=500,units='mm')
################################################################################
#
#
#     Functional Annotation enrichment for 95% credible set SNPS
#
#
################################################################################
# Use one shared dodge spec so bars + labels align perfectly
dodge <- position_dodge2(width = 0.9, preserve = "single")

p_signal <- ggplot(CS_enrichment_long, aes(x = Cat, y = Prop, fill = Method)) +
  geom_col(position = dodge, width = 0.8) +
  geom_text(
    aes(label = sprintf("%.1f", Prop)),
    position = dodge,
    vjust = -0.35,
    size = 5/14*5
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(
    values = c(
      "MESuSiE"     = "#8da0cb",
      "SuSiE"       = "#66c2a5",
      "MFD"         = "#377eb8",
      "SuSiE_Meta"  = "#9FBA95",
      "Paintor"     = "#fc8d62",
      "MultiSuSiE"  = "#e78ac3",
      "SuSiEx"      = "#E89DA0",
      "XMAP"        = "#ffd92f",
      "CARMAx"      = "#f2b56e"
    )
  ) +
  labs(x = NULL, y = "Fold Enrichment Credible Set") +
  scale_y_continuous(
    limits = c(0, ceiling(max(CS_enrichment_long$Prop, na.rm = TRUE)) + 1),
    expand = expansion(mult = c(0, 0.08))  # add headroom for labels
  ) +
  theme_bw() +
  custom_theme() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))


ggsave(paste0(plot_dir,"Supp_Enrichment_Comparison_Set_CS_updated.pdf"),p_signal,dpi=600,width=210,height =120,units = 'mm')


p_out<-p_set/p_signal+plot_annotation(tag_levels = 'a')+plot_layout(guides = "collect",heights = c(1,1))&theme(legend.position = 'bottom',plot.tag = element_text(size = 6,face="bold"))
ggsave(paste0(plot_dir,"Enrichment_Comparison_Set_CS.pdf"),p_out,dpi=500,width=280,height =120,units = 'mm')






################################################################################
#
#
#     Functional Annotation enrichment for top 100 PIP SNPs
#
#
################################################################################	
p_signal <- ggplot(data = Signal_enrichment_long,aes(x = Cat, y = Prop, fill = Method)) +
  geom_col(position = "dodge") + scale_fill_manual(values = c("MESuSiE"="#023e8a","MESuSiE_mlk" = "#0585e6","MESuSiE_glm" = "#05a9ff","SuSiE"="#2a9d8f","Paintor"="#f4a261","Paintor_annot_all"="#d17f45")) +
  geom_text(aes(x=Cat,group=Method,y=Prop,label=Prop),position = position_dodge(width = 1),vjust=-0.5,size = 5/14*5) + 
  geom_hline(yintercept = 1, linetype = "dashed") + 
  xlab("") + ylab("Fold Enrichment Top Signal") +ylim(0,round(max(Signal_enrichment_long$Prop))+1) +
  theme_bw() + custom_theme()
p_out<-p_set/p_signal+plot_annotation(tag_levels = 'a')+plot_layout(guides = "collect",heights = c(1,1))&theme(legend.position = 'bottom',plot.tag = element_text(size = 7,face="bold"))
ggsave(paste0(plot_dir,"Enrichment_Comparison_Set_Top.pdf"),p_out,dpi=500,width=180,height =120,units = 'mm')


############################################################################
#
#
#              Proportion of Signal Plot
#
#
############################################################################

Signal_number<- res_all%>%group_by(PHENONAME)%>%
  summarise(Paintor_Either_n = sum(Paintor_Signal!=0),
            #Paintor10_Either_n = sum(Paintor_10_Signal!=0),
            SuSiE_Shared_n = sum(SuSiE_Signal==1),
            SuSiE_WB_n = sum(SuSiE_Signal==2),
            SuSiE_BB_n = sum(SuSiE_Signal==3),
            MESuSiE_Shared_n = sum(MESuSiE_Signal==1),
            MESuSiE_WB_n = sum(MESuSiE_Signal==2),
            MESuSiE_BB_n = sum(MESuSiE_Signal==3),
            #MESuSiEmlk10_Shared_n = sum(MESuSiE_MLK_10_Signal==1),
            #MESuSiEmlk10_WB_n = sum(MESuSiE_MLK_10_Signal==2),
            #MESuSiEmlk10_BB_n = sum(MESuSiE_MLK_10_Signal==3),
            #MESuSiEglm10_Shared_n = sum(MESuSiE_GLM_10_Signal==1),
            #MESuSiEglm10_WB_n = sum(MESuSiE_GLM_10_Signal==2),
            #MESuSiEglm10_BB_n = sum(MESuSiE_GLM_10_Signal==3),
            MESuSiEmlk_Shared_n = sum(MESuSiE_MLK_Signal==1),
            MESuSiEmlk_WB_n = sum(MESuSiE_MLK_Signal==2),
            MESuSiEmlk_BB_n = sum(MESuSiE_MLK_Signal==3),
            MESuSiEglm_Shared_n = sum(MESuSiE_GLM_Signal==1),
            MESuSiEglm_WB_n = sum(MESuSiE_GLM_Signal==2),
            MESuSiEglm_BB_n = sum(MESuSiE_GLM_Signal==3))

Signal_number<-Signal_number%>%
  pivot_longer(cols = -c(PHENONAME), names_to = "Cat", values_to = "Num")%>%
  separate(Cat, into = c("Method", "Signal"), sep = "_", extra = "merge") %>%
  mutate(
    Method = case_when(
      #str_detect(Method, "MESuSiEmlk10") ~ "MESuSiE_mlk_10",
      #str_detect(Method, "MESuSiEglm10") ~ "MESuSiE_glm_10",
      str_detect(Method, "MESuSiEmlk") ~ "MESuSiE_mlk",
      str_detect(Method, "MESuSiEglm") ~ "MESuSiE_glm",
      str_detect(Method, "MESuSiE") ~ "MESuSiE",
      str_detect(Method, "SuSiE") ~ "SuSiE",
      #str_detect(Method, "Paintor10") ~ "Paintor_annot",
      str_detect(Method, "Paintor") ~ "Paintor",
      TRUE ~ Method
    ),
    Signal = case_when(
      str_detect(Signal, "BB_n") ~ "AFR",
      str_detect(Signal, "WB_n") ~ "EUR",
      str_detect(Signal, "Shared_n") ~ "Shared",
      str_detect(Signal, "Either") ~ "Either",
      TRUE ~ Signal
    )
  )
Signal_number<-Signal_number%>%group_by(Method,PHENONAME)%>%mutate(prop = Num/sum(Num)*100,ypos = cumsum(prop)- 0.5*prop)%>%mutate(label = paste0(Signal," ",Num))
Signal_number<-Signal_number%>%mutate(Method = factor(Method, levels = c("MESuSiE","MESuSiE_mlk","MESuSiE_glm", "SuSiE", "Paintor")),Signal = factor(Signal, levels = c("EUR","AFR","Shared","Either")))

signal_num_plot<-ggplot(Signal_number, aes(x="", y=prop, fill=Signal)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  geom_text(aes(y = ypos, label = label),  size=7/14*5,color="white") +
  scale_fill_manual(values=c("#ABDB9F","#F2C1B6","#6162B0","gray"))+
  facet_grid(vars(PHENONAME),vars(Method),labeller=label_parsed)+theme(
    strip.text.x = element_text(size = 7,face="bold"),
    strip.text.y = element_text(size = 7,face="bold"),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank())
ggsave(paste0(plot_dir,"signal_prop_pie_bmi.pdf"),signal_num_plot,dpi=500,height=130,width  =250,unit='mm')




#######################################
#
#
#   source data processing and gathering
#
#
######################################



library(ggpubr)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggpmisc)
library(VennDiagram)
library(gridExtra)
library(ggbreak)
library(DescTools)
library(coin)
library(susieR)
library(ggrepel)
library(stringr)


load("/net/fantasia/home/borang/Susie_Mult/Revision_Round_1/01_06_Real_Data/summary_res/res.RData")
############################################################################################################
#
#   Stats to report
#
#############################################################################################################
all_sets_info<-data.frame(res_all%>%group_by(PHENONAME,Region) %>% summarise(across(c("MESuSiE_cs", "SuSiE_cs","Paintor_cs"), ~ sum(.x, na.rm = TRUE))))%>%filter(MESuSiE_cs!=0, SuSiE_cs!=0, Paintor_cs!=0) ###Median Set Size across all locus
all_sets_info_long<-all_sets_info%>%pivot_longer(!(PHENONAME|Region), names_to = "Method", values_to = "Count")
all_sets_info_long$Method<-factor(all_sets_info_long$Method,levels=c("MESuSiE_cs","SuSiE_cs","Paintor_cs"))
levels(all_sets_info_long$Method)<-c("MESuSiE","SuSiE","Paintor")

MESuSiE_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(MESuSiE_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB))))
SuSiE_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(SuSiE_cs==1)%>%summarise(zmax =median(pmax(abs(zscore_WB),abs(zscore_BB))))%>%pull(zmax)
Paintor_cs_Z<-res_all%>%group_by(PHENONAME) %>%filter(Paintor_cs==1)%>%summarise(zmax = median(pmax(abs(zscore_WB),abs(zscore_BB))))%>%pull(zmax)
set_size_z_info<-data.frame(cbind(MESuSiE_cs_Z,SuSiE_cs_Z,Paintor_cs_Z))
colnames(set_size_z_info)<-c("PHENONAME",c("MESuSiE","SuSiE","Paintor"))
set_size_z_info_long<-set_size_z_info %>%pivot_longer(!(PHENONAME), names_to = "Method", values_to = "Z")%>%mutate(Method = factor(Method, levels=c("MESuSiE","SuSiE","Paintor")))

ann_col_name<-c("missense", "synonymous", "utr_comb", "promotor", "CRE","brain_ind_eQTL")
# Functions for calculating fold enrichment
calc_fold_enrichment <- function(df, cs_col, ann_col_name) {
  df %>%
    group_by(Region) %>%
    filter(sum(!!sym(cs_col)) != 0) %>%
    group_by(PHENONAME, !!sym(cs_col)) %>%
    summarise(across(ann_col_name, ~ sum(.x, na.rm = TRUE) / n())) %>%
    group_by(PHENONAME) %>%
    summarise(across(ann_col_name, ~ .x[!!sym(cs_col) == 1] / .x[!!sym(cs_col) == 0]))
}

MESuSiE_PIP_ann <- calc_fold_enrichment(res_all, "MESuSiE_cs", ann_col_name)
SuSiE_PIP_ann <- calc_fold_enrichment(res_all, "SuSiE_cs", ann_col_name)
Paintor_PIP_ann <- calc_fold_enrichment(res_all, "Paintor_cs", ann_col_name)
# Combine results
PHENONAME_CS_enrichment <- bind_rows(
  MESuSiE_PIP_ann %>% mutate(Method = "MESuSiE"),
  SuSiE_PIP_ann %>% mutate(Method = "SuSiE"),
  Paintor_PIP_ann %>% mutate(Method = "Paintor")
) %>% mutate(Method = factor(Method, levels = c("MESuSiE", "SuSiE", "Paintor")))%>%
  dplyr::select(PHENONAME,brain_ind_eQTL ,Method )%>%dplyr::rename(eQTL = brain_ind_eQTL)
# Pivot to long format
PHENONAME_CS_enrichment_long <- PHENONAME_CS_enrichment %>%
  pivot_longer(cols = -c(Method, PHENONAME), names_to = "Cat", values_to = "Prop") %>%
  mutate(Method = factor(Method, levels = c("MESuSiE", "SuSiE", "Paintor")))

PHENONAME_CS_enrichment_long%>%group_by(Method)%>%summarise(min(Prop),max(Prop),mean(Prop))


library(data.table)
library(ggpubr)
library(dplyr)
source_data_figure_dir<-paste0("/net/fantasia/home/borang/Susie_Mult/Revision_Round_1/source_data_figure/")
system(paste0("mkdir -p ",source_data_figure_dir)) 
library(xlsx)
xlsx::write.xlsx(data.frame(all_sets_info_long),paste0(source_data_figure_dir,"Fig5.xlsx"),
                 sheetName="set_size",
                 col.names = TRUE,row.names=F,showNA=TRUE,append=TRUE)
xlsx::write.xlsx(data.frame(set_size_z_info),paste0(source_data_figure_dir,"Fig5.xlsx"),
                 sheetName="Zscore",
                 col.names = TRUE,row.names=F,showNA=TRUE,append=TRUE)

xlsx::write.xlsx(data.frame(PHENONAME_CS_enrichment_long),paste0(source_data_figure_dir,"Fig5.xlsx"),
                 sheetName="Fold_Enrichment",
                 col.names = TRUE,row.names=F,showNA=TRUE,append=TRUE)


