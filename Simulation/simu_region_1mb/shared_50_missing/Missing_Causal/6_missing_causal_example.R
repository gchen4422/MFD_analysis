library(dplyr)
library(data.table)
library(MESuSiE)
library(susieR)
library(ggplot2)
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
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Causal/MFD_utility.R")

#args<-as.numeric(commandArgs(TRUE))
LD_BLOCK = 46
h2_num = 2
num_causal = 1
causal_index = 2
External_index = 1
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
EU_cov<-as.matrix(fread(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/causal_num_",num_causal,"/summary_data/CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")))
colnames(EU_cov)<-zfile$RSID
rownames(EU_cov)<-zfile$RSID
BB_cov<-as.matrix(fread(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_all/causal_num_",num_causal,"/summary_data/CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")))
colnames(BB_cov)<-zfile$RSID
rownames(BB_cov)<-zfile$RSID


if(External_index ==1 ){
  num_snp_missing = round(nrow(zfile)*0.3)
}else if(External_index ==2){
  num_snp_missing = round(nrow(zfile)*0.3+nrow(zfile)*(1-0.3)*0.3)
}


non_causal_SNP = zfile%>%filter(Signal==0)%>%pull(RSID)
EU_missing_SNP = c(sample(non_causal_SNP,num_snp_missing-1),zfile%>%filter(Signal!=0)%>%pull(RSID))
BB_missing_SNP = EU_missing_SNP
zfile<-zfile%>%mutate(EU_missing = ifelse(RSID%in%EU_missing_SNP,1,0),BB_missing = ifelse(RSID%in%BB_missing_SNP,1,0))

# DETERMINE BLOCK ID:
current_block_id <- LD_BLOCK

# Use tryCatch to skip regions where files might be missing without stopping the whole script
skip_to_next <- FALSE

cat(paste0("Processing Block: ", current_block_id, "... "))


# --- Run Decision Logic ---

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


summary_stat_1 = data.frame("SNP" = susie_EU_zfile$RSID,"CHR" = susie_EU_zfile$CHR,"POS" = susie_EU_zfile$POS,"Signal" = susie_EU_zfile$Signal,"Beta"=susie_EU_zfile$zscore_1/sqrt(num_EU),"Se"=1/sqrt(num_EU), "Z" =susie_EU_zfile$zscore_1 ,  "N" =num_EU ) %>% mutate(PVAL = 2*pnorm(abs(Z),lower.tail = F))
summary_stat_2 = data.frame("SNP" = susie_BB_zfile$RSID,"CHR" = susie_BB_zfile$CHR,"POS" = susie_BB_zfile$POS,"Signal" = susie_BB_zfile$Signal,"Beta"=susie_BB_zfile$zscore_2/sqrt(num_BB),"Se"=1/sqrt(num_BB), "Z" =susie_BB_zfile$zscore_2 ,  "N" =num_BB ) %>% mutate(PVAL = 2*pnorm(abs(Z),lower.tail = F))

setDT(summary_stat_1)
setDT(summary_stat_2)

decision <- decide_finemapping_method(summary_stat_1, summary_stat_2, susie_EU_cov, susie_BB_cov)

mfd_result = run_mf_decision_fm_test_merge_cs(summary_stat_1,summary_stat_2,susie_EU_cov,susie_BB_cov, pop_names = c("EUR", "AFR"))
save(summary_stat_1, summary_stat_2, susie_EU_cov, susie_BB_cov, file = "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/website/source_data/example.RData")


lead_SNP = zfile[zfile$Signal != 0,]$RSID
# --- EUR Handling ---
lead_SNP_index <- which(colnames(susie_EU_cov) == lead_SNP)

if (length(lead_SNP_index) == 0) {
  # If missing: create a vector of zeros matching the matrix dimensions/names
  eur_ld_vec <- setNames(rep(0, ncol(susie_EU_cov)), colnames(susie_EU_cov))
} else {
  # If found: extract column and square it
  eur_ld_vec <- setNames((susie_EU_cov[, lead_SNP_index])^2, colnames(susie_EU_cov))
}

# --- AFR Handling ---
lead_SNP_index <- which(colnames(susie_BB_cov) == lead_SNP)

if (length(lead_SNP_index) == 0) {
  # If missing: create a vector of zeros matching the matrix dimensions/names
  afr_ld_vec <- setNames(rep(0, ncol(susie_BB_cov)), colnames(susie_BB_cov))
} else {
  # If found: extract column and square it
  afr_ld_vec <- setNames((susie_BB_cov[, lead_SNP_index])^2, colnames(susie_BB_cov))
}

summary_stat_1 <- summary_stat_1 %>%
  mutate(
    # Retrieve the r2 value using the SNP ID as the key
    # Note: Replace 'SNP' with 'rsid' or 'variant' if that is the column name
    r2_EUR = eur_ld_vec[as.character(SNP)]
  ) %>%
  mutate(
    # If the SNP wasn't in the covariance matrix (result is NA), set to 0
    r2_EUR = ifelse(is.na(r2_EUR), 0, r2_EUR)
  )

summary_stat_2 <- summary_stat_2 %>%
  mutate(
    # Retrieve the r2 value using the SNP ID as the key
    # Note: Replace 'SNP' with 'rsid' or 'variant' if that is the column name
    r2_AFR = afr_ld_vec[as.character(SNP)]
  ) %>%
  mutate(
    # If the SNP wasn't in the covariance matrix (result is NA), set to 0
    r2_AFR = ifelse(is.na(r2_AFR), 0, r2_AFR)
  )

custom_theme <- function() {
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 7, face="bold"),
    axis.title.y = element_text(size = 7, face="bold"),
    strip.text.x     = element_text(size = 7, face="bold"),   # ← bump this
    strip.text.y     = element_text(size = 7, face="bold"),   # ← and/or this
    strip.background = element_blank(),
    legend.text = element_text(size=7),
    legend.title = element_text(size=7, face="bold"),
    plot.title = element_text(size=10, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.tag = element_text(size = 7),
    axis.line = element_line(color = "black")
  )
}

gwas_plot_fun <- function(data_plot, xlab_name, ylab_name, yintercept) {
  
  p_manhattan = ggplot() + geom_point(data = data_plot%>%filter(Lead_SNP==0), aes(x = POS, y = PIP, color = r2), size = 1)
  p_manhattan = p_manhattan + geom_point(data = data_plot%>%filter(Lead_SNP==1), aes(x = POS, y = PIP), size = 1.5, color = "red") +
    geom_text_repel(data = data_plot%>%filter(Lead_SNP==1), mapping = aes(x = POS, y = PIP, label = SNP), vjust = 1.2, size = 7/10*3, show.legend =FALSE) 
  p_manhattan = p_manhattan +
    scale_color_stepsn(
      colors = c("navy", "lightskyblue", "green", "orange", "red"),
      breaks = seq(0.2, 0.8, by = 0.2),
      limits = c(0, 1),
      show.limits = TRUE,
      na.value = 'grey50',
      name = expression(R^2)
    )
  p_manhattan = p_manhattan +
    geom_hline(
      yintercept = yintercept,
      linetype = "dashed",
      color = "grey50",
      size = 0.5
    ) 
  p_manhattan = p_manhattan +
    geom_vline(
      xintercept = data_plot%>%filter(lead_SNP==1)%>%pull(POS),
      linetype = "dashed",
      color = "grey50",
      size = 0.5
    ) 
  p_manhattan = p_manhattan + xlim(min(data_plot$POS),max(data_plot$POS))
  p_manhattan = p_manhattan + expand_limits(x = round(max(data_plot$POS)/1e6)*1e6)
  if(max(data_plot$POS>1e6)){
    p_manhattan = p_manhattan + scale_x_continuous(labels = function(x) paste0(x / 1e6, " MB"))
  }
  if(max(data_plot$POS<1e6)){
    p_manhattan = p_manhattan + scale_x_continuous(labels = function(x) paste0(x / 1e3, " KB"))
  }
  p_manhattan = p_manhattan + xlab(xlab_name) +ylab(ylab_name)
  p_manhattan = p_manhattan + guides(fill = guide_legend(title = as.expression(bquote(R^2))))
  p_manhattan = p_manhattan + theme_bw()+custom_theme()
  return(p_manhattan)
}


EUR_GWAS_plot_data<-summary_stat_1%>%mutate(r2 = r2_EUR,PIP = -log10(2*pnorm(-abs(Z))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
AFR_GWAS_plot_data<-summary_stat_2%>%mutate(r2 = r2_AFR,PIP = -log10(2*pnorm(-abs(Z))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
p_EUR<-gwas_plot_fun (EUR_GWAS_plot_data, "All of Us EA", "-log10(P-value)", -log10(5e-8))+ theme(legend.position = "none")
p_AFR<-gwas_plot_fun (AFR_GWAS_plot_data, "All of Us BA", "-log10(P-value)", -log10(5e-8))+ theme(legend.position = "none")



wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Causal/",External_index_name,causal_index_name,"/causal_num_",num_causal,"/")	
data_dir<-paste0(wrk_dir,"summary_data/")
data_dir_ld = paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_missing/Missing_Non_Causal/",External_index_name,causal_index_name,"/causal_num_",num_causal,"/summary_data/")
result_dir<-paste0(wrk_dir,"result/")    
setwd(result_dir)



SuSiE_res_name<-paste0(result_dir,"MESuSiE_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".RData")
load(SuSiE_res_name)
mf_decision = fread(paste0(result_dir, "Finemapping_Method_Decisions_h2_",h2_num, ".csv"))
is_mesusie <- mf_decision$Method[LD_BLOCK] == "MESuSiE"


zfile_MESuSiE_SuSiE_Paintor = zfile_MESuSiE_SuSiE_Paintor %>%
  mutate(
    # 1. Define SuSiE_PIP_WB first so it is available for use below
    SuSiE_PIP_WB = SuSiE_PIP_EU,
    
    # 2. Conditionally assign the MFD columns line-by-line
    # Use standard 'if' because 'is_mesusie' is a single value, not a vector.
    MFD_PIP_Either = if (is_mesusie) MESuSiE_PIP_Either else SuSiE_PIP_Either,
    MFD_PIP_WB     = if (is_mesusie) MESuSiE_PIP_WB     else SuSiE_PIP_WB,
    MFD_PIP_BB     = if (is_mesusie) MESuSiE_PIP_BB     else SuSiE_PIP_BB,
    MFD_PIP_Shared = if (is_mesusie) MESuSiE_PIP_Shared else SuSiE_PIP_Shared
    
  ) 
zfile_MESuSiE_SuSiE_Paintor = zfile_MESuSiE_SuSiE_Paintor %>% mutate(r2_AFR = summary_stat_2$r2_AFR, r2_EUR = 0)
zfile_MESuSiE_SuSiE_Paintor<-zfile_MESuSiE_SuSiE_Paintor%>%mutate(SuSiE_cat = case_when(SuSiE_PIP_EU>0.5&SuSiE_PIP_BB>0.5 ~ 3,
                                                                  SuSiE_PIP_EU>0.5&SuSiE_PIP_BB<0.5 ~ 1,
                                                                  SuSiE_PIP_EU<0.5&SuSiE_PIP_BB>0.5 ~ 2,
                                                                  TRUE ~ 0),
                                            MFD_cat = case_when(MFD_PIP_WB>0.5~1,
                                                                MFD_PIP_BB>0.5~2,
                                                                MFD_PIP_Shared>0.5~3,
                                                                TRUE~0),
                                            MESuSiE_cat = case_when(MESuSiE_PIP_WB>0.5~1,
                                                                    MESuSiE_PIP_BB>0.5~2,
                                                                    MESuSiE_PIP_Shared>0.5~3,
                                                                    TRUE~0),
                                            
)
MESuSiE_plot_data<-zfile_MESuSiE_SuSiE_Paintor%>%mutate(r2 = r2_AFR,PIP = MESuSiE_PIP_Either,Lead_SNP = ifelse(RSID==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MESuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(RSID,POS, r2,PIP,Lead_SNP,cat) %>% dplyr::rename(SNP=RSID)
MFD_plot_data<-zfile_MESuSiE_SuSiE_Paintor%>%mutate(r2 = r2_AFR,PIP = MFD_PIP_Either,Lead_SNP = ifelse(RSID==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MFD_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(RSID,POS, r2,PIP,Lead_SNP,cat) %>% dplyr::rename(SNP=RSID)


###Function used for PIP plot	
finemap_plot_fun<-function(data_plot,xlab_name,ylab_name,yintercept){
  pos_rng <- range(data_plot$POS, na.rm = TRUE)
  
  p_manhattan = ggplot() + geom_point(data = data_plot, aes(x = POS, y = PIP, color = r2,shape = cat))+scale_shape_manual(name="Category",drop=FALSE,values=c(20,24,25,23,22))
  p_manhattan = p_manhattan + geom_text_repel(data =data_plot%>%filter(Lead_SNP==1), mapping=aes(x=POS, y=PIP, label=SNP),vjust=1.2, size= 7/10*3,show.legend = FALSE)
  p_manhattan = p_manhattan + theme_bw()+scale_color_stepsn(
    colors = c("navy", "lightskyblue", "green", "orange", "red"),
    breaks = seq(0.2, 0.8, by = 0.2),
    limits = c(0, 1),
    show.limits = TRUE,
    na.value = 'grey50',
    name = expression(R^2)
  )
  p_manhattan = p_manhattan + geom_hline(
    yintercept =yintercept,
    linetype = "dashed",
    color = "grey50",
    size = 0.5
  ) + geom_vline(
    xintercept = data_plot%>%filter(lead_SNP==1)%>%pull(POS),
    linetype = "dashed",
    color = "grey50",
    size = 0.5
  ) 
  p_manhattan = p_manhattan + xlim(min(data_plot$POS),max(data_plot$POS))
  p_manhattan = p_manhattan + expand_limits(x = round(max(data_plot$POS)/1e6)*1e6)
  if(max(data_plot$POS>1e6)){
    p_manhattan = p_manhattan + scale_x_continuous(labels = function(x) paste0(x / 1e6, " MB"),expand = expansion(mult = c(0.05, 0.05)), limits = pos_rng)
  }
  if(max(data_plot$POS<1e6)){
    p_manhattan = p_manhattan + scale_x_continuous(labels = function(x) paste0(x / 1e3, " KB"),expand = expansion(mult = c(0.05, 0.05)), limits = pos_rng)
  }
  p_manhattan= p_manhattan+xlab(xlab_name)+ylab(ylab_name)
  p_manhattan= p_manhattan+guides(fill=guide_legend(title=as.expression(bquote(R^2))))
  p_manhattan = p_manhattan + theme_bw()+custom_theme()
  return(p_manhattan)
}	

p_MESuSiE<-finemap_plot_fun(MESuSiE_plot_data, "MESuSiE", "PIP", 0.5) + theme(legend.position = "none")
p_MFD<-finemap_plot_fun(MFD_plot_data, "MFD", "PIP", 0.5) + theme(legend.position = "none")
