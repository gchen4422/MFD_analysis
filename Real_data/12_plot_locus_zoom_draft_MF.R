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
plot_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/Figure/"
Gene_List<-fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Annotation/Gencode_genelist/Gencode_GRCh38_Genes_UniqueList2025.txt",header=T)

custom_theme <- function() {
  theme(
    axis.text.x = element_text(size = 5),
    axis.text.y = element_text(size = 5),
    axis.title.x = element_text(size = 8, face="bold"),
    axis.title.y = element_text(size = 8, face="bold"),
    strip.text.x     = element_text(size = 10, face="bold"),   # ← bump this
    strip.text.y     = element_text(size = 10, face="bold"),   # ← and/or this
    strip.background = element_blank(),
    legend.text = element_text(size=12),
    legend.title = element_text(size=12, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.tag = element_text(size = 12),
    axis.line = element_line(color = "black")
  )
}


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
##############################################################
#
#
#       Real Data Example Plotter
#
#
################################################################

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
# Function used for gene plot 
gene_range_plot_fun<-function(gene_list_data,plot.range){
  p<-ggplot(data = gene_list_data) +
    geom_linerange(aes(x = Gene, ymin = Start, ymax = End))+ylim(plot.range)+ expand_limits(y = round(max(plot.range[2])/1e6)*1e6)+scale_y_continuous(limits = plot.range,labels = function(y) paste0(y / 1e6, " MB"),expand = expansion(mult = c(0.05, 0.05)),)+coord_flip()+
    geom_text(aes(x = Gene, y = Start, label = Gene), hjust = "right", size = 5*3/10) + ylab(paste0("chr",unique(gsub("chr","",gene_list_data$Chrom))))+ xlab("Gene") + 
    theme_bw() +  theme(
      axis.text.x = element_text(size = 7),
      axis.text.y = element_blank(),  
      axis.ticks.y =  element_blank(), 
      axis.title.x = element_text(size = 7, face="bold"),
      axis.title.y = element_text(size = 7, face="bold"),
      strip.text.x     = element_text(size = 7, face="bold"),   # ← bump this
      strip.text.y     = element_text(size = 7, face="bold"),   # ← and/or this
      strip.background = element_blank(),
      legend.text = element_text(size=12),
      legend.title = element_text(size=12, face="bold"),
      plot.title = element_text(size=14, hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(), 
      axis.line.x  = element_line(color = "black"),
      axis.line.y  = element_line(color = "black")
    )
  return(p)
}







###################################################################
#
#
#   LOCI 11 (rs9821489) for DBP is an example of how MFD(susie-post hoc) can find AS-V but mesusie cannot
#   it is filtered out by taking common subset in mesusie        
#
##################################################################

drop_guides <- guides(
  color = "none",
  shape = "none",
  fill  = "none",
  linetype = "none"
)


LD_BLOCK = 11
trait = "DBP"
load(paste0(res_dir,"res_all_mf.RData"))   
res_all = res_all %>% filter(PHENONAME == "DBP")


wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/results_",trait,"/")
data_dir<-paste0(wrk_dir,"summary_data/")
out_dir<-paste0(wrk_dir,"out/")



#res_all =fread(paste0(out_dir,trait,"_REF_ALL"))

EU_cov<-as.matrix(fread(paste0(data_dir,"eur_loci_",LD_BLOCK,".ld")))
BB_cov<-as.matrix(fread(paste0(data_dir,"afr_loci_",LD_BLOCK,".ld")))
#rownames(EU_cov) = rownames(BB_cov) = colnames(EU_cov)

candidate_region<-res_all%>%filter(Region==LD_BLOCK)
#candidate_region = candidate_region %>% filter(POS > 71995699)



#EU_cov = as.matrix(EU_cov[candidate_region$SNP,candidate_region$SNP])
#BB_cov = as.matrix(BB_cov[candidate_region$SNP,candidate_region$SNP])

# rs3184504 is the eQTL of SH2B3 gene, highlighted in the paper
lead_SNP = "rs9821489"
# --- EUR Handling ---
lead_SNP_index <- which(colnames(EU_cov) == lead_SNP)

if (length(lead_SNP_index) == 0) {
  # If missing: create a vector of zeros matching the matrix dimensions/names
  eur_ld_vec <- setNames(rep(0, ncol(EU_cov)), colnames(EU_cov))
} else {
  # If found: extract column and square it
  eur_ld_vec <- setNames((EU_cov[, lead_SNP_index])^2, colnames(EU_cov))
}

# --- AFR Handling ---
lead_SNP_index <- which(colnames(BB_cov) == lead_SNP)

if (length(lead_SNP_index) == 0) {
  # If missing: create a vector of zeros matching the matrix dimensions/names
  afr_ld_vec <- setNames(rep(0, ncol(BB_cov)), colnames(BB_cov))
} else {
  # If found: extract column and square it
  afr_ld_vec <- setNames((BB_cov[, lead_SNP_index])^2, colnames(BB_cov))
}
candidate_region <- candidate_region %>%
  mutate(
    # Retrieve the r2 value using the SNP ID as the key
    # Note: Replace 'SNP' with 'rsid' or 'variant' if that is the column name
    r2_EUR = eur_ld_vec[as.character(SNP)],
    r2_AFR = afr_ld_vec[as.character(SNP)]
  ) %>%
  mutate(
    # If the SNP wasn't in the covariance matrix (result is NA), set to 0
    r2_EUR = ifelse(is.na(r2_EUR), 0, r2_EUR),
    r2_AFR = ifelse(is.na(r2_AFR), 0, r2_AFR),
    POS = as.numeric(POS)
  )

####Category Setting
candidate_region<-candidate_region%>%mutate(SuSiE_cat = case_when(SuSiE_PIP_EU>0.5&SuSiE_PIP_BB>0.5 ~ 3,
                                                                  SuSiE_PIP_EU>0.5&SuSiE_PIP_BB<0.5 ~ 1,
                                                                  SuSiE_PIP_EU<0.5&SuSiE_PIP_BB>0.5 ~ 2,
                                                                  TRUE ~ 0),
                                            MFD_cat = case_when(PIP_Ancestry_1>0.5~1,
                                                                PIP_Ancestry_2>0.5~2,
                                                                PIP_Shared>0.5~3,
                                                                TRUE~0),
                                            Paintor_cat = case_when(Paintor_PIP_Either>0.5~4,
                                                                    TRUE~0),
                                            MESuSiE_cat = case_when(MESuSiE_PIP_WB>0.5~1,
                                                                    MESuSiE_PIP_BB>0.5~2,
                                                                    MESuSiE_PIP_Shared>0.5~3,
                                                                    TRUE~0),
                                            SuSiEx_cat = case_when(SuSiEx_PIP_Either>0.5~4,
                                                                   TRUE~0),
                                            XMAP_cat = case_when(XMAP_PIP_Either>0.5~4,
                                                                 TRUE~0),
                                            CARMAx_cat = case_when(CARMAX_PIP_Shared>0.5~1,
                                                                   CARMAX_PIP_WB>0.5~2,
                                                                   CARMAX_PIP_BB>0.5~3,
                                                                   TRUE~0),
                                            SuSiE_merged_cat= case_when(SuSiE_merged_PIP_Either>0.5~4,
                                                                        TRUE~0)
                                            
)

###Stat of the region
candidate_region <- candidate_region %>%
  mutate(
    Z_eur = coalesce(Z_eur, 0),
    Z_afr = coalesce(Z_afr, 0)
  )
candidate_region%>%summarise(min(as.numeric(POS))/1024/1024,max(as.numeric(POS))/1024/1024)
candidate_region%>%filter(SNP=="rs9821489")%>%summarise(2*pnorm(-abs(Z_eur)),2*pnorm(-abs(Z_afr)))
candidate_region%>%mutate(Pvalue_BB =2*pnorm(-abs(Z_afr)) )%>%filter(Pvalue_BB<5e-8)
candidate_region %>%
  filter(SNP == "rs9821489") %>%
  dplyr::select(
    # SuSiE
    SuSiE_PIP_EU, SuSiE_PIP_BB,SuSiE_PIP_Either,SuSiE_PIP_Shared,
    # MFD
    PIP_Ancestry_1, PIP_Ancestry_2, PIP_Shared,PIP_Either,
    # Paintor
    Paintor_PIP_Either,
    # MESuSiE
    MESuSiE_PIP_WB, MESuSiE_PIP_BB, MESuSiE_PIP_Shared,
    # SuSiEx
    SuSiEx_PIP_Either,
    # XMAP
    XMAP_PIP_Either,
    # CARMAx
    CARMAX_PIP_Shared, CARMAX_PIP_WB, CARMAX_PIP_BB,
    # SuSiE Merged
    SuSiE_merged_PIP_Either
  )

#candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(n())
#candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(Z_eur))
#candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(2*pnorm(-abs(Z_eur))))

candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise(summary(Z_afr))
candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise((2*pnorm(-abs(Z_afr))))
candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise(r2_AFR)

candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%pull(Z_afr)


#annotation_results <- create_annotation_comparison(candidate_region)


###GWAS PLOT

candadate_metal = fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/results_DBP/result/Metal_DBP_result1.tbl") %>% 
  filter(MarkerName %in% candidate_region$CHR_POS)%>% 
  arrange(match(MarkerName, candidate_region$CHR_POS))

meta_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = -log10(candadate_metal$`P-value`),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
EUR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = -log10(2*pnorm(-abs(Z_eur))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
AFR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = -log10(2*pnorm(-abs(Z_afr))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
p_EUR<-gwas_plot_fun (EUR_GWAS_plot_data, "All of Us EA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")
p_AFR<-gwas_plot_fun (AFR_GWAS_plot_data, "All of Us BA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")
p_meta<-gwas_plot_fun (meta_GWAS_plot_data, "All of Us Meta", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")


###Finemap Plot
#EUR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = susie_EU,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
#AFR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = susie_BB,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
Paintor_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
SuSiE_merged_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = SuSiE_merged_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_merged_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
MESuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MESuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
MFD_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MFD_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)


EUR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = SuSiE_PIP_EU,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
AFR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = SuSiE_PIP_BB,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
SuSiEx_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = SuSiEx_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
XMAP_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = XMAP_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_merged_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
CARMAX_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = CARMAX_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MFD_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)



p_Paintor<-finemap_plot_fun(Paintor_plot_data, "Paintor", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_MESuSiE<-finemap_plot_fun(MESuSiE_plot_data, "MESuSiE", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_SuSiE_merged<-finemap_plot_fun(SuSiE_merged_plot_data, "SuSiE_Meta", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_MFD<-finemap_plot_fun(MFD_plot_data, "MFD", "PIP", 0.5)#+ drop_guides + theme(legend.position = "none")


p_EUR_SuSiE<-finemap_plot_fun(EUR_SuSiE_plot_data, "SuSiE EUR", "PIP", 0.5)
p_AFR_SuSiE<-finemap_plot_fun(AFR_SuSiE_plot_data, "SuSiE AFR", "PIP", 0.5)
p_SuSiEx<-finemap_plot_fun(SuSiEx_plot_data, "SuSiEx", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_XMAP<-finemap_plot_fun(XMAP_plot_data, "XMAP", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_CARMAX<-finemap_plot_fun(CARMAX_plot_data, "CARMA-X", "PIP", 0.5)#+ drop_guides + theme(legend.position = "none")


# Gene Plot
plot.range <- c(min(candidate_region$POS), max(candidate_region$POS))
Gene_List_sub_coding<-Gene_List%>%filter(Chrom==paste0("chr",unique(candidate_region$CHR)))%>%filter(Start<max(candidate_region$POS),End>min(candidate_region$POS))%>%filter(Coding=="proteincoding")%>%filter(!is.na(cdsLength))%>%filter(GeneLength>15000)
p2<-gene_range_plot_fun(Gene_List_sub_coding,plot.range)

##Combine Plot together
combined_plot<-(p_EUR/p_MESuSiE/p_Paintor/p2+plot_layout(heights = c(1,1,1,1))|p_AFR/p_MFD/p_SuSiE_merged/p2+plot_layout(heights = c(1,1,1,1)))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")
combined_plot
p_ASV = combined_plot
#ggsave(paste0(plot_dir,"Region_",LD_BLOCK ,"_",trait,".pdf"),combined_plot,width = 200, height = 180,units="mm",dpi=500)

combined_plot_supp = (p_meta/p_EUR_SuSiE/p_XMAP+plot_layout(heights = c(1,1,1))|p_SuSiEx/p_AFR_SuSiE/p_CARMAX+plot_layout(heights = c(1,1,1)))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")

ggsave(paste0(plot_dir,"Supp_Region_",LD_BLOCK ,"_",trait,".pdf"),combined_plot_supp,width = 210, height = 180,units="mm",dpi=500)

###################################################################
#
#   LOCI 31 (rs4932179) for 
#   DBP is an example of how joint modelings improve power of fine-mapping compared to meta-analysis
#  
#   
#
##################################################################

drop_guides <- guides(
  color = "none",
  shape = "none",
  fill  = "none",
  linetype = "none"
)


LD_BLOCK = 31
trait = "DBP"
load(paste0(res_dir,"res_all_mf.RData"))   
res_all = res_all %>% filter(PHENONAME == trait)


wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/results_",trait,"/")
data_dir<-paste0(wrk_dir,"summary_data/")
out_dir<-paste0(wrk_dir,"out/")


#res_all =fread(paste0(out_dir,trait,"_REF_ALL"))
candidate_region<-res_all%>%filter(Region==LD_BLOCK)

EU_cov<-as.matrix(fread(paste0(data_dir,"eur_loci_",LD_BLOCK,".ld")))
BB_cov<-as.matrix(fread(paste0(data_dir,"afr_loci_",LD_BLOCK,".ld")))
#rownames(EU_cov) = rownames(BB_cov) = colnames(EU_cov)

#candidate_region = candidate_region %>% filter(POS > 71995699)



#EU_cov = as.matrix(EU_cov[candidate_region$SNP,candidate_region$SNP])
#BB_cov = as.matrix(BB_cov[candidate_region$SNP,candidate_region$SNP])

# rs3184504 is the eQTL of SH2B3 gene, highlighted in the paper
lead_SNP = "rs4932179"
# --- EUR Handling ---
lead_SNP_index <- which(colnames(EU_cov) == lead_SNP)

if (length(lead_SNP_index) == 0) {
  # If missing: create a vector of zeros matching the matrix dimensions/names
  eur_ld_vec <- setNames(rep(0, ncol(EU_cov)), colnames(EU_cov))
} else {
  # If found: extract column and square it
  eur_ld_vec <- setNames((EU_cov[, lead_SNP_index])^2, colnames(EU_cov))
}

# --- AFR Handling ---
lead_SNP_index <- which(colnames(BB_cov) == lead_SNP)

if (length(lead_SNP_index) == 0) {
  # If missing: create a vector of zeros matching the matrix dimensions/names
  afr_ld_vec <- setNames(rep(0, ncol(BB_cov)), colnames(BB_cov))
} else {
  # If found: extract column and square it
  afr_ld_vec <- setNames((BB_cov[, lead_SNP_index])^2, colnames(BB_cov))
}

candidate_region <- candidate_region %>%
  mutate(
    # Retrieve the r2 value using the SNP ID as the key
    # Note: Replace 'SNP' with 'rsid' or 'variant' if that is the column name
    r2_EUR = eur_ld_vec[as.character(SNP)],
    r2_AFR = afr_ld_vec[as.character(SNP)]
  ) %>%
  mutate(
    # If the SNP wasn't in the covariance matrix (result is NA), set to 0
    r2_EUR = ifelse(is.na(r2_EUR), 0, r2_EUR),
    r2_AFR = ifelse(is.na(r2_AFR), 0, r2_AFR),
    POS = as.numeric(POS)
  )

####Category Setting
candidate_region<-candidate_region%>%mutate(SuSiE_cat = case_when(SuSiE_PIP_EU>0.5&SuSiE_PIP_BB>0.5 ~ 3,
                                                                  SuSiE_PIP_EU>0.5&SuSiE_PIP_BB<0.5 ~ 1,
                                                                  SuSiE_PIP_EU<0.5&SuSiE_PIP_BB>0.5 ~ 2,
                                                                  TRUE ~ 0),
                                            MFD_cat = case_when(PIP_Ancestry_1>0.5~1,
                                                                PIP_Ancestry_2>0.5~2,
                                                                PIP_Shared>0.5~3,
                                                                TRUE~0),
                                            Paintor_cat = case_when(Paintor_PIP_Either>0.5~4,
                                                                    TRUE~0),
                                            MESuSiE_cat = case_when(MESuSiE_PIP_WB>0.5~1,
                                                                    MESuSiE_PIP_BB>0.5~2,
                                                                    MESuSiE_PIP_Shared>0.5~3,
                                                                    TRUE~0),
                                            SuSiEx_cat = case_when(SuSiEx_PIP_Either>0.5~4,
                                                                   TRUE~0),
                                            XMAP_cat = case_when(XMAP_PIP_Either>0.5~4,
                                                                 TRUE~0),
                                            CARMAx_cat = case_when(CARMAX_PIP_Shared>0.5~1,
                                                                   CARMAX_PIP_WB>0.5~2,
                                                                   CARMAX_PIP_BB>0.5~3,
                                                                   TRUE~0),
                                            SuSiE_merged_cat= case_when(SuSiE_merged_PIP_Either>0.5~4,
                                                                        TRUE~0)
                                            
)

###Stat of the region
candidate_region <- candidate_region %>%
  mutate(
    Z_eur = coalesce(Z_eur, 0),
    Z_afr = coalesce(Z_afr, 0)
  )
candidate_region%>%summarise(min(as.numeric(POS))/1024/1024,max(as.numeric(POS))/1024/1024)
candidate_region%>%filter(SNP=="rs4932179")%>%summarise(2*pnorm(-abs(Z_eur)),2*pnorm(-abs(Z_afr)))
candidate_region%>%mutate(Pvalue_BB =2*pnorm(-abs(Z_afr)) )%>%filter(Pvalue_BB<5e-8)
candidate_region %>%
  filter(SNP == "rs4932179") %>%
  dplyr::select(
    # SuSiE
    SuSiE_PIP_EU, SuSiE_PIP_BB,SuSiE_PIP_Either,SuSiE_PIP_Shared,
    # MFD
    PIP_Ancestry_1, PIP_Ancestry_2, PIP_Shared,PIP_Either,
    # Paintor
    Paintor_PIP_Either,
    # MESuSiE
    MESuSiE_PIP_WB, MESuSiE_PIP_BB, MESuSiE_PIP_Shared,
    # SuSiEx
    SuSiEx_PIP_Either,
    # XMAP
    XMAP_PIP_Either,
    # CARMAx
    CARMAX_PIP_Shared, CARMAX_PIP_WB, CARMAX_PIP_BB,
    # SuSiE Merged
    SuSiE_merged_PIP_Either
  )

candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(n())
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(Z_eur))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(2*pnorm(-abs(Z_eur))))
candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(r2_AFR))

candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise(n())
candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise(summary(Z_afr))
candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise((2*pnorm(-abs(Z_afr))))
candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise(r2_AFR)

candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%pull(Z_afr)


#annotation_results <- create_annotation_comparison(candidate_region)


###GWAS PLOT
candadate_metal = fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/results_DBP/result/Metal_DBP_result1.tbl") %>% 
  filter(MarkerName %in% candidate_region$CHR_POS)%>% 
  arrange(match(MarkerName, candidate_region$CHR_POS))

meta_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = -log10(candadate_metal$`P-value`),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
EUR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = -log10(2*pnorm(-abs(Z_eur))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
AFR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = -log10(2*pnorm(-abs(Z_afr))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
p_EUR<-gwas_plot_fun (EUR_GWAS_plot_data, "All of Us EA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")
p_AFR<-gwas_plot_fun (AFR_GWAS_plot_data, "All of Us BA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")
p_meta<-gwas_plot_fun (meta_GWAS_plot_data, "All of Us Meta", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")

###Finemap Plot
Paintor_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
SuSiE_merged_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = SuSiE_merged_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_merged_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
MESuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MESuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
MESuSiE_plot_data[which(MESuSiE_plot_data$SNP == "rs10744783"),]$cat = "Paintor"
MFD_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MFD_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)


EUR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = SuSiE_PIP_EU,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
AFR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = SuSiE_PIP_BB,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
SuSiEx_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = SuSiEx_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
XMAP_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = XMAP_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_merged_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
CARMAX_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = CARMAX_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MFD_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)

#p_EUR_SuSiE<-finemap_plot_fun(EUR_SuSiE_plot_data, "SuSiE EUR", "PIP", 0.5)
#p_AFR_SuSiE<-finemap_plot_fun(AFR_SuSiE_plot_data, "SuSiE AFR", "PIP", 0.5)
p_Paintor<-finemap_plot_fun(Paintor_plot_data, "Paintor", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_SuSiE_merged<-finemap_plot_fun(SuSiE_merged_plot_data, "SuSiE_Meta", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")

p_MESuSiE<-finemap_plot_fun(MESuSiE_plot_data, "MESuSiE", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_MFD<-finemap_plot_fun(MFD_plot_data, "MFD", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")



p_EUR_SuSiE<-finemap_plot_fun(EUR_SuSiE_plot_data, "SuSiE EUR", "PIP", 0.5)+ theme(legend.position = "none")
p_AFR_SuSiE<-finemap_plot_fun(AFR_SuSiE_plot_data, "SuSiE AFR", "PIP", 0.5)+ theme(legend.position = "none")
p_SuSiEx<-finemap_plot_fun(SuSiEx_plot_data, "SuSiEx", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_XMAP<-finemap_plot_fun(XMAP_plot_data, "XMAP", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_CARMAX<-finemap_plot_fun(CARMAX_plot_data, "CARMA-X", "PIP", 0.5)#+ drop_guides + theme(legend.position = "none")



# Gene Plot
plot.range <- c(min(candidate_region$POS), max(candidate_region$POS))
Gene_List_sub_coding<-Gene_List%>%filter(Chrom==paste0("chr",unique(candidate_region$CHR)))%>%filter(Start<max(candidate_region$POS),End>min(candidate_region$POS))%>%filter(Coding=="proteincoding")%>%filter(!is.na(cdsLength))%>%filter(GeneLength>10000)
p2<-gene_range_plot_fun(Gene_List_sub_coding,plot.range)

##Combine Plot together
combined_plot<-(p_EUR/p_MESuSiE/p_Paintor/p2+plot_layout(heights = c(1,1,1,1))|p_AFR/p_MFD/p_SuSiE_merged/p2+plot_layout(heights = c(1,1,1,1)))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")
p_ASE = combined_plot
#ggsave(paste0(plot_dir,"Region_",LD_BLOCK ,"_",trait,".pdf"),combined_plot,width = 180, height = 210,units="mm",dpi=300)


combined_plot_supp = (p_meta/p_EUR_SuSiE/p_XMAP+plot_layout(heights = c(1,1,1))|p_SuSiEx/p_AFR_SuSiE/p_CARMAX+plot_layout(heights = c(1,1,1)))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")
ggsave(paste0(plot_dir,"Supp_Region_",LD_BLOCK ,"_",trait,".pdf"),combined_plot_supp,width = 210, height = 180,units="mm",dpi=500)


###################################################################
#
#   Combine two locus plots
#  
##################################################################


p_ASV = p_ASV +custom_theme()
p_ASE = p_ASE +custom_theme()


library(patchwork)
library(cowplot) # Required for get_legend

# --- Step 1: Extract the Legend ---
# Extract the legend from p_ASV (since p_ASE has no legend)
# Format the legend horizontally before extraction
legend_source <- p_ASV + 
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal"
  )
my_legend <- get_legend(legend_source)

# --- Step 2: Clean and Wrap the Plots ---
# 1. Turn off legends in the plots (so they don't show up inside the grid)
p_ASV_clean <- p_ASV & theme(legend.position = "none")
p_ASE_clean <- p_ASE & theme(legend.position = "none")

# 2. WRAP them using wrap_elements()
# This forces Patchwork to treat the entire 8-panel grid as ONE solid block.
# This prevents the "squishing" or distortion of the plot columns.
p_ASV_fixed <- wrap_elements(p_ASV_clean)
p_ASE_fixed <- wrap_elements(p_ASE_clean)

# --- Step 3: Combine Layout ---
# Place ASV and ASE side-by-side (Row 1)
# Place the extracted Legend below them (Row 2)
p_out <- (p_ASE_fixed + p_ASV_fixed) / my_legend +
  plot_layout(
    widths = c(1, 1),      # Force 50/50 split for the plots
    heights = c(20, 1)     # Give the plots lots of space, legend small space
  )+
  plot_annotation(
    tag_levels = 'a',
    theme = theme(plot.tag = element_text(size = 20, face = "bold")) # Customize tag style
  )

# Print
p_out

ggsave(paste0(plot_dir,"Region_combined_locus_test.pdf"),p_out,width = 260, height = 180,units="mm",dpi=600)

#ggsave(paste0(plot_dir,"Region_combined_locus.pdf"),p_out,width = 400, height = 300,units="mm",dpi=600)

