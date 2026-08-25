#res_out_name<-c("shared_50_baseline.RData","1000G_baseline.RData","BB_50000_baseline.RData","Low_cor_baseline.RData","shared_0_baseline.RData","shared_50_External_baseline.RData","shared_100_alternative.RData","Meta_N_baseline.RData")
res_out_name<-c("shared_50_external_ld_updated_xmap_meta_susie_rss.RData")
library(ggplot2)
library(ggrepel)
library(grid)
library(egg)
library(dplyr)
library(forcats)
library(gridExtra)
library(patchwork)
library(ggpattern)
library(data.table)
library(ggpubr)
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/utility.R")
res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/res_summary/")
system(paste0("mkdir -p ",res_dir))
#for(res_out in res_out_name){

res_out = res_out_name
plot_dir<-paste0(res_dir,"Figure/",unlist(strsplit(res_out,".RData")),"/")
system(paste0("mkdir -p ",plot_dir))
########################################################################
#
#           50% causal variants are shared
#
########################################################################
###################
#
#Set Size&Power
#
###################

method_levels_base <- c("MESuSiE", "SuSiE", "SuSiE-weighted","SuSiE-merged", 
                        "Paintor", "MultiSuSiE", "SuSiEx", "XMAP", "CARMAX")
method_levels_base_2 <- c("MESuSiE", "SuSiE", "SuSiE_weighted","SuSiE_merged", 
                          "Paintor", "MultiSuSiE", "SuSiEx", "XMAP", "CARMAX")


load(paste0(res_dir,res_out))   
# 1. Clean Data & Rename Methods
all_Set_data_dataframe <- all_Set_data_dataframe %>%
  mutate(Method = recode(Method,
                         "SuSiE_meta_weighted" = "SuSiE-weighted",
                         "SuSiE_meta_merged"   = "SuSiE-merged"))

set_power_summary <- set_power_summary %>%
  mutate(Method = recode(Method,
                         "SuSiE_meta_weighted" = "SuSiE-weighted",
                         "SuSiE_meta_merged"   = "SuSiE-merged"))%>% mutate(Method  = factor(Method, levels = method_levels_base))

# 2. Calculate Upper Limit for Boxplot
upper_limit <- round(
  all_Set_data_dataframe %>%
    filter(Method == "Paintor", causal_num == "Num~Causal == 5 ", h2 == "~h^2 == 10^-4") %>%
    summarise(upper = quantile(Size, 0.75)) %>%
    pull(upper)
) + 50

# ---------------------------------------------------------
# 3. CREATE BOXPLOT (Top) - NO X AXIS, NO LEGEND
# ---------------------------------------------------------
p_size_box <- Set_Size_fun_order(
  all_Set_data_dataframe %>% mutate(Size = log2(Size + 1)),
  upper_limit = log2(upper_limit)
) +
  ylab("log2(Set Size + 1)") +
  # Remove X-axis labels/ticks and suppress legend completely
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  guides(fill = "none", color = "none") # Ensure this plot adds NOTHING to the legend

# ---------------------------------------------------------
# 4. CREATE BARPLOT (Bottom) - NO X AXIS (Rely on Legend)
# ---------------------------------------------------------
p_power_bar <- Set_Power_fun_legend(set_power_summary)

# Conditional Y-limit
if (unlist(strsplit(res_out, ".RData")) == "shared_100_alternative") {
  p_power_bar <- p_power_bar + ylim(c(0, 1.1))
}

# Remove X-axis labels/ticks (assuming color legend identifies method)
p_power_bar <- p_power_bar +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )

# ---------------------------------------------------------
# 5. COMBINE
# ---------------------------------------------------------
# Logic: Stack them, collect the legend (which now only comes from p_power_bar)
size_power <- (p_size_box / p_power_bar) +
  plot_layout(guides = "collect") & 
  plot_annotation(tag_levels = 'a') &
  theme(
    plot.tag = element_text(size = 7, face = "bold"),
    legend.position = "bottom"
  )

# Print to verify
print(size_power)


ggsave(paste0(plot_dir,"Set_size_power.pdf"),size_power,height=180, width=210, units = "mm",dpi=600)


##########################################################
#
# Either ancestry 
# Set Size&Power | ROC | FDR Power 
#
##########################################################

###################
#
#ROC
#
###################
either_all_ROC_data_dataframe<-either_all_ROC_data_dataframe%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","MultiSuSiE","SuSiEx","XMAP","CARMAX"))
p_ROC_Either<-ROC_shared_fun(either_all_ROC_data_dataframe)
###################
#
#FDR&Power
#
###################
FDR_Power_either$Method <- factor(
  FDR_Power_either$Method,
  levels = method_levels_base_2
)

power_upper_limit<-FDR_Power_either%>%filter(FDR!=0.5)%>%ungroup(Method,h2,causal_num)%>%summarise(upper_limit = min(ceiling(max(Power)*10)/10+0.1,1))%>%pull(upper_limit)
p_FDR_Power_Either<-FDR_Power_shared_fun_noannot(FDR_Power_either%>%filter(FDR!=0.5))+ylim(0,power_upper_limit)
#ROC_FDR_Power_Either <- (p_ROC_Either / p_FDR_Power_Either) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))

####################
#
#PIP calibration
#
####################
# ROC_FDR_Power_Either <- (p_ROC_Either / p_FDR_Power_Either) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
#  ROC_FDR_Power_Either<-ROC_FDR_Power_Either+ plot_layout(guides = "collect",heights = c(2, 1))&theme(legend.position = 'bottom')

PIP_calibration_either_byh2<-create_obs_frq_byh2(data_all%>%select(Signal,h2,causal_num,MESuSiE_Either,SuSiE_Either, SuSiE_weighted_PIP , SuSiE_merged_PIP , Paintor_PIP, SuSiEx_PIP, XMAP_PIP,MultiSuSiE_PIP,CARMAX_Either),c(1,2,3),c("MESuSiE_Either","SuSiE_Either","SuSiE_weighted_PIP","SuSiE_merged_PIP","Paintor_PIP","SuSiEx_PIP","XMAP_PIP","MultiSuSiE_PIP","CARMAX_Either"))
PIP_calibration_either_byh2<- PIP_calibration_either_byh2%>%mutate(Method = fct_recode(Method, "MESuSiE" = "MESuSiE_Either","SuSiE" = "SuSiE_Either","SuSiE_weighted" = "SuSiE_weighted_PIP","SuSiE_merged" = "SuSiE_merged_PIP","Paintor" = "Paintor_PIP", "SuSiEx" = "SuSiEx_PIP","XMAP" = "XMAP_PIP","MultiSuSiE" = "MultiSuSiE_PIP","CARMAX" = "CARMAX_Either"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor", "MultiSuSiE","SuSiEx","XMAP","CARMAX"))

p_calibration_byh2<-PIP_calibration_shared_byh2_fun(PIP_calibration_either_byh2)

ROC_FDR_Power_Calibration_Either_Plot<-ggpubr::ggarrange(p_ROC_Either,p_FDR_Power_Either,p_calibration_byh2,nrow = 3,ncol=1,
                                                         common.legend = TRUE, legend="bottom",labels = c("a","b","c","d","e","f"),font.label=list(color="black",size=7))
ggsave(paste0(plot_dir,"ROC_FDR_Power_Calibration_Either.pdf"),ROC_FDR_Power_Calibration_Either_Plot,
       height = 180, width = 210,units = "mm",dpi=600)

##########################################################
#
# Shared Signal 
#  ROC | FDR Power | PIP calibration | PIP scatter plot 
#
##########################################################
###################
#
#ROC
#
###################
shared_all_ROC_data_dataframe<-shared_all_ROC_data_dataframe%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","MultiSuSiE","SuSiEx","XMAP","CARMAX"))
p_ROC_shared<-ROC_shared_fun(shared_all_ROC_data_dataframe)
###################
#
#FDR&Power
#
###################
FDR_Power_shared$Method <- factor(
  FDR_Power_shared$Method,
  levels = method_levels_base_2
)
power_upper_limit<-FDR_Power_shared%>%filter(FDR!=0.5)%>%ungroup(Method,h2,causal_num)%>%summarise(upper_limit = min(ceiling(max(Power)*10)/10+0.1,1))%>%pull(upper_limit)
p_FDR_Power_shared<-FDR_Power_shared_fun_noannot(FDR_Power_shared%>%filter(FDR!=0.5))+ylim(0,power_upper_limit)

#  p_FDR_Power_shared<-FDR_Power_shared_fun(FDR_Power_shared)+ylim(0, min(ceiling(max(FDR_Power_shared$Power)*10)/10+0.1,1))
#  ROC_FDR_Power_shared <- (p_ROC_shared / p_FDR_Power_shared) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
#  ROC_FDR_Power_shared<-ROC_FDR_Power_shared+ plot_layout(guides = "collect",heights = c(2, 1))&theme(legend.position = 'bottom')
#  ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_Shared.pdf"),ROC_FDR_Power_shared,height=180, width=180, units = "mm",dpi=500)

####################
#
#PIP calibration
#
####################


PIP_calibration_shared_byh2<-create_obs_frq_byh2(data_all%>%select(Signal,h2,causal_num,MESuSiE_Shared,SuSiE_Shared, SuSiE_weighted_PIP , SuSiE_merged_PIP , Paintor_PIP, SuSiEx_PIP, XMAP_PIP,MultiSuSiE_PIP,CARMAX_Shared),c(3),c("MESuSiE_Shared","SuSiE_Shared","SuSiE_weighted_PIP","SuSiE_merged_PIP","Paintor_PIP","SuSiEx_PIP","XMAP_PIP","MultiSuSiE_PIP","CARMAX_Shared"))
PIP_calibration_shared_byh2<- PIP_calibration_shared_byh2%>%mutate(Method = fct_recode(Method, "MESuSiE" = "MESuSiE_Shared","SuSiE" = "SuSiE_Shared","SuSiE_weighted" = "SuSiE_weighted_PIP","SuSiE_merged" = "SuSiE_merged_PIP","Paintor" = "Paintor_PIP", "SuSiEx" = "SuSiEx_PIP","XMAP" = "XMAP_PIP","MultiSuSiE" = "MultiSuSiE_PIP","CARMAX" = "CARMAX_Shared"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","MultiSuSiE", "SuSiEx","XMAP","CARMAX"))


p_calibration_shared_byh2<-PIP_calibration_shared_byh2_fun(PIP_calibration_shared_byh2)

ROC_FDR_Power_Calibration_shared_Plot<-ggpubr::ggarrange(p_ROC_shared,p_FDR_Power_shared,p_calibration_shared_byh2,nrow = 3,ncol=1,
                                                         common.legend = TRUE, legend="bottom",labels = c("a","b","c","d","e","f"),font.label=list(color="black",size=7))
ggsave(paste0(plot_dir,"ROC_FDR_Power_Calibration_shared.pdf"),ROC_FDR_Power_Calibration_shared_Plot,
       height = 180, width = 210,units = "mm",dpi=600)


# p_calibration_shared<-PIP_calibration_shared_fun(PIP_calibration_Shared)
# ggsave(paste0(plot_dir,"PIP_Calibration_Shared.pdf"),p_calibration_shared,height = 121, width = 121,units = "mm",dpi=500)

##########################################################
#
# Ancestry-specific Signal 
#  ROC | FDR Power | PIP calibration  
#
##########################################################
###################
#
#ROC
#
###################

ancestry_all_ROC_data_dataframe <- ancestry_all_ROC_data_dataframe %>% 
  mutate(Method = as.character(Method))

split_list <- strsplit(ancestry_all_ROC_data_dataframe %>% pull(Method), " +")

ancestry_all_ROC_data_dataframe <- ancestry_all_ROC_data_dataframe %>%mutate(
  Method = sapply(split_list, `[`, 1),
  Ancestry = sapply(split_list, `[`, 2)
)%>%mutate(Method = fct_relevel(Method,c(
  "MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged",
  "Paintor","MultiSuSiE","SuSiEx","XMAP","CARMAX"
)),Ancestry = fct_relevel(Ancestry, "WB","BB"))%>%mutate(Ancestry = fct_recode(Ancestry, "White British" = "WB"   , "Black British" = "BB" ))


p_ROC_ancestry<-ROC_ancestry_fun(ancestry_all_ROC_data_dataframe)

###################
#
#FDR&Power
#
###################
FDR_Power_ancestry <- FDR_Power_ancestry %>% mutate(Method = as.character(Method))
split_list <- strsplit(FDR_Power_ancestry %>% pull(Method), " +")
FDR_Power_ancestry <- FDR_Power_ancestry%>%ungroup(h2,causal_num,Method) %>%mutate(
  Method = sapply(split_list, `[`, 1),
  Ancestry = sapply(split_list, `[`, 2)
)%>%mutate(Method = fct_relevel(Method, "MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","MultiSuSiE","SuSiEx","XMAP","CARMAX"),Ancestry = fct_relevel(Ancestry, "WB","BB"))%>%mutate(Ancestry = fct_recode(Ancestry, "White British" = "WB"   , "Black British" = "BB" ))

FDR_Power_ancestry<-FDR_Power_ancestry%>%filter(FDR!=0.5)
power_upper_limit<-FDR_Power_ancestry%>%summarise(upper_limit = min(ceiling(max(Power)*10)/10+0.1,1))%>%pull(upper_limit)

dyn.load("/apps/spack/negishi/apps/curl/7.85.0-gcc-12.2.0-4phaxqw/lib/libcurl.so.4",
         local = FALSE)

dyn.load("/apps/spack/negishi/apps/r-sf/1.0-9-gcc-12.2.0-sftozdw/rlib/R/library/sf/libs/sf.so")

p_FDR_Power_ancestry<-FDR_Power_ancestry_fun(FDR_Power_ancestry)+ylim(0, power_upper_limit)
ROC_FDR_Power_ancestry<- (p_ROC_ancestry / p_FDR_Power_ancestry) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
ROC_FDR_Power_ancestry<-ROC_FDR_Power_ancestry+ plot_layout(heights = c(1, 1))
ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_Ancestry.pdf"),ROC_FDR_Power_ancestry,height=180, width=210, units = "mm",dpi=600)


###################
#
#PIP calibration
#
###################   

PIP_calibration_ancestry<- PIP_calibration_ancestry%>%group_by(causal_num)%>%mutate(Method = fct_recode(Method, "MESuSiE White British" = "MESuSiE~WB", "MESuSiE Black British" = "MESuSiE~BB","SuSiE White British" = "SuSiE~WB","SuSiE Black British" = "SuSiE~BB","SuSiE_weighted White British" = "SuSiE_weighted~WB","SuSiE_weighted Black British" = "SuSiE_weighted~BB","SuSiE_merged White British" = "SuSiE_merged~WB","SuSiE_merged Black British" = "SuSiE_merged~BB","Paintor White British" = "Paintor~WB","Paintor Black British" = "Paintor~BB","SuSiEx White British" = "SuSiEx~WB","SuSiEx Black British" = "SuSiEx~BB","XMAP White British" = "XMAP~WB","XMAP Black British" = "XMAP~BB","MultiSuSiE White British" = "MultiSuSiE~WB","MultiSuSiE Black British" = "MultiSuSiE~BB","CARMAX White British" = "CARMAX~WB","CARMAX Black British" = "CARMAX~BB"))
levels(PIP_calibration_ancestry$Method)<-c(paste0("MESuSiE~","White~","British"),paste0("MESuSiE~","Black~","British"),paste0("SuSiE~","White~","British"),paste0("SuSiE~","Black~","British"),paste0("SuSiE_weighted~","White~","British"),paste0("SuSiE_weighted~","Black~","British"),paste0("SuSiE_merged~","White~","British"),paste0("SuSiE_merged~","Black~","British"),paste0("Paintor~","White~","British"),paste0("Paintor~","Black~","British"),paste0("MultiSuSiE~","White~","British"),paste0("MultiSuSiE~","Black~","British"),paste0("SuSiEx~","White~","British"),paste0("SuSiEx~","Black~","British"),paste0("XMAP~","White~","British"),paste0("XMAP~","Black~","British"),paste0("CARMAX~","White~","British"),paste0("CARMAX~","Black~","British"))
PIP_calibration_ancestry_WB = PIP_calibration_ancestry %>% 
  filter(grepl("White~British", Method))

# Filter for Black British methods
PIP_calibration_ancestry_BB = PIP_calibration_ancestry %>% 
  filter(grepl("Black~British", Method))

p_calibration_ancestry_WB<-PIP_calibration_ancestry_fun(PIP_calibration_ancestry_WB)
p_calibration_ancestry_BB<-PIP_calibration_ancestry_fun(PIP_calibration_ancestry_BB)
p_out = (p_calibration_ancestry_WB+p_calibration_ancestry_BB)+plot_annotation(tag_levels = 'a') &
  theme(
    plot.tag = element_text(size = 7, face = "bold"))
ggsave(paste0(plot_dir,"PIP_Calibration_Ancestry.pdf"),p_out,height = 180, width = 210,units = "mm",dpi=500)
#}      



####################################SLALOM##################################################
#res_out_name<-c("shared_50_baseline.RData","1000G_baseline.RData","BB_50000_baseline.RData","Low_cor_baseline.RData","shared_0_baseline.RData","shared_50_External_baseline.RData","shared_100_alternative.RData","Meta_N_baseline.RData")
res_out_name<-c("shared_50_external_ld_updated_xmap_meta_slalom.RData")
library(ggplot2)
library(ggrepel)
library(grid)
library(egg)
library(dplyr)
library(forcats)
library(gridExtra)
library(patchwork)
library(ggpattern)
library(data.table)
library(ggpubr)
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/utility.R")
res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/res_summary/")
system(paste0("mkdir -p ",res_dir))
#for(res_out in res_out_name){

res_out = res_out_name
plot_dir<-paste0(res_dir,"Figure/",unlist(strsplit(res_out,".RData")),"/")
system(paste0("mkdir -p ",plot_dir))
########################################################################
#
#           50% causal variants are shared
#
########################################################################
###################
#
#Set Size&Power
#
###################

method_levels_base <- c("MESuSiE", "SuSiE", "SuSiE-weighted","SuSiE-merged", 
                        "Paintor", "MultiSuSiE", "SuSiEx", "XMAP", "CARMAX")
method_levels_base_2 <- c("MESuSiE", "SuSiE", "SuSiE_weighted","SuSiE_merged", 
                          "Paintor", "MultiSuSiE", "SuSiEx", "XMAP", "CARMAX")


load(paste0(res_dir,res_out))   
# 1. Clean Data & Rename Methods
all_Set_data_dataframe <- all_Set_data_dataframe %>%
  mutate(Method = recode(Method,
                         "SuSiE_meta_weighted" = "SuSiE-weighted",
                         "SuSiE_meta_merged"   = "SuSiE-merged"))

set_power_summary <- set_power_summary %>%
  mutate(Method = recode(Method,
                         "SuSiE_meta_weighted" = "SuSiE-weighted",
                         "SuSiE_meta_merged"   = "SuSiE-merged"))%>% mutate(Method  = factor(Method, levels = method_levels_base))

# 2. Calculate Upper Limit for Boxplot
upper_limit <- round(
  all_Set_data_dataframe %>%
    filter(Method == "Paintor", causal_num == "Num~Causal == 5 ", h2 == "~h^2 == 10^-4") %>%
    summarise(upper = quantile(Size, 0.75)) %>%
    pull(upper)
) + 50

# ---------------------------------------------------------
# 3. CREATE BOXPLOT (Top) - NO X AXIS, NO LEGEND
# ---------------------------------------------------------
p_size_box <- Set_Size_fun_order(
  all_Set_data_dataframe %>% mutate(Size = log2(Size + 1)),
  upper_limit = log2(upper_limit)
) +
  ylab("log2(Set Size + 1)") +
  # Remove X-axis labels/ticks and suppress legend completely
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  guides(fill = "none", color = "none") # Ensure this plot adds NOTHING to the legend

# ---------------------------------------------------------
# 4. CREATE BARPLOT (Bottom) - NO X AXIS (Rely on Legend)
# ---------------------------------------------------------
p_power_bar <- Set_Power_fun_legend(set_power_summary)

# Conditional Y-limit
if (unlist(strsplit(res_out, ".RData")) == "shared_100_alternative") {
  p_power_bar <- p_power_bar + ylim(c(0, 1.1))
}

# Remove X-axis labels/ticks (assuming color legend identifies method)
p_power_bar <- p_power_bar +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )

# ---------------------------------------------------------
# 5. COMBINE
# ---------------------------------------------------------
# Logic: Stack them, collect the legend (which now only comes from p_power_bar)
size_power <- (p_size_box / p_power_bar) +
  plot_layout(guides = "collect") & 
  plot_annotation(tag_levels = 'a') &
  theme(
    plot.tag = element_text(size = 7, face = "bold"),
    legend.position = "bottom"
  )

# Print to verify
print(size_power)


ggsave(paste0(plot_dir,"Set_size_power.pdf"),size_power,height=180, width=210, units = "mm",dpi=600)


##########################################################
#
# Either ancestry 
# Set Size&Power | ROC | FDR Power 
#
##########################################################

###################
#
#ROC
#
###################
either_all_ROC_data_dataframe<-either_all_ROC_data_dataframe%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","MultiSuSiE","SuSiEx","XMAP","CARMAX"))
p_ROC_Either<-ROC_shared_fun(either_all_ROC_data_dataframe)
###################
#
#FDR&Power
#
###################
FDR_Power_either$Method <- factor(
  FDR_Power_either$Method,
  levels = method_levels_base_2
)

power_upper_limit<-FDR_Power_either%>%filter(FDR!=0.5)%>%ungroup(Method,h2,causal_num)%>%summarise(upper_limit = min(ceiling(max(Power)*10)/10+0.1,1))%>%pull(upper_limit)
p_FDR_Power_Either<-FDR_Power_shared_fun_noannot(FDR_Power_either%>%filter(FDR!=0.5))+ylim(0,power_upper_limit)
#ROC_FDR_Power_Either <- (p_ROC_Either / p_FDR_Power_Either) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))

####################
#
#PIP calibration
#
####################
# ROC_FDR_Power_Either <- (p_ROC_Either / p_FDR_Power_Either) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
#  ROC_FDR_Power_Either<-ROC_FDR_Power_Either+ plot_layout(guides = "collect",heights = c(2, 1))&theme(legend.position = 'bottom')

PIP_calibration_either_byh2<-create_obs_frq_byh2(data_all%>%select(Signal,h2,causal_num,MESuSiE_Either,SuSiE_Either, SuSiE_weighted_PIP , SuSiE_merged_PIP , Paintor_PIP, SuSiEx_PIP, XMAP_PIP,MultiSuSiE_PIP,CARMAX_Either),c(1,2,3),c("MESuSiE_Either","SuSiE_Either","SuSiE_weighted_PIP","SuSiE_merged_PIP","Paintor_PIP","SuSiEx_PIP","XMAP_PIP","MultiSuSiE_PIP","CARMAX_Either"))
PIP_calibration_either_byh2<- PIP_calibration_either_byh2%>%mutate(Method = fct_recode(Method, "MESuSiE" = "MESuSiE_Either","SuSiE" = "SuSiE_Either","SuSiE_weighted" = "SuSiE_weighted_PIP","SuSiE_merged" = "SuSiE_merged_PIP","Paintor" = "Paintor_PIP", "SuSiEx" = "SuSiEx_PIP","XMAP" = "XMAP_PIP","MultiSuSiE" = "MultiSuSiE_PIP","CARMAX" = "CARMAX_Either"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor", "MultiSuSiE","SuSiEx","XMAP","CARMAX"))

p_calibration_byh2<-PIP_calibration_shared_byh2_fun(PIP_calibration_either_byh2)

ROC_FDR_Power_Calibration_Either_Plot<-ggpubr::ggarrange(p_ROC_Either,p_FDR_Power_Either,p_calibration_byh2,nrow = 3,ncol=1,
                                                         common.legend = TRUE, legend="bottom",labels = c("a","b","c","d","e","f"),font.label=list(color="black",size=7))
ggsave(paste0(plot_dir,"ROC_FDR_Power_Calibration_Either.pdf"),ROC_FDR_Power_Calibration_Either_Plot,
       height = 180, width = 210,units = "mm",dpi=600)

##########################################################
#
# Shared Signal 
#  ROC | FDR Power | PIP calibration | PIP scatter plot 
#
##########################################################
###################
#
#ROC
#
###################
shared_all_ROC_data_dataframe<-shared_all_ROC_data_dataframe%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","MultiSuSiE","SuSiEx","XMAP","CARMAX"))
p_ROC_shared<-ROC_shared_fun(shared_all_ROC_data_dataframe)
###################
#
#FDR&Power
#
###################
FDR_Power_shared$Method <- factor(
  FDR_Power_shared$Method,
  levels = method_levels_base_2
)
power_upper_limit<-FDR_Power_shared%>%filter(FDR!=0.5)%>%ungroup(Method,h2,causal_num)%>%summarise(upper_limit = min(ceiling(max(Power)*10)/10+0.1,1))%>%pull(upper_limit)
p_FDR_Power_shared<-FDR_Power_shared_fun_noannot(FDR_Power_shared%>%filter(FDR!=0.5))+ylim(0,power_upper_limit)

#  p_FDR_Power_shared<-FDR_Power_shared_fun(FDR_Power_shared)+ylim(0, min(ceiling(max(FDR_Power_shared$Power)*10)/10+0.1,1))
#  ROC_FDR_Power_shared <- (p_ROC_shared / p_FDR_Power_shared) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
#  ROC_FDR_Power_shared<-ROC_FDR_Power_shared+ plot_layout(guides = "collect",heights = c(2, 1))&theme(legend.position = 'bottom')
#  ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_Shared.pdf"),ROC_FDR_Power_shared,height=180, width=180, units = "mm",dpi=500)

####################
#
#PIP calibration
#
####################


PIP_calibration_shared_byh2<-create_obs_frq_byh2(data_all%>%select(Signal,h2,causal_num,MESuSiE_Shared,SuSiE_Shared, SuSiE_weighted_PIP , SuSiE_merged_PIP , Paintor_PIP, SuSiEx_PIP, XMAP_PIP,MultiSuSiE_PIP,CARMAX_Shared),c(3),c("MESuSiE_Shared","SuSiE_Shared","SuSiE_weighted_PIP","SuSiE_merged_PIP","Paintor_PIP","SuSiEx_PIP","XMAP_PIP","MultiSuSiE_PIP","CARMAX_Shared"))
PIP_calibration_shared_byh2<- PIP_calibration_shared_byh2%>%mutate(Method = fct_recode(Method, "MESuSiE" = "MESuSiE_Shared","SuSiE" = "SuSiE_Shared","SuSiE_weighted" = "SuSiE_weighted_PIP","SuSiE_merged" = "SuSiE_merged_PIP","Paintor" = "Paintor_PIP", "SuSiEx" = "SuSiEx_PIP","XMAP" = "XMAP_PIP","MultiSuSiE" = "MultiSuSiE_PIP","CARMAX" = "CARMAX_Shared"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","MultiSuSiE", "SuSiEx","XMAP","CARMAX"))


p_calibration_shared_byh2<-PIP_calibration_shared_byh2_fun(PIP_calibration_shared_byh2)

ROC_FDR_Power_Calibration_shared_Plot<-ggpubr::ggarrange(p_ROC_shared,p_FDR_Power_shared,p_calibration_shared_byh2,nrow = 3,ncol=1,
                                                         common.legend = TRUE, legend="bottom",labels = c("a","b","c","d","e","f"),font.label=list(color="black",size=7))
ggsave(paste0(plot_dir,"ROC_FDR_Power_Calibration_shared.pdf"),ROC_FDR_Power_Calibration_shared_Plot,
       height = 180, width = 210,units = "mm",dpi=600)


# p_calibration_shared<-PIP_calibration_shared_fun(PIP_calibration_Shared)
# ggsave(paste0(plot_dir,"PIP_Calibration_Shared.pdf"),p_calibration_shared,height = 121, width = 121,units = "mm",dpi=500)

##########################################################
#
# Ancestry-specific Signal 
#  ROC | FDR Power | PIP calibration  
#
##########################################################
###################
#
#ROC
#
###################

ancestry_all_ROC_data_dataframe <- ancestry_all_ROC_data_dataframe %>% 
  mutate(Method = as.character(Method))

split_list <- strsplit(ancestry_all_ROC_data_dataframe %>% pull(Method), " +")

ancestry_all_ROC_data_dataframe <- ancestry_all_ROC_data_dataframe %>%mutate(
  Method = sapply(split_list, `[`, 1),
  Ancestry = sapply(split_list, `[`, 2)
)%>%mutate(Method = fct_relevel(Method,c(
  "MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged",
  "Paintor","MultiSuSiE","SuSiEx","XMAP","CARMAX"
)),Ancestry = fct_relevel(Ancestry, "WB","BB"))%>%mutate(Ancestry = fct_recode(Ancestry, "White British" = "WB"   , "Black British" = "BB" ))


p_ROC_ancestry<-ROC_ancestry_fun(ancestry_all_ROC_data_dataframe)

###################
#
#FDR&Power
#
###################
FDR_Power_ancestry <- FDR_Power_ancestry %>% mutate(Method = as.character(Method))
split_list <- strsplit(FDR_Power_ancestry %>% pull(Method), " +")
FDR_Power_ancestry <- FDR_Power_ancestry%>%ungroup(h2,causal_num,Method) %>%mutate(
  Method = sapply(split_list, `[`, 1),
  Ancestry = sapply(split_list, `[`, 2)
)%>%mutate(Method = fct_relevel(Method, "MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","MultiSuSiE","SuSiEx","XMAP","CARMAX"),Ancestry = fct_relevel(Ancestry, "WB","BB"))%>%mutate(Ancestry = fct_recode(Ancestry, "White British" = "WB"   , "Black British" = "BB" ))

FDR_Power_ancestry<-FDR_Power_ancestry%>%filter(FDR!=0.5)
power_upper_limit<-FDR_Power_ancestry%>%summarise(upper_limit = min(ceiling(max(Power)*10)/10+0.1,1))%>%pull(upper_limit)

dyn.load("/apps/spack/negishi/apps/curl/7.85.0-gcc-12.2.0-4phaxqw/lib/libcurl.so.4",
         local = FALSE)

dyn.load("/apps/spack/negishi/apps/r-sf/1.0-9-gcc-12.2.0-sftozdw/rlib/R/library/sf/libs/sf.so")

p_FDR_Power_ancestry<-FDR_Power_ancestry_fun(FDR_Power_ancestry)+ylim(0, power_upper_limit)
ROC_FDR_Power_ancestry<- (p_ROC_ancestry / p_FDR_Power_ancestry) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
ROC_FDR_Power_ancestry<-ROC_FDR_Power_ancestry+ plot_layout(heights = c(1, 1))
ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_Ancestry.pdf"),ROC_FDR_Power_ancestry,height=180, width=210, units = "mm",dpi=600)


###################
#
#PIP calibration
#
###################   

PIP_calibration_ancestry<- PIP_calibration_ancestry%>%group_by(causal_num)%>%mutate(Method = fct_recode(Method, "MESuSiE White British" = "MESuSiE~WB", "MESuSiE Black British" = "MESuSiE~BB","SuSiE White British" = "SuSiE~WB","SuSiE Black British" = "SuSiE~BB","SuSiE_weighted White British" = "SuSiE_weighted~WB","SuSiE_weighted Black British" = "SuSiE_weighted~BB","SuSiE_merged White British" = "SuSiE_merged~WB","SuSiE_merged Black British" = "SuSiE_merged~BB","Paintor White British" = "Paintor~WB","Paintor Black British" = "Paintor~BB","SuSiEx White British" = "SuSiEx~WB","SuSiEx Black British" = "SuSiEx~BB","XMAP White British" = "XMAP~WB","XMAP Black British" = "XMAP~BB","MultiSuSiE White British" = "MultiSuSiE~WB","MultiSuSiE Black British" = "MultiSuSiE~BB","CARMAX White British" = "CARMAX~WB","CARMAX Black British" = "CARMAX~BB"))
levels(PIP_calibration_ancestry$Method)<-c(paste0("MESuSiE~","White~","British"),paste0("MESuSiE~","Black~","British"),paste0("SuSiE~","White~","British"),paste0("SuSiE~","Black~","British"),paste0("SuSiE_weighted~","White~","British"),paste0("SuSiE_weighted~","Black~","British"),paste0("SuSiE_merged~","White~","British"),paste0("SuSiE_merged~","Black~","British"),paste0("Paintor~","White~","British"),paste0("Paintor~","Black~","British"),paste0("MultiSuSiE~","White~","British"),paste0("MultiSuSiE~","Black~","British"),paste0("SuSiEx~","White~","British"),paste0("SuSiEx~","Black~","British"),paste0("XMAP~","White~","British"),paste0("XMAP~","Black~","British"),paste0("CARMAX~","White~","British"),paste0("CARMAX~","Black~","British"))
PIP_calibration_ancestry_WB = PIP_calibration_ancestry %>% 
  filter(grepl("White~British", Method))

# Filter for Black British methods
PIP_calibration_ancestry_BB = PIP_calibration_ancestry %>% 
  filter(grepl("Black~British", Method))

p_calibration_ancestry_WB<-PIP_calibration_ancestry_fun(PIP_calibration_ancestry_WB)
p_calibration_ancestry_BB<-PIP_calibration_ancestry_fun(PIP_calibration_ancestry_BB)
p_out = (p_calibration_ancestry_WB+p_calibration_ancestry_BB)+plot_annotation(tag_levels = 'a') &
  theme(
    plot.tag = element_text(size = 7, face = "bold"))
ggsave(paste0(plot_dir,"PIP_Calibration_Ancestry.pdf"),p_out,height = 180, width = 210,units = "mm",dpi=500)
#}      

#####################################################################################################################
#
#
#                     Compare result of BB = 50000 versus BB = 300,000/ In sample LD versus External LD
#
#
#####################################################################################################################

FDR_Power_Compare_fun<-function(FDR_Power,pattern_input){
  p = ggplot(FDR_Power, aes(FDR, Power,fill = Method,pattern = Pattern)) + 
    geom_col_pattern(position = "dodge", pattern_density = .01, pattern_spacing = 0.04,pattern_fill ="#DEDBD2",pattern_colour ="#d8c99b" ) +
    scale_fill_manual(name = "Method",values = c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","SuSiE_weighted"="#B2D3A4","SuSiE_merged"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAX"="#f2b56e"),guide=FALSE)+
    scale_pattern_manual(values=pattern_input)+
    facet_grid(vars(h2),vars(causal_num),labeller=label_parsed)+
    ylab("Power")+xlab("FDR")+
    theme_bw() + custom_theme() + theme(legend.position = "bottom")+
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2)+
    guides(pattern = guide_legend(override.aes = list(fill = "white", order = 1)))
  return(p)
}

ROC_ancestry_Compare_fun <- function(all_ROC_data_dataframe, pattern_input) {
  
  # 1. Define Base Colors
  base_colors <- c("MESuSiE"="#8da0cb", "SuSiE"="#66c2a5", "SuSiE_weighted"="#B2D3A4",
                   "SuSiE_merged"="#9FBA95", "Paintor"="#fc8d62", "MultiSuSiE"="#e78ac3",
                   "SuSiEx"="#E89DA0", "XMAP"="#ffd92f", "CARMAX"="#f2b56e", "MFD" = "#377eb8")
  
  # 2. Expand Colors for Ancestry Labels (WB/BB)
  # This ensures "MFD WB" and "MFD BB" get the same color as "MFD"
  ancestry_colors <- c(base_colors, 
                       setNames(base_colors, paste(names(base_colors), "WB")), 
                       setNames(base_colors, paste(names(base_colors), "BB")))
  
  p <- ggplot(all_ROC_data_dataframe, aes(x = Power, y = 1 - FDR)) + 
    
    # Lines for Methods
    geom_line(aes(linetype = Pattern, color = Method), size = 0.8) +
    
    # Points for Cutoff (PIP > 0.5)
    geom_point(data = all_ROC_data_dataframe %>% filter(Cutoff == 0.5), 
               aes(x = Power, y = 1 - FDR, color = Method, shape = 'PIP > 0.5'), 
               size = 3, 
               show.legend = c(color = FALSE)) + # Hide these points from the color legend line
    
    # Scales
    scale_color_manual(values = ancestry_colors) +
    scale_linetype_manual(values = pattern_input) +
    scale_shape_manual(name = ' ', values = c('PIP > 0.5' = 21)) + 
    
    # Axes
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) + 
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
    xlab("Recall") + 
    ylab("Precision") + 
    
    # Facets
    facet_grid(h2 ~ causal_num, labeller = label_parsed) + 
    
    # Theme & Layout
    theme_bw() + 
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(size = 7, face = "bold"),
      axis.text.y = element_text(size = 7),  
      axis.title.x = element_text(size = 7, face = "bold"),
      axis.title.y = element_text(size = 7, face = "bold"),
      strip.text.x = element_text(size = 7),
      strip.text.y = element_text(size = 7),
      strip.background = element_blank(),
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 5, face = "bold"),
      plot.title = element_text(size = 7, hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(), 
      axis.line = element_line(color = "black")
    ) + 
    
    # Border
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2) +
    
    # --- LEGEND CONFIGURATION ---
    guides(
      color = guide_legend(
        title = "Method",
        order = 1, 
        nrow = 2   # Display methods in two rows
      ),
      linetype = guide_legend(
        title = "Pattern",
        order = 2, 
        nrow = 2   # Align pattern legend height
      ), 
      shape = guide_legend(order = 3)
    )
  
  return(p)
}

# Function to load data and add a pattern column
load_data_with_pattern <- function(data_frame, pattern_name) {
  data_frame %>% mutate(Pattern = pattern_name)
}


source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/utility_missing.R")

res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/res_summary/")
plot_dir_name<-"shared_50_external_ld_control_compare"
res_out_1<-"shared_50_external_ld_updated_xmap_meta_susie_rss.RData"
res_out_2<-"shared_50_external_ld_updated_xmap_meta_slalom.RData"
pattern_name_1 = "Controlled by slalom"
pattern_name_2 = "Controlled by susie rss"
pattern_name_3 = "No Control"

Compare_Plot <- function(res_dir, res_out, pattern_name_1, pattern_name_2) {
  
  # Load the first dataset
  load(paste0(res_dir, res_out_1))  
  all_Set_data_dataframe_1<-load_data_with_pattern(all_Set_data_dataframe, pattern_name_1)
  set_power_summary_1<-load_data_with_pattern(set_power_summary, pattern_name_1)
  
  either_all_ROC_data_dataframe_1 <- load_data_with_pattern(either_all_ROC_data_dataframe, pattern_name_1)
  shared_all_ROC_data_dataframe_1 <- load_data_with_pattern(shared_all_ROC_data_dataframe, pattern_name_1)
  ancestry_all_ROC_data_dataframe_1 <- load_data_with_pattern(ancestry_all_ROC_data_dataframe, pattern_name_1)
  
  FDR_Power_either_1 <- load_data_with_pattern(FDR_Power_either, pattern_name_1)
  FDR_Power_shared_1 <- load_data_with_pattern(FDR_Power_shared, pattern_name_1)
  FDR_Power_ancestry_1 <- load_data_with_pattern(FDR_Power_ancestry, pattern_name_1)
  
  # Load the second dataset
  load(paste0(res_dir, res_out_2))  
  all_Set_data_dataframe_2<-load_data_with_pattern(all_Set_data_dataframe, pattern_name_2)
  set_power_summary_2<-load_data_with_pattern(set_power_summary, pattern_name_2)
  
  either_all_ROC_data_dataframe_2 <- load_data_with_pattern(either_all_ROC_data_dataframe, pattern_name_2)
  shared_all_ROC_data_dataframe_2 <- load_data_with_pattern(shared_all_ROC_data_dataframe, pattern_name_2)
  ancestry_all_ROC_data_dataframe_2 <- load_data_with_pattern(ancestry_all_ROC_data_dataframe, pattern_name_2)
  
  FDR_Power_either_2 <- load_data_with_pattern(FDR_Power_either, pattern_name_2)
  FDR_Power_shared_2 <- load_data_with_pattern(FDR_Power_shared, pattern_name_2)
  FDR_Power_ancestry_2 <- load_data_with_pattern(FDR_Power_ancestry, pattern_name_2)
  
  # Load the third dataset
  
  load(paste0(res_dir, "shared_50_external_ld_updated_xmap_meta.RData")) 
  all_Set_data_dataframe_3<-load_data_with_pattern(all_Set_data_dataframe, pattern_name_3)
  set_power_summary_3<-load_data_with_pattern(set_power_summary, pattern_name_3)
  
  either_all_ROC_data_dataframe_3 <- load_data_with_pattern(either_all_ROC_data_dataframe, pattern_name_3)
  shared_all_ROC_data_dataframe_3 <- load_data_with_pattern(shared_all_ROC_data_dataframe, pattern_name_3)
  ancestry_all_ROC_data_dataframe_3 <- load_data_with_pattern(ancestry_all_ROC_data_dataframe, pattern_name_3)
  
  FDR_Power_either_3 <- load_data_with_pattern(FDR_Power_either, pattern_name_3)
  FDR_Power_shared_3 <- load_data_with_pattern(FDR_Power_shared, pattern_name_3)
  FDR_Power_ancestry_3 <- load_data_with_pattern(FDR_Power_ancestry, pattern_name_3)
  
  # Combine Data Together
  Set_Data <- rbind(all_Set_data_dataframe_1,all_Set_data_dataframe_2,all_Set_data_dataframe_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  Set_Power_Data<-rbind(set_power_summary_1,set_power_summary_2,set_power_summary_3)%>%mutate(Pattern = factor(Pattern, levels = c(pattern_name_1, pattern_name_2, pattern_name_3)))
  
  ROC_either <- rbind(either_all_ROC_data_dataframe_1, either_all_ROC_data_dataframe_2,either_all_ROC_data_dataframe_3)
  ROC_shared <- rbind(shared_all_ROC_data_dataframe_1, shared_all_ROC_data_dataframe_2,shared_all_ROC_data_dataframe_3)
  ROC_ancestry <- rbind(ancestry_all_ROC_data_dataframe_1, ancestry_all_ROC_data_dataframe_2,ancestry_all_ROC_data_dataframe_3)
  
  FDR_Power_either <- rbind(FDR_Power_either_1, FDR_Power_either_2,FDR_Power_either_3)
  FDR_Power_shared <- rbind(FDR_Power_shared_1, FDR_Power_shared_2,FDR_Power_shared_3)
  FDR_Power_ancestry <- rbind(FDR_Power_ancestry_1, FDR_Power_ancestry_2,FDR_Power_ancestry_3)
  
  
  library(openxlsx)
  library(stringr)
  
  res_dir_table <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/supp_tables/"
  res_out_xlsx <- "External_LD_control_compare.xlsx"
  full_path <- paste0(res_dir_table, res_out_xlsx)
  sheet_name_raw <- str_remove(res_out_xlsx, ".xlsx")
  sheet_name <- substr(sheet_name_raw, 1, 31) 
  
  # 1. Define ALL Orders (Method, Pattern, Signal Type)
  method_order_base <- c("MFD", "MESuSiE", "SuSiE", "SuSiE_weighted", "SuSiE_merged", 
                         "Paintor", "MultiSuSiE", "SuSiEx", "XMAP", "CARMAX")
  
  # Generate Ancestry-Specific Levels (e.g., "MESuSiE WB", "MESuSiE BB")
  method_order_ancestry <- unlist(lapply(method_order_base, function(m) paste0(m, c(" WB", " BB"))))
  master_method_levels <- c(method_order_base, method_order_ancestry)
  
  pattern_order <- c("No Control", "Controlled by slalom", "Controlled by susie rss")
  signal_type_order <- c("Either", "Shared", "Ancestry")
  
  # 2. Process and Sort the Data
  combined_data <- bind_rows(
    "Either"   = FDR_Power_either,
    "Shared"   = FDR_Power_shared,
    "Ancestry" = FDR_Power_ancestry,
    .id = "Signal_Type" 
  ) %>%
    select(Method, h2, causal_num, Signal_Type, Pattern,FDR, Power) %>%
    
    # Apply Factors for ALL sorting columns
    mutate(
      Method = factor(Method, levels = master_method_levels),
      Pattern = factor(Pattern, levels = pattern_order),
      Signal_Type = factor(Signal_Type, levels = signal_type_order) # Enforces Either -> Shared -> Ancestry
    ) %>%
    
    # Order methods for reporting
    arrange(Signal_Type, Pattern, Method)
  
  # 3. View the result
  print(head(combined_data))
  
  
  wb <- createWorkbook()
  addWorksheet(wb, sheet_name)
  
  # Create a bold style for the Title
  title_style <- createStyle(textDecoration = "bold", fontSize = 12)
  
  # Write the Main Title at Row 1
  writeData(wb, sheet_name, "Supplementary Table: FDR Controlled Power Analysis (Combined)", startRow = 1)
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1)
  
  # Write the Combined Data Table starting at Row 2
  writeData(wb, sheet_name, combined_data, startRow = 2)
  
  # 5. Save
  saveWorkbook(wb, full_path, overwrite = TRUE)
  
  print(paste("Saved combined file to:", full_path))
  print(paste("Saved combined table to:", full_path))
  
  
  
  
  
  
  # Define a directory to save the plots
  plot_dir <- paste0(res_dir, "Figure/",plot_dir_name,"/")
  system(paste0("mkdir -p ",plot_dir))
  
  # Define patterns for plotting
  pattern_input_Power <- c("stripe","circle", "none")
  
  names(pattern_input_Power) <- c(pattern_name_1, pattern_name_2, pattern_name_3)
  
  pattern_input_line <- c("dashed","dotted", "solid")
  names(pattern_input_line) <- c(pattern_name_1, pattern_name_2, pattern_name_3)
  
  
  #####################################################################################
  rename_map <- c(
    "SuSiE_weighted" = "SuSiE_meta_weighted",
    "SuSiE_merged"   = "SuSiE_meta_merged"
  )
  
  Set_Power_Data <- Set_Power_Data %>%
    mutate(Method = fct_recode(Method, !!!rename_map))
  
  Set_Data <- Set_Data %>%
    mutate(Method = fct_recode(Method, !!!rename_map))
  
  Set_Data$Pattern<-factor(Set_Data$Pattern,levels = c(pattern_name_1,pattern_name_2,pattern_name_3))
  p_size_box<-Set_Size_Compare_fun(Set_Data%>%mutate(Size = log2(Size+1)),upper_limit = log2(max(Set_Data$Size))+1,pattern_input_Power)
  p_size_box<-p_size_box+ ylab("log2(Set Size + 1)")
  p_size_box <- p_size_box + 
    theme(
      axis.text.x = element_blank(),  # Removes the label text (e.g., Method names)
      axis.ticks.x = element_blank(), # Removes the tick marks
      axis.title.x = element_blank()  # Optional: Removes the axis title "Method" if present
    )
  Set_Power_Data$Pattern<-factor(Set_Power_Data$Pattern,levels = c(pattern_name_1,pattern_name_2,pattern_name_3))
  p_power_bar<-Set_Power_Compare_fun_supp(Set_Power_Data,pattern_input_Power)
  p_power_bar <- p_power_bar + 
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
  p_size_power <- (p_size_box / p_power_bar) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  p_size_power<-p_size_power+ plot_layout(guides = "collect",heights=c(2,1))&theme(legend.position = 'bottom')
  
  ggsave(paste0(plot_dir,"Set_Size_Power_Compare.pdf"),p_size_power,height=180, width=210, units = "mm",dpi=600)
  
  # Function to generate FDR power comparison plot (defined elsewhere)
  p_FDR_Power_either <- FDR_Power_Compare_fun(FDR_Power_either, pattern_input_Power)
  
  # Function to generate ROC ancestry comparison plot (defined elsewhere)
  p_ROC_either <- ROC_ancestry_Compare_fun(ROC_either, pattern_input_line)
  
  # Create and save the combined plot for either PIP
  p_ROC_FDR_Power_either <- (p_ROC_either + p_FDR_Power_either) + 
    plot_annotation(tag_levels = 'a') & 
    theme(plot.tag = element_text(size = 7, face = "bold"))
  p_ROC_FDR_Power_either <- p_ROC_FDR_Power_either + plot_layout(heights = c(2, 1))
  
  ggsave(paste0(plot_dir, "PIP_ROC_FDR_Power_Either_Compare.pdf"), p_ROC_FDR_Power_either, height = 180, width = 210, units = "mm", dpi = 500)
  
  ###Plot for Shared PIP
  
  p_FDR_Power_shared<-FDR_Power_Compare_fun(FDR_Power_shared,pattern_input_Power)
  p_ROC_shared<-ROC_ancestry_Compare_fun(ROC_shared,pattern_input_line )
  
  p_ROC_FDR_Power_shared<-(p_ROC_shared + p_FDR_Power_shared) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  p_ROC_FDR_Power_shared<-p_ROC_FDR_Power_shared+ plot_layout(heights = c(2, 1))
  ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_Shared_Compare.pdf"),p_ROC_FDR_Power_shared,height=180, width=210, units = "mm",dpi=500)
  
  ###Plot for WB specific
  p_FDR_Power_WB <- FDR_Power_Compare_fun(
    FDR_Power_ancestry[grep("WB", FDR_Power_ancestry$Method), ] %>%
      mutate(Method = case_when(
        str_detect(Method, "^MESuSiE(\\b|_)")                          ~ "MESuSiE",
        str_detect(Method, "^MultiSuSiE(\\b|_)")                       ~ "MultiSuSiE",
        str_detect(Method, "^SuSiEx(\\b|_)")                           ~ "SuSiEx",
        str_detect(Method, "^SuSiE(_meta_weighted|_weighted)(\\b|$)")  ~ "SuSiE_weighted",
        str_detect(Method, "^SuSiE(_meta_merged|_merged)(\\b|$)")      ~ "SuSiE_merged",
        str_detect(Method, "^SuSiE(\\b|_)")                            ~ "SuSiE",
        str_detect(Method, "^Paintor(\\b|_)")                          ~ "Paintor",
        str_detect(Method, "^XMAP(\\b|_)")                             ~ "XMAP",
        str_detect(Method, "^CARMAX(\\b|_)")                           ~ "CARMAX",
        TRUE ~ Method
      )),
    pattern_input_Power
  )
  
  p_ROC_WB <- ROC_ancestry_Compare_fun(
    ROC_ancestry[grep("WB", ROC_ancestry$Method), ] %>%
      mutate(Method = case_when(
        str_detect(Method, "^MESuSiE(\\b|_)")                          ~ "MESuSiE",
        str_detect(Method, "^MultiSuSiE(\\b|_)")                       ~ "MultiSuSiE",
        str_detect(Method, "^SuSiEx(\\b|_)")                           ~ "SuSiEx",
        str_detect(Method, "^SuSiE(_meta_weighted|_weighted)(\\b|$)")  ~ "SuSiE_weighted",
        str_detect(Method, "^SuSiE(_meta_merged|_merged)(\\b|$)")      ~ "SuSiE_merged",
        str_detect(Method, "^SuSiE(\\b|_)")                            ~ "SuSiE",
        str_detect(Method, "^Paintor(\\b|_)")                          ~ "Paintor",
        str_detect(Method, "^XMAP(\\b|_)")                             ~ "XMAP",
        str_detect(Method, "^CARMAX(\\b|_)")                           ~ "CARMAX",
        TRUE ~ Method
      )),
    pattern_input_line
  )
  
  p_ROC_FDR_Power_WB<-(p_ROC_WB + p_FDR_Power_WB) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  p_ROC_FDR_Power_WB<-p_ROC_FDR_Power_WB+ plot_layout(heights = c(2, 1))
  ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_WB_Compare.pdf"),p_ROC_FDR_Power_WB,height=180, width=210, units = "mm",dpi=500)
  
  
  ###Plot for BB specific
  p_FDR_Power_BB <- FDR_Power_Compare_fun(
    FDR_Power_ancestry[grep("BB", FDR_Power_ancestry$Method), ] %>%
      mutate(Method = case_when(
        str_detect(Method, "^MESuSiE(\\b|_)")                          ~ "MESuSiE",
        str_detect(Method, "^MultiSuSiE(\\b|_)")                       ~ "MultiSuSiE",
        str_detect(Method, "^SuSiEx(\\b|_)")                           ~ "SuSiEx",
        str_detect(Method, "^SuSiE(_meta_weighted|_weighted)(\\b|$)")  ~ "SuSiE_weighted",
        str_detect(Method, "^SuSiE(_meta_merged|_merged)(\\b|$)")      ~ "SuSiE_merged",
        str_detect(Method, "^SuSiE(\\b|_)")                            ~ "SuSiE",
        str_detect(Method, "^Paintor(\\b|_)")                          ~ "Paintor",
        str_detect(Method, "^XMAP(\\b|_)")                             ~ "XMAP",
        str_detect(Method, "^CARMAX(\\b|_)")                           ~ "CARMAX",
        TRUE ~ Method
      )),
    pattern_input_Power
  )
  
  p_ROC_BB <- ROC_ancestry_Compare_fun(
    ROC_ancestry[grep("BB", ROC_ancestry$Method), ] %>%
      mutate(Method = case_when(
        str_detect(Method, "^MESuSiE(\\b|_)")                          ~ "MESuSiE",
        str_detect(Method, "^MultiSuSiE(\\b|_)")                       ~ "MultiSuSiE",
        str_detect(Method, "^SuSiEx(\\b|_)")                           ~ "SuSiEx",
        str_detect(Method, "^SuSiE(_meta_weighted|_weighted)(\\b|$)")  ~ "SuSiE_weighted",
        str_detect(Method, "^SuSiE(_meta_merged|_merged)(\\b|$)")      ~ "SuSiE_merged",
        str_detect(Method, "^SuSiE(\\b|_)")                            ~ "SuSiE",
        str_detect(Method, "^Paintor(\\b|_)")                          ~ "Paintor",
        str_detect(Method, "^XMAP(\\b|_)")                             ~ "XMAP",
        str_detect(Method, "^CARMAX(\\b|_)")                           ~ "CARMAX",
        TRUE ~ Method
      )),
    pattern_input_line
  )
  
  p_ROC_FDR_Power_BB<-(p_ROC_BB + p_FDR_Power_BB) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  p_ROC_FDR_Power_BB<-p_ROC_FDR_Power_BB+ plot_layout(heights = c(2, 1))
  ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_BB_Compare.pdf"),p_ROC_FDR_Power_BB,height=180, width=210, units = "mm",dpi=500)
  
}



#Compare_Plot(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/res_summary/"), "shared_50_external_ld_updated_xmap_meta_susie_rss.RData", "SuSiE RSS", "No Control")
#Compare_Plot(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/res_summary/"), "shared_50_external_ld_updated.RData", "External LD", "In sample LD")
