#res_out_name<-c("shared_50_baseline.RData","1000G_baseline.RData","BB_50000_baseline.RData","Low_cor_baseline.RData","shared_0_baseline.RData","shared_50_External_baseline.RData","shared_100_alternative.RData","Meta_N_baseline.RData")
res_out_name<-"admixed_population_sensitivity.RData"
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

# Use scenario-specific facet labels
# utility.R assumes baseline params (h2=1e-4/2e-4, causal=1/3/5)
# admix uses h2_locus=0.005 and num_causal_SNP=c(2,4,6)
relabel_admix <- function(df) {
  if ("h2" %in% names(df) && is.factor(df$h2)) {
    lvls <- levels(df$h2)
    lvls[lvls == "~h^2 == 10^-4"] <- "~h^2 == 5%*%10^-3"
    levels(df$h2) <- lvls
  }
  if ("causal_num" %in% names(df) && is.factor(df$causal_num)) {
    lvls <- levels(df$causal_num)
    lvls <- sub("== 1 ", "== 2 ", lvls)
    lvls <- sub("== 3 ", "== 4 ", lvls)
    lvls <- sub("== 5 ", "== 6 ", lvls)
    levels(df$causal_num) <- lvls
  }
  df
}

res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_admix/res_summary/")
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
                        "Paintor", "MultiSuSiE", "SuSiEx", "XMAP", "CARMA-X")
method_levels_base_2 <- c("MESuSiE", "SuSiE", "SuSiE_weighted","SuSiE_merged",
                        "Paintor","MultiSuSiE", "SuSiEx", "XMAP", "CARMA-X")


load(paste0(res_dir,res_out))

# Relabel loaded data to admix-specific h2 and causal_num values
all_Set_data_dataframe <- relabel_admix(all_Set_data_dataframe)
set_power_summary <- relabel_admix(set_power_summary)
either_all_ROC_data_dataframe <- relabel_admix(either_all_ROC_data_dataframe)
shared_all_ROC_data_dataframe <- relabel_admix(shared_all_ROC_data_dataframe)
ancestry_all_ROC_data_dataframe <- relabel_admix(ancestry_all_ROC_data_dataframe)
FDR_Power_either <- relabel_admix(FDR_Power_either)
FDR_Power_shared <- relabel_admix(FDR_Power_shared)
FDR_Power_ancestry <- relabel_admix(FDR_Power_ancestry)
PIP_calibration_ancestry <- relabel_admix(PIP_calibration_ancestry)

# 1. Clean Data & Rename Methods (including CARMAX -> CARMA-X display label)
rename_carmax <- function(df) {
  if ("Method" %in% names(df)) {
    df$Method <- gsub("^CARMAX$", "CARMA-X", as.character(df$Method))
    df$Method <- gsub("^CARMAX ", "CARMA-X ", df$Method)
  }
  df
}
all_Set_data_dataframe <- rename_carmax(all_Set_data_dataframe)
set_power_summary <- rename_carmax(set_power_summary)
either_all_ROC_data_dataframe <- rename_carmax(either_all_ROC_data_dataframe)
shared_all_ROC_data_dataframe <- rename_carmax(shared_all_ROC_data_dataframe)
ancestry_all_ROC_data_dataframe <- rename_carmax(ancestry_all_ROC_data_dataframe)
FDR_Power_either <- rename_carmax(FDR_Power_either)
FDR_Power_shared <- rename_carmax(FDR_Power_shared)
FDR_Power_ancestry <- rename_carmax(FDR_Power_ancestry)
PIP_calibration_ancestry <- rename_carmax(PIP_calibration_ancestry)

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
    filter(Method == "Paintor", causal_num == "Num~Causal  == 6 ", h2 == "~h^2 == 5%*%10^-3") %>%
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
either_all_ROC_data_dataframe<-either_all_ROC_data_dataframe%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","MultiSuSiE","SuSiEx","XMAP","CARMA-X"))
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

PIP_calibration_either_byh2<-create_obs_frq_byh2(data_all%>%select(Signal,h2,causal_num,MESuSiE_Either,SuSiE_Either, SuSiE_weighted_PIP , SuSiE_merged_PIP , Paintor_PIP,SuSiEx_PIP, XMAP_PIP,MultiSuSiE_PIP,CARMAX_Either),c(1,2,3),c("MESuSiE_Either","SuSiE_Either","SuSiE_weighted_PIP","SuSiE_merged_PIP","Paintor_PIP","SuSiEx_PIP","XMAP_PIP","MultiSuSiE_PIP","CARMAX_Either"))
PIP_calibration_either_byh2<- PIP_calibration_either_byh2%>%mutate(Method = fct_recode(Method, "MESuSiE" = "MESuSiE_Either","SuSiE" = "SuSiE_Either","SuSiE_weighted" = "SuSiE_weighted_PIP","SuSiE_merged" = "SuSiE_merged_PIP","Paintor" = "Paintor_PIP", "SuSiEx" = "SuSiEx_PIP","XMAP" = "XMAP_PIP","MultiSuSiE" = "MultiSuSiE_PIP","CARMA-X" = "CARMAX_Either"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","MultiSuSiE","SuSiEx","XMAP","CARMA-X"))
PIP_calibration_either_byh2 <- relabel_admix(PIP_calibration_either_byh2)

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
shared_all_ROC_data_dataframe<-shared_all_ROC_data_dataframe%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","MultiSuSiE","SuSiEx","XMAP","CARMA-X"))
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
PIP_calibration_shared_byh2<- PIP_calibration_shared_byh2%>%mutate(Method = fct_recode(Method, "MESuSiE" = "MESuSiE_Shared","SuSiE" = "SuSiE_Shared","SuSiE_weighted" = "SuSiE_weighted_PIP","SuSiE_merged" = "SuSiE_merged_PIP","Paintor" = "Paintor_PIP", "SuSiEx" = "SuSiEx_PIP","XMAP" = "XMAP_PIP","MultiSuSiE" = "MultiSuSiE_PIP","CARMA-X" = "CARMAX_Shared"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","MultiSuSiE", "SuSiEx","XMAP","CARMA-X"))
PIP_calibration_shared_byh2 <- relabel_admix(PIP_calibration_shared_byh2)

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
  "Paintor","MultiSuSiE","SuSiEx","XMAP","CARMA-X"
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
)%>%mutate(Method = fct_relevel(Method, "MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","MultiSuSiE","SuSiEx","XMAP","CARMA-X"),Ancestry = fct_relevel(Ancestry, "WB","BB"))%>%mutate(Ancestry = fct_recode(Ancestry, "White British" = "WB"   , "Black British" = "BB" ))

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

PIP_calibration_ancestry<- PIP_calibration_ancestry%>%group_by(causal_num)%>%mutate(Method = fct_recode(Method, "MESuSiE White British" = "MESuSiE~WB", "MESuSiE Black British" = "MESuSiE~BB","SuSiE White British" = "SuSiE~WB","SuSiE Black British" = "SuSiE~BB","SuSiE_weighted White British" = "SuSiE_weighted~WB","SuSiE_weighted Black British" = "SuSiE_weighted~BB","SuSiE_merged White British" = "SuSiE_merged~WB","SuSiE_merged Black British" = "SuSiE_merged~BB","Paintor White British" = "Paintor~WB","Paintor Black British" = "Paintor~BB","SuSiEx White British" = "SuSiEx~WB","SuSiEx Black British" = "SuSiEx~BB","XMAP White British" = "XMAP~WB","XMAP Black British" = "XMAP~BB","MultiSuSiE White British" = "MultiSuSiE~WB","MultiSuSiE Black British" = "MultiSuSiE~BB","CARMA-X White British" = "CARMA-X~WB","CARMA-X Black British" = "CARMA-X~BB"))
levels(PIP_calibration_ancestry$Method)<-c(paste0("MESuSiE~","White~","British"),paste0("MESuSiE~","Black~","British"),paste0("SuSiE~","White~","British"),paste0("SuSiE~","Black~","British"),paste0("SuSiE_weighted~","White~","British"),paste0("SuSiE_weighted~","Black~","British"),paste0("SuSiE_merged~","White~","British"),paste0("SuSiE_merged~","Black~","British"),paste0("Paintor~","White~","British"),paste0("Paintor~","Black~","British"),paste0("MultiSuSiE~","White~","British"),paste0("MultiSuSiE~","Black~","British"),paste0("SuSiEx~","White~","British"),paste0("SuSiEx~","Black~","British"),paste0("XMAP~","White~","British"),paste0("XMAP~","Black~","British"),paste0("CARMA-X~","White~","British"),paste0("CARMA-X~","Black~","British"))
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
ggsave(paste0(plot_dir,"PIP_Calibration_WB.pdf"),p_calibration_ancestry_WB,height = 250, width = 210,units = "mm",dpi=500)
ggsave(paste0(plot_dir,"PIP_Calibration_BB.pdf"),p_calibration_ancestry_BB,height = 250, width = 210,units = "mm",dpi=500)

#ggsave(paste0(plot_dir,"PIP_Calibration_Ancestry.pdf"),p_out,height = 180, width = 210,units = "mm",dpi=500)
#}
