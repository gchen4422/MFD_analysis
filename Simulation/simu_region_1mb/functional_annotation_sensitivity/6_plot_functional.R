#res_out_name<-c("shared_50_baseline.RData","1000G_baseline.RData","BB_50000_baseline.RData","Low_cor_baseline.RData","shared_0_baseline.RData","shared_50_External_baseline.RData","shared_100_alternative.RData","Meta_N_baseline.RData")
res_out_name<-c("functional_annotation_sensitivity_w2.RData")
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
res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
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
                        "Paintor","Paintor_fun","Paintor_all_fun", "MultiSuSiE", "SuSiEx", "XMAP", "CARMA-X")
method_levels_base_2 <- c("MESuSiE", "SuSiE", "SuSiE_weighted","SuSiE_merged",
                        "Paintor","Paintor_fun","Paintor_all_fun","MultiSuSiE", "SuSiEx", "XMAP", "CARMA-X")


load(paste0(res_dir,res_out))
# 1. Clean Data & Rename Methods
all_Set_data_dataframe <- all_Set_data_dataframe %>%
  mutate(Method = recode(Method,
                         "SuSiE_meta_weighted" = "SuSiE-weighted",
                         "SuSiE_meta_merged"   = "SuSiE-merged",
                         "CARMAX"   = "CARMA-X"))

set_power_summary <- set_power_summary %>%
  mutate(Method = recode(Method,
                         "SuSiE_meta_weighted" = "SuSiE-weighted",
                         "SuSiE_meta_merged"   = "SuSiE-merged",
                         "CARMAX"   = "CARMA-X"))%>% mutate(Method  = factor(Method, levels = method_levels_base))

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
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    fill  = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal"
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
either_all_ROC_data_dataframe<-either_all_ROC_data_dataframe%>%mutate(Method = recode(Method, "CARMAX" = "CARMA-X"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","Paintor_fun","Paintor_all_fun","MultiSuSiE","SuSiEx","XMAP","CARMA-X"))
p_ROC_Either<-ROC_shared_fun(either_all_ROC_data_dataframe)
p_ROC_Either <- p_ROC_Either + guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
###################
#
#FDR&Power
#
###################
FDR_Power_either <- FDR_Power_either %>% mutate(Method = recode(Method, "CARMAX" = "CARMA-X"))
FDR_Power_either$Method <- factor(
  FDR_Power_either$Method,
  levels = method_levels_base_2
)

FDR_Power_either %>%
  filter(FDR == 0.05) %>%
  select(Method, h2, causal_num, Power) %>%
  bind_rows(
    FDR_Power_either %>%
      filter(FDR == 0.05) %>%
      group_by(Method) %>%
      summarise(Power = mean(Power)) %>%
      mutate(h2 = "Average", causal_num = "All")
  ) %>%
  print(n = Inf)

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

PIP_calibration_either_byh2<-create_obs_frq_byh2(data_all%>%select(Signal,h2,causal_num,MESuSiE_Either,SuSiE_Either, SuSiE_weighted_PIP , SuSiE_merged_PIP , Paintor_PIP, Paintor_fun_PIP,Paintor_all_fun_PIP,SuSiEx_PIP, XMAP_PIP,MultiSuSiE_PIP,CARMAX_Either),c(1,2,3),c("MESuSiE_Either","SuSiE_Either","SuSiE_weighted_PIP","SuSiE_merged_PIP","Paintor_PIP","Paintor_fun_PIP","Paintor_all_fun_PIP","SuSiEx_PIP","XMAP_PIP","MultiSuSiE_PIP","CARMAX_Either"))
PIP_calibration_either_byh2<- PIP_calibration_either_byh2%>%mutate(Method = fct_recode(Method, "MESuSiE" = "MESuSiE_Either","SuSiE" = "SuSiE_Either","SuSiE_weighted" = "SuSiE_weighted_PIP","SuSiE_merged" = "SuSiE_merged_PIP","Paintor" = "Paintor_PIP","Paintor_fun" = "Paintor_fun_PIP","Paintor_all_fun" = "Paintor_all_fun_PIP", "SuSiEx" = "SuSiEx_PIP","XMAP" = "XMAP_PIP","MultiSuSiE" = "MultiSuSiE_PIP","CARMA-X" = "CARMAX_Either"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","Paintor_fun","Paintor_all_fun","MultiSuSiE","SuSiEx","XMAP","CARMA-X"))

p_calibration_byh2<-PIP_calibration_shared_byh2_fun(PIP_calibration_either_byh2)

ROC_FDR_Power_Calibration_Either_Plot<-ggpubr::ggarrange(p_ROC_Either,p_FDR_Power_Either,p_calibration_byh2,nrow = 3,ncol=1,
                                                         common.legend = TRUE, legend="bottom",labels = c("a","b","c","d","e","f"),font.label=list(color="black",size=7))
ggsave(paste0(plot_dir,"ROC_FDR_Power_Calibration_Either.pdf"),ROC_FDR_Power_Calibration_Either_Plot,
       height = 180, width = 240,units = "mm",dpi=600)

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
shared_all_ROC_data_dataframe<-shared_all_ROC_data_dataframe%>%mutate(Method = recode(Method, "CARMAX" = "CARMA-X"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","Paintor_fun","Paintor_all_fun","MultiSuSiE","SuSiEx","XMAP","CARMA-X"))
p_ROC_shared<-ROC_shared_fun(shared_all_ROC_data_dataframe)
p_ROC_shared <- p_ROC_shared + guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
###################
#
#FDR&Power
#
###################
FDR_Power_shared <- FDR_Power_shared %>% mutate(Method = recode(Method, "CARMAX" = "CARMA-X"))
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


PIP_calibration_shared_byh2<-create_obs_frq_byh2(data_all%>%select(Signal,h2,causal_num,MESuSiE_Shared,SuSiE_Shared, SuSiE_weighted_PIP , SuSiE_merged_PIP , Paintor_PIP, Paintor_fun_PIP, Paintor_all_fun_PIP, SuSiEx_PIP, XMAP_PIP,MultiSuSiE_PIP,CARMAX_Shared),c(3),c("MESuSiE_Shared","SuSiE_Shared","SuSiE_weighted_PIP","SuSiE_merged_PIP","Paintor_PIP","Paintor_fun_PIP","Paintor_all_fun_PIP","SuSiEx_PIP","XMAP_PIP","MultiSuSiE_PIP","CARMAX_Shared"))
PIP_calibration_shared_byh2<- PIP_calibration_shared_byh2%>%mutate(Method = fct_recode(Method, "MESuSiE" = "MESuSiE_Shared","SuSiE" = "SuSiE_Shared","SuSiE_weighted" = "SuSiE_weighted_PIP","SuSiE_merged" = "SuSiE_merged_PIP","Paintor" = "Paintor_PIP","Paintor_fun" = "Paintor_fun_PIP","Paintor_all_fun" = "Paintor_all_fun_PIP", "SuSiEx" = "SuSiEx_PIP","XMAP" = "XMAP_PIP","MultiSuSiE" = "MultiSuSiE_PIP","CARMA-X" = "CARMAX_Shared"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","Paintor_fun","Paintor_all_fun","MultiSuSiE", "SuSiEx","XMAP","CARMA-X"))


p_calibration_shared_byh2<-PIP_calibration_shared_byh2_fun(PIP_calibration_shared_byh2)

ROC_FDR_Power_Calibration_shared_Plot<-ggpubr::ggarrange(p_ROC_shared,p_FDR_Power_shared,p_calibration_shared_byh2,nrow = 3,ncol=1,
                                                         common.legend = TRUE, legend="bottom",labels = c("a","b","c","d","e","f"),font.label=list(color="black",size=7))
ggsave(paste0(plot_dir,"ROC_FDR_Power_Calibration_shared.pdf"),ROC_FDR_Power_Calibration_shared_Plot,
       height = 180, width = 240,units = "mm",dpi=600)


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
  mutate(Method = as.character(Method)) %>%
  mutate(Method = recode(Method, "CARMAX WB" = "CARMA-X WB", "CARMAX BB" = "CARMA-X BB"))

split_list <- strsplit(ancestry_all_ROC_data_dataframe %>% pull(Method), " +")

ancestry_all_ROC_data_dataframe <- ancestry_all_ROC_data_dataframe %>%mutate(
  Method = sapply(split_list, `[`, 1),
  Ancestry = sapply(split_list, `[`, 2)
)%>%mutate(Method = fct_relevel(Method,c(
  "MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged",
  "Paintor","Paintor_fun","Paintor_all_fun","MultiSuSiE","SuSiEx","XMAP","CARMA-X"
)),Ancestry = fct_relevel(Ancestry, "WB","BB"))%>%mutate(Ancestry = fct_recode(Ancestry, "White British" = "WB"   , "Black British" = "BB" ))


p_ROC_ancestry<-ROC_ancestry_fun(ancestry_all_ROC_data_dataframe)

###################
#
#FDR&Power
#
###################
FDR_Power_ancestry <- FDR_Power_ancestry %>% mutate(Method = as.character(Method)) %>%
  mutate(Method = recode(Method, "CARMAX WB" = "CARMA-X WB", "CARMAX BB" = "CARMA-X BB"))
split_list <- strsplit(FDR_Power_ancestry %>% pull(Method), " +")
FDR_Power_ancestry <- FDR_Power_ancestry%>%ungroup(h2,causal_num,Method) %>%mutate(
  Method = sapply(split_list, `[`, 1),
  Ancestry = sapply(split_list, `[`, 2)
)%>%mutate(Method = fct_relevel(Method, "MESuSiE","SuSiE","SuSiE_weighted","SuSiE_merged","Paintor","Paintor_fun","Paintor_all_fun","MultiSuSiE","SuSiEx","XMAP","CARMA-X"),Ancestry = fct_relevel(Ancestry, "WB","BB"))%>%mutate(Ancestry = fct_recode(Ancestry, "White British" = "WB"   , "Black British" = "BB" ))

FDR_Power_ancestry<-FDR_Power_ancestry%>%filter(FDR!=0.5)
power_upper_limit<-FDR_Power_ancestry%>%summarise(upper_limit = min(ceiling(max(Power)*10)/10+0.1,1))%>%pull(upper_limit)

dyn.load("/apps/spack/negishi/apps/curl/7.85.0-gcc-12.2.0-4phaxqw/lib/libcurl.so.4",
         local = FALSE)

dyn.load("/apps/spack/negishi/apps/r-sf/1.0-9-gcc-12.2.0-sftozdw/rlib/R/library/sf/libs/sf.so")

p_FDR_Power_ancestry<-FDR_Power_ancestry_fun(FDR_Power_ancestry)+ylim(0, power_upper_limit)
ROC_FDR_Power_ancestry<- (p_ROC_ancestry / p_FDR_Power_ancestry) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
ROC_FDR_Power_ancestry<-ROC_FDR_Power_ancestry+ plot_layout(heights = c(1, 1))
ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_Ancestry.pdf"),ROC_FDR_Power_ancestry,height=180, width=240, units = "mm",dpi=600)


###################
#
#PIP calibration
#
###################

PIP_calibration_ancestry<- PIP_calibration_ancestry%>%group_by(causal_num)%>%mutate(Method = fct_recode(Method, "MESuSiE White British" = "MESuSiE~WB", "MESuSiE Black British" = "MESuSiE~BB","SuSiE White British" = "SuSiE~WB","SuSiE Black British" = "SuSiE~BB","SuSiE_weighted White British" = "SuSiE_weighted~WB","SuSiE_weighted Black British" = "SuSiE_weighted~BB","SuSiE_merged White British" = "SuSiE_merged~WB","SuSiE_merged Black British" = "SuSiE_merged~BB","Paintor White British" = "Paintor~WB","Paintor Black British" = "Paintor~BB","Paintor_fun White British" = "Paintor_fun~WB","Paintor_fun Black British" = "Paintor_fun~BB","Paintor_all_fun White British" = "Paintor_all_fun~WB","Paintor_all_fun Black British" = "Paintor_all_fun~BB","SuSiEx White British" = "SuSiEx~WB","SuSiEx Black British" = "SuSiEx~BB","XMAP White British" = "XMAP~WB","XMAP Black British" = "XMAP~BB","MultiSuSiE White British" = "MultiSuSiE~WB","MultiSuSiE Black British" = "MultiSuSiE~BB","CARMA-X White British" = "CARMAX~WB","CARMA-X Black British" = "CARMAX~BB"))
levels(PIP_calibration_ancestry$Method)<-c(paste0("MESuSiE~","White~","British"),paste0("MESuSiE~","Black~","British"),paste0("SuSiE~","White~","British"),paste0("SuSiE~","Black~","British"),paste0("SuSiE_weighted~","White~","British"),paste0("SuSiE_weighted~","Black~","British"),paste0("SuSiE_merged~","White~","British"),paste0("SuSiE_merged~","Black~","British"),paste0("Paintor~","White~","British"),paste0("Paintor~","Black~","British"),paste0("Paintor_fun~","White~","British"),paste0("Paintor_fun~","Black~","British"),paste0("Paintor_all_fun~","White~","British"),paste0("Paintor_all_fun~","Black~","British"),paste0("MultiSuSiE~","White~","British"),paste0("MultiSuSiE~","Black~","British"),paste0("SuSiEx~","White~","British"),paste0("SuSiEx~","Black~","British"),paste0("XMAP~","White~","British"),paste0("XMAP~","Black~","British"),paste0("CARMA-X~","White~","British"),paste0("CARMA-X~","Black~","British"))
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
ggsave(paste0(plot_dir,"PIP_Calibration_WB.pdf"),p_calibration_ancestry_WB,height = 300, width = 210,units = "mm",dpi=500)
ggsave(paste0(plot_dir,"PIP_Calibration_BB.pdf"),p_calibration_ancestry_BB,height = 300, width = 210,units = "mm",dpi=500)

#ggsave(paste0(plot_dir,"PIP_Calibration_Ancestry.pdf"),p_out,height = 180, width = 210,units = "mm",dpi=500)
#}


###Time Plot
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(openxlsx)
library(stringr)

source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/utility.R")
res_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/"
simulation_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/functional_annotation_sensitivity/"

# ============================================================
# 1. Paintor (baseline) from time_all
# ============================================================
# ============================================================
# 1. Paintor (baseline) time + memory from shared_50 Runtime_memory files
# ============================================================
baseline_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/"

paintor_base_list <- list()

for(causal_num in 1:3) {
  for(h2 in 1:2) {
    wrk_dir <- paste0(baseline_dir, "causal_num_", causal_num, "/")
    result_dir <- paste0(wrk_dir, "result/")
    data_dir <- paste0(wrk_dir, "summary_data/")

    for(LOCI_num in 1:100) {
      f_main <- paste0(result_dir, "Runtime_memory_CAUSAL_", causal_num, "_LOCI_", LOCI_num, "_h2_", h2, ".txt")
      f_zfile <- paste0(data_dir, "CAUSAL_", causal_num, "_LOCI_", LOCI_num, "_h2_", h2)

      if(!file.exists(f_main) || !file.exists(f_zfile)) next

      d <- fread(f_main) %>% filter(Method == "Paintor")
      nsnp <- as.integer(system(paste0("wc -l < ", f_zfile), intern = TRUE)) - 1
      d$NUM_SNP <- nsnp

      paintor_base_list[[length(paintor_base_list) + 1]] <- d
    }
  }
  cat(paste("Done baseline Causal:", causal_num, "\n"))
}

paintor_base_all <- rbindlist(paintor_base_list)

summary_paintor_base <- paintor_base_all %>%
  group_by(NUM_SNP, Method) %>%
  summarise(
    Mean_Time = mean(Time_Min, na.rm = TRUE),
    Mean_RAM  = mean(Peak_RAM_Used_MiB, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 2. Paintor_fun & Paintor_all_fun from Runtime_memory files
# ============================================================
paintor_fun_list <- list()

for(causal_num in 1:3) {
  for(h2 in 1:2) {
    wrk_dir <- paste0(simulation_dir, "causal_num_", causal_num, "/")
    result_dir <- paste0(wrk_dir, "result/")
    data_dir <- paste0(wrk_dir, "summary_data/w2_startover/")

    for(LOCI_num in 1:100) {
      f_mem <- paste0(result_dir, "Runtime_memory_CAUSAL_", causal_num, "_LOCI_", LOCI_num, "_h2_", h2, ".txt")
      f_zfile <- paste0(data_dir, "CAUSAL_", causal_num, "_LOCI_", LOCI_num, "_h2_", h2)

      if(!file.exists(f_mem) || !file.exists(f_zfile)) next

      d <- fread(f_mem)
      nsnp <- as.integer(system(paste0("wc -l < ", f_zfile), intern = TRUE)) - 1
      d$NUM_SNP <- nsnp

      paintor_fun_list[[length(paintor_fun_list) + 1]] <- d
    }
  }
  cat(paste("Done Causal:", causal_num, "\n"))
}

paintor_fun_all <- rbindlist(paintor_fun_list)

summary_paintor_fun <- paintor_fun_all %>%
  group_by(NUM_SNP, Method) %>%
  summarise(
    Mean_Time = mean(Time_Min, na.rm = TRUE),
    Mean_RAM  = mean(Peak_RAM_Used_MiB, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 3. Combine Paintor methods only
# ============================================================
summary_perf <- bind_rows(summary_paintor_base, summary_paintor_fun)
summary_perf$Method <- factor(summary_perf$Method, levels = c("Paintor", "Paintor_fun", "Paintor_all_fun"))

# ============================================================
# 4. Define Colors & Plot
# ============================================================
my_colors <- c("Paintor"="#fc8d62", "Paintor_fun"="#e5703e", "Paintor_all_fun"="#d4533a")

p_time <- ggplot(summary_perf, aes(x = NUM_SNP, y = Mean_Time, color = Method)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "loess", se = FALSE, span = 0.8) +
  scale_y_log10() +
  scale_color_manual(values = my_colors) +
  labs(x = "Number of SNPs", y = "Time (minutes, log scale)") +
  theme_bw() + custom_theme()

p_mem <- ggplot(summary_perf %>% filter(!is.na(Mean_RAM)),
                aes(x = NUM_SNP, y = Mean_RAM, color = Method)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "loess", se = FALSE, span = 0.8) +
  scale_y_log10() +
  scale_color_manual(values = my_colors) +
  labs(x = "Number of SNPs", y = "Peak Memory (MiB, log scale)") +
  theme_bw() + custom_theme()

p_perf_out <- (p_time / p_mem) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a') &
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold"))

print(p_perf_out)
ggsave(paste0(res_dir, "Figure/shared_50_baseline_updated_meta_updated_xmap/Time_paintor.pdf"),
       p_perf_out, height = 180, width = 180, units = "mm", dpi = 600)

# ============================================================
# 5. Summary table
# ============================================================
perf_method_summary <- summary_perf %>%
  group_by(Method) %>%
  summarise(
    n_regions    = n(),
    time_mean    = mean(Mean_Time, na.rm = TRUE),
    time_median  = median(Mean_Time, na.rm = TRUE),
    time_q25     = quantile(Mean_Time, 0.25, na.rm = TRUE),
    time_q75     = quantile(Mean_Time, 0.75, na.rm = TRUE),
    time_p95     = quantile(Mean_Time, 0.95, na.rm = TRUE),
    ram_mean     = mean(Mean_RAM, na.rm = TRUE),
    ram_median   = median(Mean_RAM, na.rm = TRUE),
    ram_q25      = quantile(Mean_RAM, 0.25, na.rm = TRUE),
    ram_q75      = quantile(Mean_RAM, 0.75, na.rm = TRUE),
    ram_p95      = quantile(Mean_RAM, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

res_dir_table <- paste0(res_dir, "supp_tables/")
system(paste0("mkdir -p ", res_dir_table))
full_path <- paste0(res_dir_table, "Time_memory_paintor.xlsx")

wb <- createWorkbook()
addWorksheet(wb, "Time_memory")
title_style <- createStyle(textDecoration = "bold", fontSize = 12)
writeData(wb, "Time_memory", "Supplementary Table: Runtime and Memory Summary (Paintor variants)", startRow = 1)
addStyle(wb, "Time_memory", title_style, rows = 1, cols = 1)
writeData(wb, "Time_memory", perf_method_summary, startRow = 2)
saveWorkbook(wb, full_path, overwrite = TRUE)
print(paste("Saved combined table to:", full_path))
