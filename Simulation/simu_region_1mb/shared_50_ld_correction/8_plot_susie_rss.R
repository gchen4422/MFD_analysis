
# containers to collect data across all 12 files
all_Set_data_list        <- list()
set_power_list           <- list()
either_ROC_list          <- list()
FDR_Power_either_list    <- list()
shared_ROC_list          <- list()
FDR_Power_shared_list    <- list()
ancestry_ROC_list        <- list()
FDR_Power_ancestry_list  <- list()



res_out_name_susie_rss <- c("flipped_1_rho_1_susie_rss_updated.RData","flipped_1_rho_2_susie_rss_updated.RData","flipped_1_rho_3_susie_rss_updated.RData","flipped_1_rho_4_susie_rss_updated.RData","flipped_2_rho_1_susie_rss_updated.RData","flipped_2_rho_2_susie_rss_updated.RData","flipped_2_rho_3_susie_rss_updated.RData","flipped_2_rho_4_susie_rss_updated.RData","flipped_3_rho_1_susie_rss_updated.RData","flipped_3_rho_2_susie_rss_updated.RData","flipped_3_rho_3_susie_rss_updated.RData","flipped_3_rho_4_susie_rss_updated.RData")

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
source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/utility.R")
res_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
system(paste0("mkdir -p ",res_dir))
for(res_out in res_out_name_susie_rss){
  
  load(paste0(res_dir,res_out))   
  
  # parse flip / rho from file name, e.g. "flipped_1_rho_2.RData"
  base_name <- sub("\\.RData$", "", res_out)
  parts     <- strsplit(base_name, "_")[[1]]  # c("flipped","1","rho","2")
  flip_id   <- parts[2]
  rho_id    <- parts[4]
  
  # tag and store for later combined plots
  all_Set_data_list[[base_name]]       <- all_Set_data_dataframe %>%
    mutate(flip = flip_id, rho = rho_id)
  
  set_power_list[[base_name]]          <- set_power_summary %>%
    mutate(flip = flip_id, rho = rho_id)
  
  either_ROC_list[[base_name]]         <- either_all_ROC_data_dataframe %>%
    mutate(flip = flip_id, rho = rho_id)
  
  FDR_Power_either_list[[base_name]]   <- FDR_Power_either %>%
    mutate(flip = flip_id, rho = rho_id)
  
  shared_ROC_list[[base_name]]         <- shared_all_ROC_data_dataframe %>%
    mutate(flip = flip_id, rho = rho_id)
  
  FDR_Power_shared_list[[base_name]]   <- FDR_Power_shared %>%
    mutate(flip = flip_id, rho = rho_id)
  
  ancestry_ROC_list[[base_name]]       <- ancestry_all_ROC_data_dataframe %>%
    mutate(flip = flip_id, rho = rho_id)
  
  FDR_Power_ancestry_list[[base_name]] <- FDR_Power_ancestry %>%
    mutate(flip = flip_id, rho = rho_id)
  
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
  upper_limit<-round(all_Set_data_dataframe%>%filter(Method == "Paintor",causal_num=="Num~Causal  == 1 ",h2=="~h^2 == 10^-4")%>%summarise(upper = quantile(Size,0.75))%>%pull(upper))+50
  p_size_box<-Set_Size_fun(all_Set_data_dataframe%>%mutate(Size = log2(Size+1)),upper_limit = log2(upper_limit))
  p_size_box<-p_size_box+ ylab("log2(Set Size + 1)")
  
  p_power_bar<-Set_Power_fun(set_power_summary)
  
  #p_power_bar<-p_power_bar+ylim(0,min(max(set_power_summary$Power_name)+0.05,1))
  size_power<-p_size_box/p_power_bar+plot_annotation(tag_levels = 'a')& 
    theme(plot.tag = element_text(size = 7,face="bold"))
  
  ggsave(paste0(plot_dir,"Set_size_power.pdf"),size_power,height=120, width=150, units = "mm",dpi=500)
  
  
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
  either_all_ROC_data_dataframe<-either_all_ROC_data_dataframe%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX"))
  p_ROC_Either<-ROC_shared_fun(either_all_ROC_data_dataframe)
  ###################
  #
  #FDR&Power
  #
  ###################
  power_upper_limit<-FDR_Power_either%>%filter(FDR!=0.5)%>%ungroup(Method,h2,causal_num)%>%summarise(upper_limit = min(ceiling(max(Power)*10)/10+0.1,1))%>%pull(upper_limit)
  p_FDR_Power_Either<-FDR_Power_shared_fun(FDR_Power_either%>%filter(FDR!=0.5))+ylim(0,power_upper_limit)
  #ROC_FDR_Power_Either <- (p_ROC_Either / p_FDR_Power_Either) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  
  ####################
  #
  #PIP calibration
  #
  ####################
  # ROC_FDR_Power_Either <- (p_ROC_Either / p_FDR_Power_Either) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  #  ROC_FDR_Power_Either<-ROC_FDR_Power_Either+ plot_layout(guides = "collect",heights = c(2, 1))&theme(legend.position = 'bottom')
  
  PIP_calibration_either_byh2<-create_obs_frq_byh2(data_all%>%select(Signal,h2,MESuSiE_Either,SuSiE_Either, Paintor_PIP,SuSiEx_PIP,XMAP_PIP,MultiSuSiE_PIP),c(1,2,3),c("MESuSiE_Either","SuSiE_Either","Paintor_PIP","SuSiEx_PIP","XMAP_PIP","MultiSuSiE_PIP"))
  PIP_calibration_either_byh2<- PIP_calibration_either_byh2%>%mutate(Method = fct_recode(Method, "MESuSiE" = "MESuSiE_Either","SuSiE" = "SuSiE_Either","Paintor" = "Paintor_PIP", "SuSiEx" = "SuSiEx_PIP","XMAP" = "XMAP_PIP","MultiSuSiE" = "MultiSuSiE_PIP"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","Paintor", "SuSiEx","XMAP","MultiSuSiE"))
  
  p_calibration_byh2<-PIP_calibration_shared_byh2_fun(PIP_calibration_either_byh2)
  
  ROC_FDR_Power_Calibration_Either_Plot<-ggarrange(p_ROC_Either,p_FDR_Power_Either,p_calibration_byh2,nrow = 3,ncol=1,
                                                   common.legend = TRUE, legend="bottom",labels = c("a","b","c","d","e","f"),font.label=list(color="black",size=7))
  ggsave(paste0(plot_dir,"ROC_FDR_Power_Calibration_Either.pdf"),ROC_FDR_Power_Calibration_Either_Plot,
         height = 210, width = 300,units = "mm",dpi=500)
  
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
  shared_all_ROC_data_dataframe<-shared_all_ROC_data_dataframe%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX"))
  p_ROC_shared<-ROC_shared_fun(shared_all_ROC_data_dataframe)
  ###################
  #
  #FDR&Power
  #
  ###################
  
  power_upper_limit<-FDR_Power_shared%>%filter(FDR!=0.5)%>%ungroup(Method,h2,causal_num)%>%summarise(upper_limit = min(ceiling(max(Power)*10)/10+0.1,1))%>%pull(upper_limit)
  p_FDR_Power_shared<-FDR_Power_shared_fun(FDR_Power_shared%>%filter(FDR!=0.5))+ylim(0,power_upper_limit)
  
  #  p_FDR_Power_shared<-FDR_Power_shared_fun(FDR_Power_shared)+ylim(0, min(ceiling(max(FDR_Power_shared$Power)*10)/10+0.1,1))
  #  ROC_FDR_Power_shared <- (p_ROC_shared / p_FDR_Power_shared) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  #  ROC_FDR_Power_shared<-ROC_FDR_Power_shared+ plot_layout(guides = "collect",heights = c(2, 1))&theme(legend.position = 'bottom')
  #  ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_Shared.pdf"),ROC_FDR_Power_shared,height=180, width=180, units = "mm",dpi=500)
  
  ####################
  #
  #PIP calibration
  #
  ####################
  
  PIP_calibration_shared_byh2<-create_obs_frq_byh2(data_all%>%select(Signal,h2,causal_num,MESuSiE_Shared,SuSiE_Shared,Paintor_PIP, SuSiEx_PIP, XMAP_PIP, MultiSuSiE_PIP,CARMAX_Shared),c(3),c("MESuSiE_Shared","SuSiE_Shared","Paintor_PIP","SuSiEx_PIP","XMAP_PIP","MultiSuSiE_PIP","CARMAX_Shared"))
  PIP_calibration_shared_byh2<- PIP_calibration_shared_byh2%>%mutate(Method = fct_recode(Method, "MESuSiE" = "MESuSiE_Shared","SuSiE" = "SuSiE_Shared","Paintor" = "Paintor_PIP","SuSiEx" = "SuSiEx_PIP","XMAP" = "XMAP_PIP","MultiSuSiE" = "MultiSuSiE_PIP","CARMAX" = "CARMAX_Shared"))%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX"))
  
  p_calibration_shared_byh2<-PIP_calibration_shared_byh2_fun(PIP_calibration_shared_byh2)
  
  ROC_FDR_Power_Calibration_shared_Plot<-ggarrange(p_ROC_shared,p_FDR_Power_shared,p_calibration_shared_byh2,nrow = 3,ncol=1,
                                                   common.legend = TRUE, legend="bottom",labels = c("a","b","c","d","e","f"),font.label=list(color="black",size=7))
  ggsave(paste0(plot_dir,"ROC_FDR_Power_Calibration_shared.pdf"),ROC_FDR_Power_Calibration_shared_Plot,
         height = 210, width = 300,units = "mm",dpi=500)
  
  
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
    "MESuSiE","SuSiE","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX"
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
  )%>%mutate(Method = fct_relevel(Method,"MESuSiE","SuSiE","Paintor","SuSiEx","XMAP","MultiSuSiE","CARMAX"),Ancestry = fct_relevel(Ancestry, "WB","BB"))%>%mutate(Ancestry = fct_recode(Ancestry, "White British" = "WB"   , "Black British" = "BB" ))
  
  FDR_Power_ancestry<-FDR_Power_ancestry%>%filter(FDR!=0.5)
  power_upper_limit<-FDR_Power_ancestry%>%summarise(upper_limit = min(ceiling(max(Power)*10)/10+0.1,1))%>%pull(upper_limit)
  
  p_FDR_Power_ancestry<-FDR_Power_ancestry_fun(FDR_Power_ancestry)+ylim(0, power_upper_limit)
  ROC_FDR_Power_ancestry<- (p_ROC_ancestry / p_FDR_Power_ancestry) +plot_annotation(tag_levels = 'a')&theme(plot.tag = element_text(size = 7, face = "bold"))
  ROC_FDR_Power_ancestry<-ROC_FDR_Power_ancestry+ plot_layout(heights = c(1, 1))
  ggsave(paste0(plot_dir,"PIP_ROC_FDR_Power_Ancestry.pdf"),ROC_FDR_Power_ancestry,height=180, width=300, units = "mm",dpi=500)
  
  ###################
  #
  #PIP calibration
  #
  ###################   
  
  PIP_calibration_ancestry<- PIP_calibration_ancestry%>%group_by(causal_num)%>%mutate(Method = fct_recode(Method, "MESuSiE White British" = "MESuSiE~WB", "MESuSiE Black British" = "MESuSiE~BB","SuSiE White British" = "SuSiE~WB","SuSiE Black British" = "SuSiE~BB","Paintor White British" = "Paintor~WB","Paintor Black British" = "Paintor~BB","SuSiEx White British" = "SuSiEx~WB","SuSiEx Black British" = "SuSiEx~BB","XMAP White British" = "XMAP~WB","XMAP Black British" = "XMAP~BB","MultiSuSiE White British" = "MultiSuSiE~WB","MultiSuSiE Black British" = "MultiSuSiE~BB","CARMAX White British" = "CARMAX~WB","CARMAX Black British" = "CARMAX~BB"))
  levels(PIP_calibration_ancestry$Method)<-c(paste0("MESuSiE~","White~","British"),paste0("MESuSiE~","Black~","British"),paste0("SuSiE~","White~","British"),paste0("SuSiE~","Black~","British"),paste0("SuSiE_weighted~","White~","British"),paste0("SuSiE_weighted~","Black~","British"),paste0("SuSiE_merged~","White~","British"),paste0("SuSiE_merged~","Black~","British"),paste0("Paintor~","White~","British"),paste0("Paintor~","Black~","British"),paste0("SuSiEx~","White~","British"),paste0("SuSiEx~","Black~","British"),paste0("XMAP~","White~","British"),paste0("XMAP~","Black~","British"),paste0("MultiSuSiE~","White~","British"),paste0("MultiSuSiE~","Black~","British"),paste0("CARMAX~","White~","British"),paste0("CARMAX~","Black~","British"))
  p_calibration_ancestry<-PIP_calibration_ancestry_fun(PIP_calibration_ancestry)
  ggsave(paste0(plot_dir,"PIP_Calibration_Ancestry.pdf"),p_calibration_ancestry,height = 900, width = 242,units = "mm",dpi=500)
}      



library(dplyr)
library(patchwork)
library(ggpubr)

# bind all data
all_Set_data_all       <- bind_rows(all_Set_data_list)
set_power_all          <- bind_rows(set_power_list)
either_ROC_all         <- bind_rows(either_ROC_list)
FDR_Power_either_all   <- bind_rows(FDR_Power_either_list)
shared_ROC_all         <- bind_rows(shared_ROC_list)
FDR_Power_shared_all   <- bind_rows(FDR_Power_shared_list)
ancestry_ROC_all       <- bind_rows(ancestry_ROC_list)
FDR_Power_ancestry_all <- bind_rows(FDR_Power_ancestry_list)

# global upper limit for set size
upper_limit_all <- all_Set_data_all %>%
  filter(Method == "Paintor",
         causal_num == "Num~Causal  == 1 ",
         h2 == "~h^2 == 10^-4") %>%
  summarise(upper = quantile(Size, 0.75)) %>%
  pull(upper) %>%
  round() + 50

# make sure flip / rho are factors or characters
all_Set_data_all <- all_Set_data_all %>%
  mutate(
    flip = factor(flip),
    rho  = factor(rho)
  )

set_power_all <- set_power_all %>%
  mutate(
    flip = factor(flip),
    rho  = factor(rho)
  )

# mapping: flip(1,2,3) -> (0.001, 0.01, 0.1)
flip_map <- c("1" = "0.001",
              "2" = "0.01",
              "3" = "0.1")

# mapping: rho(1,2,3,4) -> (0.5, 0, -0.5, -1)
rho_map  <- c("1" = "0.5",
              "2" = "0",
              "3" = "-0.5",
              "4" = "-1")

facet_labeller <- labeller(
  flip = function(x) paste0("Flip prop. = ", flip_map[as.character(x)]),
  rho  = function(x) paste0("Rho sig. = ",  rho_map[as.character(x)])
)

p_size_box_all <- Set_Size_fun(
  all_Set_data_all %>% mutate(Size = log2(Size + 1)),
  upper_limit = log2(upper_limit_all)
) +
  ylab("log2(Set Size + 1)") +
  facet_grid(flip ~ rho, labeller = facet_labeller)

p_power_bar_all <- Set_Power_fun(set_power_all) +
  facet_grid(flip ~ rho, labeller = facet_labeller)


#size_power_all <- p_size_box_all / p_power_bar_all +
#  plot_annotation(tag_levels = "a") &
#  theme(plot.tag = element_text(size = 7, face = "bold"))

ggsave(
  file.path(res_dir, "Figure", "LD_susie_rss_Set_size_all_flips_rhos_combined.pdf"),
  p_size_box_all,
  height = 210, width = 300, units = "mm", dpi = 500
)


ggsave(
  file.path(res_dir, "Figure", "LD_susie_rss_power_all_flips_rhos_combined.pdf"),
  p_power_bar_all,
  height = 210, width = 300, units = "mm", dpi = 500
)

# -------- Either ancestry: ROC + FDR/Power (combined) --------
either_ROC_all <- either_ROC_all %>%
  mutate(Method = fct_relevel(Method,
                              "MESuSiE","SuSiE","Paintor",
                              "SuSiEx","XMAP","MultiSuSiE","CARMAX"))

p_ROC_Either_all <- ROC_shared_fun(either_ROC_all) +
  facet_grid(flip ~ rho, labeller = facet_labeller)

power_upper_limit_either_all <- FDR_Power_either_all %>%
  filter(FDR != 0.5) %>%
  ungroup() %>%
  summarise(upper_limit = min(ceiling(max(Power) * 10) / 10 + 0.1, 1)) %>%
  pull(upper_limit)

p_FDR_Power_Either_all <- FDR_Power_shared_fun(
  FDR_Power_either_all %>% filter(FDR != 0.5)
) +
  ylim(0, power_upper_limit_either_all) +
  facet_grid(flip ~ rho, labeller = facet_labeller)

#ROC_FDR_Either_all <- ggarrange(
#  p_ROC_Either_all, p_FDR_Power_Either_all,
#  nrow = 2, ncol = 1,
#  common.legend = TRUE, legend = "bottom",
#  labels = c("a","b"), font.label = list(color = "black", size = 7)
#)


ggsave(
  file.path(res_dir, "Figure", "LD_susie_rss_roc_either_all_flips_rhos_combined.pdf"),
  p_ROC_Either_all,
  height = 210, width = 300, units = "mm", dpi = 500
)


ggsave(
  file.path(res_dir, "Figure", "LD_susie_rss_fdr_power_either_all_flips_rhos_combined.pdf"),
  p_FDR_Power_Either_all,
  height = 210, width = 300, units = "mm", dpi = 500
)


# -------- Shared signal: ROC + FDR/Power (combined) --------
shared_ROC_all <- shared_ROC_all %>%
  mutate(Method = fct_relevel(Method,
                              "MESuSiE","SuSiE","Paintor",
                              "SuSiEx","XMAP","MultiSuSiE","CARMAX"))

p_ROC_shared_all <- ROC_shared_fun(shared_ROC_all) +
  facet_grid(flip ~ rho, labeller = facet_labeller)

power_upper_limit_shared_all <- FDR_Power_shared_all %>%
  filter(FDR != 0.5) %>%
  ungroup() %>%
  summarise(upper_limit = min(ceiling(max(Power) * 10) / 10 + 0.1, 1)) %>%
  pull(upper_limit)

p_FDR_Power_shared_all <- FDR_Power_shared_fun(
  FDR_Power_shared_all %>% filter(FDR != 0.5)
) +
  ylim(0, power_upper_limit_shared_all) +
  facet_grid(flip ~ rho, labeller = facet_labeller)

ROC_FDR_shared_all <- ggarrange(
  p_ROC_shared_all, p_FDR_Power_shared_all,
  nrow = 2, ncol = 1,
  common.legend = TRUE, legend = "bottom",
  labels = c("a","b"), font.label = list(color = "black", size = 7)
)

ggsave(
  file.path(res_dir, "Figure", "LD_susie_rss_roc_shared_all_flips_rhos_combined.pdf"),
  p_ROC_shared_all,
  height = 210, width = 300, units = "mm", dpi = 500
)


ggsave(
  file.path(res_dir, "Figure", "LD_susie_rss_fdr_power_shared_all_flips_rhos_combined.pdf"),
  p_FDR_Power_shared_all,
  height = 210, width = 300, units = "mm", dpi = 500
)



# -------- Ancestry-specific: ROC + FDR/Power (combined) --------

## --- Split Method into Method + Ancestry for ancestry_ROC_all ---
split_list_ROC <- strsplit(as.character(ancestry_ROC_all$Method), " +")

ancestry_ROC_all <- ancestry_ROC_all %>%
  ungroup() %>%  # make sure no grouping
  mutate(
    Method   = sapply(split_list_ROC, `[`, 1),
    Ancestry = sapply(split_list_ROC, `[`, 2)
  ) %>%
  mutate(
    Method   = fct_relevel(Method,
                           "MESuSiE","SuSiE","Paintor",
                           "SuSiEx","XMAP","MultiSuSiE","CARMAX"),
    Ancestry = fct_relevel(Ancestry, "WB","BB"),
    Ancestry = fct_recode(Ancestry,
                          "White British" = "WB",
                          "Black British" = "BB")
  )

p_ROC_ancestry_all <- ROC_ancestry_fun(ancestry_ROC_all) +
  facet_grid(flip ~ rho, labeller = facet_labeller)

## --- Split Method into Method + Ancestry for FDR_Power_ancestry_all ---
FDR_Power_ancestry_all <- FDR_Power_ancestry_all %>% ungroup()

split_list_FDR <- strsplit(as.character(FDR_Power_ancestry_all$Method), " +")

FDR_Power_ancestry_all <- FDR_Power_ancestry_all %>%
  mutate(
    Method   = sapply(split_list_FDR, `[`, 1),
    Ancestry = sapply(split_list_FDR, `[`, 2)
  ) %>%
  mutate(
    Method   = fct_relevel(Method,
                           "MESuSiE","SuSiE","Paintor",
                           "SuSiEx","XMAP","MultiSuSiE","CARMAX"),
    Ancestry = fct_relevel(Ancestry, "WB","BB"),
    Ancestry = fct_recode(Ancestry,
                          "White British" = "WB",
                          "Black British" = "BB")
  ) %>%
  filter(FDR != 0.5)

power_upper_limit_ancestry_all <- FDR_Power_ancestry_all %>%
  summarise(upper_limit = min(ceiling(max(Power) * 10) / 10 + 0.1, 1)) %>%
  pull(upper_limit)

p_FDR_Power_ancestry_all <- FDR_Power_ancestry_fun(FDR_Power_ancestry_all) +
  ylim(0, power_upper_limit_ancestry_all) +
  facet_grid(flip ~ rho, labeller = facet_labeller)

#ROC_FDR_ancestry_all <- (p_ROC_ancestry_all / p_FDR_Power_ancestry_all) +
#  plot_annotation(tag_levels = "a") &
#  theme(plot.tag = element_text(size = 7, face = "bold"))

#ROC_FDR_ancestry_all <- ROC_FDR_ancestry_all +
#  plot_layout(heights = c(1, 1))

#ggsave(
#  file.path(res_dir, "Figure", "ROC_FDR_Ancestry_all_flips_rhos.pdf"),
#  ROC_FDR_ancestry_all,
#  height = 210, width = 300, units = "mm", dpi = 500
#)

ggsave(
  file.path(res_dir, "Figure", "LD_susie_rss_roc_ancestry_all_flips_rhos_combined.pdf"),
  p_ROC_ancestry_all,
  height = 210, width = 300, units = "mm", dpi = 500
)


ggsave(
  file.path(res_dir, "Figure", "LD_susie_rss_fdr_power_ancestry_all_flips_rhos_combined.pdf"),
  p_FDR_Power_ancestry_all,
  height = 210, width = 300, units = "mm", dpi = 500
)


