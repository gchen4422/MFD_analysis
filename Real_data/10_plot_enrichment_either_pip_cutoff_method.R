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
res_all<-res_all%>%mutate(utr_comb = ifelse((utr_3+utr_5)>0,1,0))




ann_col_name<-c( "non_synonymous","synonymous", "promotor","utr_comb","CRE","heart_ind_eQTL","artery_ind_eQTL","brain_ind_eQTL")




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

# Define thresholds from 0.5 to 0.9 (adjust step as needed)
thresholds <- seq(0.5, 0.9, by = 0.05)

# Define binary annotations (all except CADD and GTEx_eQTL_MaxCPP_common)
binary_ann <- ann_col_name

cont_ann <- c("Z_eur")


# 2) PIP columns
pip_cols <- c(
  MFD = "PIP_Either",
  MESuSiE     = "MESuSiE_PIP_Either",
  SuSiE       = "SuSiE_PIP_Either",
  Paintor     = "Paintor_PIP_Either",
  SuSiEx     = "SuSiEx_PIP_Either",
  XMAP     = "XMAP_PIP_Either",
  CARMAX     = "CARMAX_PIP_Either",
  SuSiE_merged     = "SuSiE_merged_PIP_Either"
  
)


# 4) compute enrichment + continuous means
enrichment_results <- map_dfr(thresholds, function(pip_thresh) {
  map_dfr(names(pip_cols), function(method) {
    pip_col <- pip_cols[[method]]
    
    # signal SNPs & background size
    sig   <- filter(res_all, .data[[pip_col]] > pip_thresh)
    n_sig <- nrow(sig)
    
    # background freq for binary annots
    bg_bin <- res_all %>%
      summarise(across(all_of(binary_ann),
                       ~ sum(.x, na.rm = TRUE) / (n() - n_sig)))
    
    # signal freq for binary annots
    sig_bin <- sig %>%
      summarise(across(all_of(binary_ann),
                       ~ sum(.x, na.rm = TRUE) / n_sig))
    
    # enrichment = signal / background
    enrich_bin <- sig_bin / bg_bin
    
    # mean of continuous annots
    enrich_cont <- sig %>%
      summarise(across(all_of(cont_ann),
                       ~ mean(.x, na.rm = TRUE)))
    
    # stack into long format
    bin_long  <- enrich_bin  %>%
      pivot_longer(everything(),
                   names_to  = "annotation",
                   values_to = "enrichment")
    cont_long <- enrich_cont %>%
      pivot_longer(everything(),
                   names_to  = "annotation",
                   values_to = "enrichment")
    
    bind_rows(bin_long, cont_long) %>%
      mutate(
        method    = method,
        threshold = pip_thresh
      )
  })
})

# Reorder the 'annotation' factor in enrichment_results
enrichment_results <- enrichment_results %>%
  mutate(annotation = factor(annotation, 
                             levels = c("missense",
                                        "synonymous", 
                                        "promotor","utr_comb",
                                        "CRE",
                                        "heart_ind_eQTL",
                                        "artery_ind_eQTL","brain_ind_eQTL","Z_eur"))) %>%
  mutate(annotation = recode(annotation,
                             "missense" = "Missense",
                             "synonymous" = "Synonymous",
                             "promotor" = "Promotor",
                             "utr_comb" = "UTR",
                             "CRE" = "CRE",
                             "heart_ind_eQTL" = "Heart_ind_eQTL",
                             "artery_ind_eQTL" = "Artery_ind_eQTL",
                             "brain_ind_eQTL" = "Brain_ind_eQTL",
                             "Z_eur" = "Z_eur"))

enrichment_results$method <- factor(
  enrichment_results$method, 
  # 1. 'levels' must match the CURRENT strings in the data frame
  levels = c("MFD", "MESuSiE", "SuSiE", "Paintor", "SuSiEx", "XMAP", "CARMAX", "SuSiE_merged"),
  # 2. 'labels' are the names assigned to the levels
  labels = c("MFD", "MESuSiE", "SuSiE", "Paintor", "SuSiEx", "XMAP", "CARMAx", "SuSiE_Meta")
)

# Create the plot with colorblind-friendly colors
enrichment_plot <- ggplot(enrichment_results, aes(x = threshold, y = enrichment, color = method)) +
  geom_line() +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title = "Enrichment Results Across MESuSiE_PIP_Either_all_mlk Cutoffs for BMI",
    x = "MESuSiE_PIP_Either_all_mlk Threshold",
    y = "Enrichment (Signal / Background)"
  ) +
  # Manually set colors for each method
  scale_color_manual(values = 
    
    c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#B2D3A4","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e")
    
  ) +
  theme_bw()

enrichment_plot
plot_dir<-"/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/Figure/"
ggsave(paste0(plot_dir,"enrichment_results_plot_bmi_tidy_mesusiemlk_cadd_Either.png"),enrichment_plot,width=225,height =150,dpi=500,units='mm')






######################################


#Barplot at PIP threshold 0.5


######################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

fill_colors =c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#B2D3A4","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e")

# Create Panel A plot: Enrichment (Signal / Background)
enrichment_panelA = enrichment_results %>% filter(threshold==0.5) %>% filter(annotation != "Z_eur")
plotA <- ggplot(enrichment_panelA, aes(x = method, y = enrichment, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title = "Panel A: Enrichment Results (Signal / Background) at Threshold = 0.5 for BMI",
    x = "Method",
    y = "Enrichment (Signal / Background)",
    color    = "Method",   # legend title for color
  ) +
  scale_fill_manual(name   = "Method", values = fill_colors) +
  theme_bw()+custom_theme() + theme(legend.position = "bottom",    plot.title      = element_blank())

plotA_annot = plotA+ geom_text(label = round(enrichment_panelA %>% arrange(annotation) %>% pull(enrichment),2),position = position_dodge(width = 1),vjust=-0.5,size = 5*5/20)

# Create Panel B plot: Average Score for continuous scores
plotB <- ggplot(enrichment_panelB, aes(x = method, y = score, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ annotation, scales = "free_y") +
  labs(
    title = "Panel B: Average Scores at Threshold = 0.5 for BMI",
    x = "Method",
    y = "Average Score"
  ) +
  scale_fill_manual(values = fill_colors) +
  theme_bw()

# Combine the two panels vertically using patchwork
combined_plot <- plotA / plotB

# Display the combined plot
combined_plot

# Save the combined plot to the designated directory
plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/present/"
system(paste0("mkdir -p ",plot_dir))
ggsave(paste0(plot_dir, "combined_enrichment_and_score_plot_bmi_tidy_mesusiemlk_cadd_Either_bar_noGLM.png"),
       combined_plot, width = 300, height = 300, dpi = 500, units = 'mm')




################################################

# Size and Z-score plot


################################################
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

#p_z = ggplot(data = set_size_z_info_long,aes(x = PHENONAME, y=Z,fill=Method))+geom_bar( stat = "identity",position="dodge")+scale_fill_manual(values=c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#B2D3A4","SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e"))
#p_z = p_z + theme_bw() + xlab("") +ylab("Median |Z|")+ ylim(0,max(round(set_size_z_info_long$Z,2)+1)) +custom_theme()
#p_z = p_z + geom_text(label = round(set_size_z_info_long$Z,2),position = position_dodge(width = 1),vjust=-0.5,size = 5*5/14)

p_z = ggplot(data = set_size_z_info_long, aes(x = PHENONAME, y = Z, fill = Method)) +
  # 1. Be explicit about the width in geom_bar to match it later
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) + 
  scale_fill_manual(values = c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","MFD"="#B2D3A4",
                               "SuSiE_Meta"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3",
                               "SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMAx"="#f2b56e"))

# 2. Move label INSIDE aes() and match the position_dodge width
p_z = p_z + geom_text(
  aes(label = round(Z, 2)),  # Display rounded Z values
  position = position_dodge(width = 0.9), # Match the geom_bar width
  vjust = -0.5,
  size = 5 * 5 /20
)

p_z = p_z + theme_bw() + xlab("") + ylab("Median |Z|") + 
  ylim(0, max(round(set_size_z_info_long$Z, 2) + 1))
p_z = p_z + custom_theme()



p_out<-p_set/p_z+plot_annotation(tag_levels = 'a')+plot_layout(guides = "collect",heights = c(1,1))+  guides(fill = "none", colour = "none", size = "none") +
  theme(legend.position = "none")&theme(legend.position = 'none',plot.tag = element_text(size = 6,face="bold"))
p_out

plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/present/"
ggsave(paste0(plot_dir, "set_size_bmi.png"),
       p_out, width = 200, height = 300, dpi = 500, units = 'mm')



###plotA is barplot 0.5 with first 6 annotations




plotA_annot <- plotA_annot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
# Suppress unused legends
drop_guides <- guides(
  fill = "none",
  color = "none",
  shape = "none",
  plot.title = "none",
  linetype = "none"
)

# 2. Apply to p_out
p_out_no_legend <- p_out + drop_guides

# 3. Assemble the plot
p_total <- (
  p_out_no_legend |
    (plotA_annot + 
       labs(title = NULL) +
       # Apply the 1-row formatting HERE, specifically for plotA
       guides(fill = guide_legend(nrow = 1)) 
    )
) +
  plot_layout(
    guides  = "collect",
    heights = c(0.3, 1),
    widths  = c(1, 1)
  ) +
  plot_annotation(
    title = "Multi-ancestry fine-mapping on All of Us data",
    tag_levels = "a"
  )  & 
  theme(
    plot.title       = element_text(hjust = 0.5),
    plot.tag         = element_text(size = 12),
    legend.position  = "bottom",
    legend.justification = "center",
    legend.box       = "horizontal",
    legend.box.just  = "center",
    strip.text.x     = element_text(size = 7,face = "bold"),
    strip.text.y     = element_text(size = 7,face = "bold"),
    strip.background = element_blank(),
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 7, face = "bold"),
    axis.line        = element_line(color = "black")
  )


p_total


plot_dir <- "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/Figure/"
ggsave(paste0(plot_dir, "set_size_enrichment_pip_all.pdf"),
       p_total, width = 11, height = 6, dpi = 500, units = 'in')









