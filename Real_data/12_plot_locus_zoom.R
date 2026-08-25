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
    axis.text.x = element_text(size = 7.5),
    axis.text.y = element_text(size = 7.5),
    axis.title.x = element_text(size = 10, face="bold"),
    axis.title.y = element_text(size = 10, face="bold"),
    strip.text.x     = element_text(size = 12, face="bold"),   # ← bump this
    strip.text.y     = element_text(size = 12, face="bold"),   # ← and/or this
    strip.background = element_blank(),
    legend.text = element_text(size=10),
    legend.title = element_text(size=10, face="bold"),
    plot.title = element_text(size=14, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.tag = element_text(size = 12),
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
    geom_text_repel(data = data_plot%>%filter(Lead_SNP==1), mapping = aes(x = POS, y = PIP, label = SNP), vjust = 1.2, size = 7/14*5, show.legend =FALSE) 
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
  p_manhattan = p_manhattan + geom_text_repel(data =data_plot%>%filter(Lead_SNP==1), mapping=aes(x=POS, y=PIP, label=SNP),vjust=1.2, size= 7/14*5,show.legend = FALSE)
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
    geom_text(aes(x = Gene, y = Start, label = Gene), hjust = "right", size = 5/14*5) + ylab(paste0("chr",unique(gsub("chr","",gene_list_data$Chrom))))+ xlab("Gene") + 
    theme_bw() +  theme(
      axis.text.x = element_text(size = 7.5),
      axis.text.y = element_blank(),  
      axis.ticks.y =  element_blank(), 
      axis.title.x = element_text(size = 10, face="bold"),
      axis.title.y = element_text(size = 10, face="bold"),
      strip.text.x     = element_text(size = 12, face="bold"),   # ← bump this
      strip.text.y     = element_text(size = 12, face="bold"),   # ← and/or this
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

create_annotation_comparison <- function(candidate_region, 
                                         mafpc_pip_col = "MESuSiE_PIP_Shared_all_mlk",
                                         mesusie_pip_col = "MESuSiE_PIP_Shared",
                                         custom_annotations = NULL) {
  
  # Default key annotations if not provided
  if (is.null(custom_annotations)) {
    key_annotations <- c(
      # Conserved elements
      "Conserved_LindbladToh_common",
      
      # Chromatin accessibility & marks
      "DHS_Trynka_common", 
      "H3K27ac_Hnisz_common", 
      "H3K4me3_Trynka_common",
      "Enhancer_Andersson_common",
      
      # Transcription factor binding
      "TFBS_ENCODE_common",
      
      # Regulatory regions
      "promotor",
      
      # Coding impact (for comparison)
      "Coding_UCSC_common",
      "Intron_UCSC_common",
      "missense",
      "synonymous"
    )
  } else {
    key_annotations <- custom_annotations
  }
  
  # Get top SNPs based on PIP columns
  top_mafpc_snp <- candidate_region %>% 
    slice_max(.data[[mafpc_pip_col]], n = 1) %>%
    pull(SNP)
  
  top_mesusie_snp <- candidate_region %>% 
    slice_max(.data[[mesusie_pip_col]], n = 1) %>%
    pull(SNP)
  
  # Check if annotations exist in the data
  available_annotations <- intersect(key_annotations, colnames(candidate_region))
  if (length(available_annotations) == 0) {
    stop("None of the specified annotations are found in the candidate_region data.")
  }
  
  if (length(available_annotations) < length(key_annotations)) {
    missing_annotations <- setdiff(key_annotations, available_annotations)
    warning(paste("Missing annotations:", paste(missing_annotations, collapse = ", ")))
  }
  
  # Create raw annotation comparison
  annotation_comparison_raw <- candidate_region %>%
    filter(SNP %in% c(top_mafpc_snp, top_mesusie_snp)) %>%
    select(SNP, all_of(available_annotations))
  
  # Create cleaned annotation comparison
  annotation_comparison <- annotation_comparison_raw %>%
    # Convert to long format first
    pivot_longer(cols = -SNP, names_to = "Annotation", values_to = "Present") %>%
    # Clean annotation names
    mutate(
      Annotation = str_replace_all(Annotation, "_", " "),
      Annotation = str_replace_all(Annotation, "common", ""),
      Annotation = str_trim(Annotation),
      Present = as.numeric(Present),  # Convert to 0/1 for heatmap
      SNP_Label = case_when(
        SNP == top_mafpc_snp ~ paste0(SNP, " (MAFPC Top)"),
        SNP == top_mesusie_snp ~ paste0(SNP, " (MESuSiE Top)"),
        TRUE ~ SNP
      )
    )
  
  # Return both raw and processed data along with top SNPs
  return(list(
    annotation_comparison = annotation_comparison,
    annotation_comparison_raw = annotation_comparison_raw,
    top_mafpc_snp = top_mafpc_snp,
    top_mesusie_snp = top_mesusie_snp,
    available_annotations = available_annotations,
    mafpc_pip = candidate_region %>% filter(SNP == top_mafpc_snp) %>% pull(.data[[mafpc_pip_col]]),
    mesusie_pip = candidate_region %>% filter(SNP == top_mesusie_snp) %>% pull(.data[[mesusie_pip_col]])
  ))
}



####################################
#
# For shared causal variants for DBP
#
####################################


load(paste0(res_dir,"res_all_mf.RData"))   
res_all = res_all %>% filter(PHENONAME == "BMI")

# 1) set threshold
pip_thr <- 0.5

# 2) define column names
mfd_shared_col <- "PIP_Either"
mesusie_col      <- "MESuSiE_PIP_Either"
SuSiE_merged_col      <- "SuSiE_merged_PIP_Either"

# 3) filter for novel shared hits
novel_shared <- res_all %>%
  filter(
    # called by MAFPC as shared
    .data[[mfd_shared_col]] > pip_thr,
    # NOT called by MESuSiE
    .data[[mesusie_col]]        > pip_thr,
    # NOT called by Paintor
    .data[[SuSiE_merged_col]]        <= pip_thr
  ) %>%
  select(Region, SNP, Z_eur,Z_afr,all_of(c(mfd_shared_col, mesusie_col, SuSiE_merged_col))) %>% filter(abs(Z_eur)>6|abs(Z_afr)>6)

# 4) get the set of regions with ≥1 novel shared SNP
novel_regions <- novel_shared %>%
  pull(Region) %>%
  unique()

# 5) optional: a summary table
novel_summary <- novel_shared %>%
  group_by(Region) %>%
  summarise(
    n_novel_snps = n(),
    top_novel_pip = max(.data[[mfd_shared_col]])
  ) %>%
  arrange(desc(n_novel_snps), desc(top_novel_pip))

# inspect
head(novel_shared)
head(novel_regions)
head(novel_summary)




####################################
#
# For WB causal variants for DBP
#
####################################

load(paste0(res_dir,"res_all_mf.RData"))   
res_all = res_all %>% filter(PHENONAME == "DBP")


# 1) set threshold
pip_thr <- 0.5

# 2) define column names
mafpc_col <- "MESuSiE_PIP_WB_all_mlk"
mesusie_col      <- "MESuSiE_PIP_WB"
paintor_col      <- "Paintor_PIP"

# 3) filter for novel shared hits
novel_wb <- res_all %>%
  filter(
    # called by MAFPC as shared
    .data[[mafpc_col]] > pip_thr,
    # NOT called by MESuSiE
    .data[[mesusie_col]]        <= pip_thr,
    # NOT called by Paintor
    .data[[paintor_col]]        <= pip_thr
  ) %>%
  select(Region, SNP, all_of(c(mafpc_col, mesusie_col, paintor_col)))

# 4) get the set of regions with ≥1 novel shared SNP
novel_regions <- novel_wb %>%
  pull(Region) %>%
  unique()

# 5) optional: a summary table
novel_summary <- novel_wb %>%
  group_by(Region) %>%
  summarise(
    n_novel_snps = n(),
    top_novel_pip = max(.data[[mafpc_col]])
  ) %>%
  arrange(desc(n_novel_snps), desc(top_novel_pip))

# inspect
head(novel_wb)
head(novel_regions)
head(novel_summary)



####################################
#
# For BB causal variants for DBP
#
####################################

load(paste0(res_dir,"res_bp_pca_bic_modified_llk_updated.RData"))
res_all = res_all %>% filter(Trait == "DBP")


# 1) set threshold
pip_thr <- 0.5

# 2) define column names
mafpc_col <- "MESuSiE_PIP_BB_all_mlk"
mesusie_col      <- "MESuSiE_PIP_BB"
paintor_col      <- "Paintor_PIP"

# 3) filter for novel shared hits
novel_bb <- res_all %>%
  filter(
    # called by MAFPC as shared
    .data[[mafpc_col]] > pip_thr,
    # NOT called by MESuSiE
    .data[[mesusie_col]]        <= pip_thr,
    # NOT called by Paintor
    .data[[paintor_col]]        <= pip_thr
  ) %>%
  select(Region, SNP, all_of(c(mafpc_col, mesusie_col, paintor_col)))

# 4) get the set of regions with ≥1 novel shared SNP
novel_regions <- novel_bb %>%
  pull(Region) %>%
  unique()

# 5) optional: a summary table
novel_summary <- novel_bb %>%
  group_by(Region) %>%
  summarise(
    n_novel_snps = n(),
    top_novel_pip = max(.data[[mafpc_col]])
  ) %>%
  arrange(desc(n_novel_snps), desc(top_novel_pip))

# inspect
head(novel_bb)
head(novel_regions)
head(novel_summary)



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
EUR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = -log10(2*pnorm(-abs(Z_eur))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
AFR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = -log10(2*pnorm(-abs(Z_afr))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
p_EUR<-gwas_plot_fun (EUR_GWAS_plot_data, "All of Us EA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")
p_AFR<-gwas_plot_fun (AFR_GWAS_plot_data, "All of Us BA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")

###Finemap Plot
#EUR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = susie_EU,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
#AFR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = susie_BB,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
Paintor_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
MESuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MESuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
MFD_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MFD_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)


#p_EUR_SuSiE<-finemap_plot_fun(EUR_SuSiE_plot_data, "SuSiE EUR", "PIP", 0.5)
#p_AFR_SuSiE<-finemap_plot_fun(AFR_SuSiE_plot_data, "SuSiE AFR", "PIP", 0.5)
p_Paintor<-finemap_plot_fun(Paintor_plot_data, "Paintor", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_MESuSiE<-finemap_plot_fun(MESuSiE_plot_data, "MESuSiE", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_MFD<-finemap_plot_fun(MFD_plot_data, "MFD", "PIP", 0.5)#+ drop_guides + theme(legend.position = "none")

# Gene Plot
plot.range <- c(min(candidate_region$POS), max(candidate_region$POS))
Gene_List_sub_coding<-Gene_List%>%filter(Chrom==paste0("chr",unique(candidate_region$CHR)))%>%filter(Start<max(candidate_region$POS),End>min(candidate_region$POS))%>%filter(Coding=="proteincoding")%>%filter(!is.na(cdsLength))%>%filter(GeneLength>15000)
p2<-gene_range_plot_fun(Gene_List_sub_coding,plot.range)

##Combine Plot together
combined_plot<-(p_EUR/p_MESuSiE/p2+plot_layout(heights = c(1,1,1))|p_AFR/p_MFD/p2+plot_layout(heights = c(1,1,1)))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")
combined_plot
ggsave(paste0(plot_dir,"Region_",LD_BLOCK ,"_",trait,".pdf"),combined_plot,width = 200, height = 180,units="mm",dpi=500)


annotation_results <- create_annotation_comparison(candidate_region,"MESuSiE_PIP_Either_all_mlk")


annotation_colors <- c(
  "Coding_UCSC_common"   = "#0072B2", # A nice blue
  "H3K27ac_Hnisz_common" = "#FFB6C1", # A light pink
  "Promoter_UCSC_common" = "#009E73", # A bluish green
  "default"              = "gray30"  # A dark gray for other annotations
)

plot_annotation_highlight <- function(data, 
                                      annotation_col, 
                                      highlight_snp_id,
                                      snp_col = "SNP",
                                      pos_col = "POS",
                                      chr_col = "CHR",
                                      plot_range = NULL,
                                      custom_theme = NULL) {
  
  # --- 1. PREPARE VARIABLES ---
  
  # Select the color for the '1' values from the palette
  track_color <- annotation_colors[[annotation_col]]
  if (is.null(track_color)) {
    track_color <- annotation_colors[["default"]]
  }
  
  plot_data <- data[data[[annotation_col]] == 1, ]
  
  
  # Convert annotation column to a factor for discrete coloring
  plot_data[[annotation_col]] <- as.factor(plot_data[[annotation_col]])
  
  chromosome_label <- paste0("chr", unique(plot_data[[chr_col]])[1])
  
  if (is.null(plot_range)) {
    pos_rng <- range(data[[pos_col]], na.rm = TRUE)
  } else {
    pos_rng <- plot_range
  }
  
  
  highlight_snp <- data[data[[snp_col]] == highlight_snp_id, ]
  
  # --- 2. BUILD THE PLOT ---
  
  p <- ggplot(data = plot_data, aes(x = .data[[pos_col]], y = 1)) +
    # This layer's color now depends on the value (0 or 1)
    geom_point(
      aes(color = .data[[annotation_col]]), # Color is now an aesthetic
      shape = "|", 
      size = 5, 
      alpha = 0.7,
      show.legend = FALSE # Hide the color legend
    ) +
    
    # Manually set the colors for '0' and '1'
    scale_color_manual(values = c("1" = track_color)) +
    
    # Layer to Highlight the specified SNP
    geom_point(
      data = highlight_snp,
      color = "red",
      size = 1.6
    ) +
    
    # Set the genomic position axis
    scale_x_continuous(
      limits = pos_rng,
      labels = function(x) paste0(x / 1e6, " MB")
    ) +
    
    # NEW: Force the y-axis to only show labels for 0 and 1
    scale_y_discrete(breaks = c(0, 1)) +
    
    labs(
      y = annotation_col, 
      x = chromosome_label
    ) +
    theme_bw()
  
  if (!is.null(custom_theme) && is.function(custom_theme)) {
    p <- p + custom_theme()
  }
  
  return(p)
}

p_coding = plot_annotation_highlight(candidate_region, "Coding_UCSC_common", "rs3184504",custom_theme = custom_theme() ) + custom_theme() +    theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.y = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x =   element_blank(),
  axis.line.x = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank()
)

p_H3K27ac= plot_annotation_highlight(candidate_region, "H3K27ac_Hnisz_common", "rs3184504",custom_theme = custom_theme() ) + custom_theme()+    theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.y = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x =   element_blank(),
  axis.line.x = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank()
)


combined_plot<-(p_EUR/p_coding/p_H3K27ac/p_MAFPC/p2+plot_layout(heights = c(1,0.1,0.1,1,1))|p_AFR/p_coding/p_H3K27ac/p_MESuSiE/p2+plot_layout(heights = c(1,0.1,0.1,1,1)))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")

ggsave(paste0(plot_dir,"Region_",LD_BLOCK ,"_",trait,"_annoted.pdf"),combined_plot,width = 7.42, height = 210,units="mm",dpi=300)

plot_dir = "/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/format_summary_stats/Figure/Poster_ashg_2025/"
ggsave(paste0(plot_dir,"Region_",LD_BLOCK ,"_",trait,"_annoted.png"),combined_plot,width = 11.02, height = 7.42,dpi=600)
ggsave(paste0(plot_dir,"Region_",LD_BLOCK ,"_",trait,"_annoted.pdf"),combined_plot,width = 11.02, height = 7.42,dpi=600)



###################################################################
#
#   LOCI 18 (rs7582535), LOCI 23 (rs1861151), LOCI 52 (rs9825254) LOCI 164 (rs4776304) 170 (rs201934469) for BMI is an example of how MFD(susie-post hoc) can find AS-V but mesusie cannot
#   LOCI 115 (rs56663592) LOCI 145 (rs1265564)  180 (rs1317867)
#   LOCI 31 (rs13400734) for BMI is an example of how MFD(susie-post hoc) can find AS-V but mesusie cannot
#
##################################################################

drop_guides <- guides(
  color = "none",
  shape = "none",
  fill  = "none",
  linetype = "none"
)


LD_BLOCK = 115
trait = "BMI"
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
lead_SNP = "rs56663592"
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
candidate_region%>%filter(SNP=="rs56663592")%>%summarise(2*pnorm(-abs(Z_eur)),2*pnorm(-abs(Z_afr)))
candidate_region%>%mutate(Pvalue_BB =2*pnorm(-abs(Z_afr)) )%>%filter(Pvalue_BB<5e-8)
candidate_region %>%
  filter(SNP == "rs56663592") %>%
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
EUR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = -log10(2*pnorm(-abs(Z_eur))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
AFR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = -log10(2*pnorm(-abs(Z_afr))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
p_EUR<-gwas_plot_fun (EUR_GWAS_plot_data, "All of Us EA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")
p_AFR<-gwas_plot_fun (AFR_GWAS_plot_data, "All of Us BA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")

###Finemap Plot
#EUR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = susie_EU,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
#AFR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = susie_BB,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
Paintor_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
SuSiE_merged_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = SuSiE_merged_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_merged_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
MESuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MESuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
MFD_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MFD_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)


#p_EUR_SuSiE<-finemap_plot_fun(EUR_SuSiE_plot_data, "SuSiE EUR", "PIP", 0.5)
#p_AFR_SuSiE<-finemap_plot_fun(AFR_SuSiE_plot_data, "SuSiE AFR", "PIP", 0.5)
p_Paintor<-finemap_plot_fun(Paintor_plot_data, "Paintor", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_SuSiE_merged<-finemap_plot_fun(SuSiE_merged_plot_data, "SuSiE_Meta", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")

p_MESuSiE<-finemap_plot_fun(MESuSiE_plot_data, "MESuSiE", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_MFD<-finemap_plot_fun(MFD_plot_data, "MFD", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")

# Gene Plot
plot.range <- c(min(candidate_region$POS), max(candidate_region$POS))
Gene_List_sub_coding<-Gene_List%>%filter(Chrom==paste0("chr",unique(candidate_region$CHR)))%>%filter(Start<max(candidate_region$POS),End>min(candidate_region$POS))%>%filter(Coding=="proteincoding")%>%filter(!is.na(cdsLength))%>%filter(GeneLength>15000)
p2<-gene_range_plot_fun(Gene_List_sub_coding,plot.range)

##Combine Plot together
combined_plot<-(p_EUR/p_MESuSiE/p_Paintor/p2+plot_layout(heights = c(1,1,1,1.5))|p_AFR/p_MFD/p_SuSiE_merged/p2+plot_layout(heights = c(1,1,1,1.5)))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")
combined_plot
ggsave(paste0(plot_dir,"Region_",LD_BLOCK ,"_",trait,".pdf"),combined_plot,width = 180, height = 210,units="mm",dpi=300)



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
candidate_region%>%filter(SNP=="rs201934469")%>%summarise(2*pnorm(-abs(Z_eur)),2*pnorm(-abs(Z_afr)))
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

#candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(n())
#candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(Z_eur))
#candidate_region%>%filter(abs(candidate_region$r2_EUR)>0.9)%>%summarise(summary(2*pnorm(-abs(Z_eur))))

candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise(summary(Z_afr))
candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise((2*pnorm(-abs(Z_afr))))
candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%summarise(r2_AFR)

candidate_region%>%filter(abs(candidate_region$r2_AFR)>0.9)%>%pull(Z_afr)


#annotation_results <- create_annotation_comparison(candidate_region)


###GWAS PLOT
EUR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = -log10(2*pnorm(-abs(Z_eur))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
AFR_GWAS_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = -log10(2*pnorm(-abs(Z_afr))),Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
p_EUR<-gwas_plot_fun (EUR_GWAS_plot_data, "All of Us EA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")
p_AFR<-gwas_plot_fun (AFR_GWAS_plot_data, "All of Us BA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")

###Finemap Plot
#EUR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = susie_EU,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
#AFR_SuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_AFR,PIP = susie_BB,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
Paintor_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = Paintor_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(Paintor_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
SuSiE_merged_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = SuSiE_merged_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(SuSiE_merged_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
MESuSiE_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = MESuSiE_PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MESuSiE_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)
MFD_plot_data<-candidate_region%>%mutate(r2 = r2_EUR,PIP = PIP_Either,Lead_SNP = ifelse(SNP==lead_SNP,1,0),POS= as.numeric(POS),cat = factor(MFD_cat,levels = c("0", "1", "2", "3", "4"), labels = c("Non", "EUR", "AFR", "Shared", "Paintor")))%>%select(SNP,POS, r2,PIP,Lead_SNP,cat)


#p_EUR_SuSiE<-finemap_plot_fun(EUR_SuSiE_plot_data, "SuSiE EUR", "PIP", 0.5)
#p_AFR_SuSiE<-finemap_plot_fun(AFR_SuSiE_plot_data, "SuSiE AFR", "PIP", 0.5)
p_Paintor<-finemap_plot_fun(Paintor_plot_data, "Paintor", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_SuSiE_merged<-finemap_plot_fun(SuSiE_merged_plot_data, "SuSiE_Meta", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")

p_MESuSiE<-finemap_plot_fun(MESuSiE_plot_data, "MESuSiE", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")
p_MFD<-finemap_plot_fun(MFD_plot_data, "MFD", "PIP", 0.5)+ drop_guides + theme(legend.position = "none")

# Gene Plot
plot.range <- c(min(candidate_region$POS), max(candidate_region$POS))
Gene_List_sub_coding<-Gene_List%>%filter(Chrom==paste0("chr",unique(candidate_region$CHR)))%>%filter(Start<max(candidate_region$POS),End>min(candidate_region$POS))%>%filter(Coding=="proteincoding")%>%filter(!is.na(cdsLength))%>%filter(GeneLength>10000)
p2<-gene_range_plot_fun(Gene_List_sub_coding,plot.range)

##Combine Plot together
combined_plot<-(p_EUR/p_MESuSiE/p_Paintor/p2+plot_layout(heights = c(1,1,1,1.5))|p_AFR/p_MFD/p_SuSiE_merged/p2+plot_layout(heights = c(1,1,1,1.5)))+plot_layout(guides = 'collect')&theme(legend.position = "bottom")
combined_plot
ggsave(paste0(plot_dir,"Region_",LD_BLOCK ,"_",trait,".pdf"),combined_plot,width = 180, height = 210,units="mm",dpi=300)


###############


test = fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/results_DBP/result/Metal_DBP_result1.tbl")

metal_results_loci_raw <- test %>% filter(MarkerName %in% candidate_region$CHR_POS) %>% 
  separate(MarkerName, into = c("CHR","POS"), sep = "_", remove = FALSE, convert = TRUE) %>% arrange(CHR, POS) 

chr = 15
start = 90361425
end = 91391644
lead_SNP = "rs9821489"
test_dbp_loci_31 = metal_results_loci_raw %>% filter(CHR == chr) %>% filter(POS > start, POS < end)

AFR_GWAS_plot_data<-test_dbp_loci_31%>%mutate(SNP = rsid,r2 = 0,PIP = -log10(p_value),Lead_SNP = ifelse(rsid==lead_SNP,1,0),POS= as.numeric(base_pair_location))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
p_AFR<-gwas_plot_fun (AFR_GWAS_plot_data, "All of Us BA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")


#################



test = fread("/scratch/negishi/chen4422/hapnest/multi-ans-sum-stats/Lipids_all_of_us_v8/formatted/test_MF_pipleline/summstats_to_finemap/dbp_afr_verma_2024.tsv.gz")
chr = 3
start = 53041701
end = 54041701
lead_SNP = "rs9821489"
test_dbp_loci_11 = test %>% filter(chromosome == 3) %>% filter(base_pair_location > start, base_pair_location < end)

AFR_GWAS_plot_data<-test_dbp_loci_11%>%mutate(SNP = rsid,r2 = 0,PIP = -log10(p_value),Lead_SNP = ifelse(rsid==lead_SNP,1,0),POS= as.numeric(base_pair_location))%>%dplyr::select(SNP,POS, r2,PIP,Lead_SNP)
p_AFR<-gwas_plot_fun (AFR_GWAS_plot_data, "All of Us BA", "-log10(P-value)", -log10(5e-8))+ drop_guides + theme(legend.position = "none")
