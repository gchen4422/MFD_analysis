custom_theme <- function() {
  theme_bw() +  # Start with the standard black & white theme
    theme(
      # 1. Remove all grids
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Configure text sizes and styles
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      axis.title.x = element_text(size = 9, face = "bold"),
      axis.title.y = element_text(size = 9, face = "bold"),
      strip.text.x = element_text(size = 9),
      strip.text.y = element_text(size = 9),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 9, hjust = 0.5)
    )
}


custom_theme_supp <- function() {
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),  
    axis.title.x = element_text(size = 7, face="bold"),
    axis.title.y = element_text(size = 7, face="bold"),
    strip.text.x = element_text(size = 7),
    strip.text.y = element_text(size = 7),
    strip.background = element_blank(),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7, face="bold"),
    plot.title = element_text(size = 7, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line = element_line(color = "black")
  )
}

#####################################################################################################################
#
#
#                     Compare result of BB = 50000 versus BB = 300,000/ In sample LD versus External LD
#
#
#####################################################################################################################


FDR_Power_Compare_fun_supp<-function(FDR_Power,pattern_input){
  p = ggplot(FDR_Power, aes(FDR, Power,fill = Method,pattern = Pattern)) + 
    geom_col_pattern(position = "dodge", pattern_density = .01, pattern_spacing = 0.04,pattern_key_scale_factor = 0.3,pattern_fill ="#DEDBD2",pattern_colour ="#d8c99b" ) +
    scale_fill_manual(name = "Method",values = c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","SuSiE_weighted"="#B2D3A4","SuSiE_merged"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMA-X"="#f2b56e","MFD" = "#377eb8","MFD_MESuSiE"="#5A7EC0","MFD_SuSiEx"="#9467BD","MFD_XMAP"="#2D5E9E"),guide=FALSE)+
    scale_pattern_manual(values=pattern_input)+
    facet_grid(vars(h2),vars(causal_num),labeller=label_parsed)+
    ylab("Power")+xlab("FDR")+
    theme_bw() + custom_theme_supp() + theme(legend.position = "bottom")+
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2)+
    guides(pattern = guide_legend(override.aes = list(fill = "white", order = 1)))
  return(p)
}

FDR_Power_Compare_fun <- function(FDR_Power, pattern_input) {
  
  # CRITICAL: Arrange data first so the underlying text order matches the plot order
  # Sort by Method (outer group) then Pattern (inner group)
  require(dplyr)
  FDR_Power <- FDR_Power %>% arrange(Method, Pattern)
  
  # CRITICAL: Define 'group' in the main ggplot call.
  # group = interaction(Pattern, Method) ensures:
  # 1. 'Pattern' varies fast (inner loop), 'Method' varies slow (outer loop).
  #    This keeps the 3 patterns for Method A together, then Method B, etc.
  # 2. Both geom_col_pattern and geom_text use this EXACT same grouping.
  
  p <- ggplot(FDR_Power, aes(x = Cat, y = Power, 
                             fill = Method, 
                             pattern = Pattern,
                             group = interaction(Pattern, Method))) + 
    
    # 1. Bars
    geom_col_pattern(
      position = position_dodge(width = 0.9), 
      pattern_density = 0.01, 
      pattern_spacing = 0.02,
      pattern_key_scale_factor = 0.3,
      pattern_fill = "#DEDBD2",
      pattern_colour = "#d8c99b"
    ) +
    
    # 2. Text
    # The group mapping is inherited from ggplot().
    # It will automatically follow the exact same dodge positions as the bars.
    geom_text(
      aes(label = round(Power, 3)),
      position = position_dodge(width = 0.9), 
      vjust = -0.5,
      size = 1.5
    ) +
    
    # 3. Colors and Scales
    scale_fill_manual(
      name = "Method",
      values = c("MESuSiE"="#8da0cb", "SuSiE"="#66c2a5", "SuSiE_weighted"="#B2D3A4",
                 "SuSiE_merged"="#9FBA95", "Paintor"="#fc8d62", "MultiSuSiE"="#e78ac3",
                 "SuSiEx"="#E89DA0", "XMAP"="#ffd92f", "CARMA-X"="#f2b56e"),
      guide = "none"
    ) +
    
    scale_pattern_manual(values = pattern_input) +
    
    # 4. Axes and Themes
    ylab("Power") + 
    xlab("FDR = 0.05") +
    ylim(0, 0.3) + 
    
    theme_bw() + 
    custom_theme() + 
    theme(legend.position = "bottom") +
    
    # 5. Border and Guides
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2) +
    guides(pattern = guide_legend(override.aes = list(fill = "white", order = 1)))
  
  return(p)
}

FDR_Power_Compare_fun_noannot <- function(FDR_Power, pattern_input) {
  
  # CRITICAL: Arrange data first so the underlying text order matches the plot order
  # Sort by Method (outer group) then Pattern (inner group)
  require(dplyr)
  FDR_Power <- FDR_Power %>% arrange(Method, Pattern)
  
  # CRITICAL: Define 'group' in the main ggplot call.
  # group = interaction(Pattern, Method) ensures:
  # 1. 'Pattern' varies fast (inner loop), 'Method' varies slow (outer loop).
  #    This keeps the 3 patterns for Method A together, then Method B, etc.
  # 2. Both geom_col_pattern and geom_text use this EXACT same grouping.
  
  p <- ggplot(FDR_Power, aes(x = Cat, y = Power, 
                             fill = Method, 
                             pattern = Pattern,
                             group = interaction(Pattern, Method))) + 
    
    # 1. Bars
    geom_col_pattern(
      position = position_dodge(width = 0.9), 
      pattern_density = 0.01, 
      pattern_spacing = 0.02,
      pattern_key_scale_factor = 0.3,
      pattern_fill = "#DEDBD2",
      pattern_colour = "#d8c99b"
    ) +
    
    # 2. Text
    # The group mapping is inherited from ggplot().
    # It will automatically follow the exact same dodge positions as the bars.
    #geom_text(
    #  aes(label = round(Power, 3)),
    #  position = position_dodge(width = 0.9), 
    #  vjust = -0.5,
    #  size = 1.5
    #) +
    
    # 3. Colors and Scales
    scale_fill_manual(
      name = "Method",
      values = c("MESuSiE"="#8da0cb", "SuSiE"="#66c2a5", "SuSiE_weighted"="#B2D3A4",
                 "SuSiE_merged"="#9FBA95", "Paintor"="#fc8d62", "MultiSuSiE"="#e78ac3",
                 "SuSiEx"="#E89DA0", "XMAP"="#ffd92f", "CARMA-X"="#f2b56e", "MFD" = "#377eb8","MFD_MESuSiE"="#5A7EC0","MFD_SuSiEx"="#9467BD","MFD_XMAP"="#2D5E9E"),
      guide = "none"
    ) +
    
    scale_pattern_manual(values = pattern_input) +
    
    # 4. Axes and Themes
    ylab("Power") + 
    xlab("FDR = 0.05") +
    ylim(0, 0.3) + 
    
    theme_bw() + 
    custom_theme() + 
    theme(legend.position = "bottom") +
    
    # 5. Border and Guides
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2) +
    guides(pattern = guide_legend(override.aes = list(fill = "white", order = 1)))
  
  return(p)
}

ROC_ancestry_Compare_fun <- function(all_ROC_data_dataframe, pattern_input) {
  
  # 1. Define Base Colors
  base_colors <- c("MESuSiE"="#8da0cb", "SuSiE"="#66c2a5", "SuSiE_weighted"="#B2D3A4",
                   "SuSiE_merged"="#9FBA95", "Paintor"="#fc8d62", "MultiSuSiE"="#e78ac3",
                   "SuSiEx"="#E89DA0", "XMAP"="#ffd92f", "CARMA-X"="#f2b56e", "MFD" = "#377eb8","MFD_MESuSiE"="#5A7EC0","MFD_SuSiEx"="#9467BD","MFD_XMAP"="#2D5E9E")
  
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



##################################################
#
#         Set Size Plot
#
##################################################
Set_Size_Compare_fun <- function(all_Set_data_dataframe, upper_limit, pattern_input_Power) {
  
  # 1. Define the exact order
  method_order <- c("MFD","MESuSiE", "SuSiE", "SuSiE_weighted", "SuSiE_merged", 
                    "Paintor", "MultiSuSiE", "SuSiEx", "XMAP", "CARMA-X")
  
  # 2. Plot
  p_size_box <- ggplot(data = all_Set_data_dataframe, 
                       aes(x = factor(Method, levels = method_order), # Apply method order
                           y = Size, 
                           fill = Method, 
                           pattern = Pattern)) +
    
    # Use geom_boxplot_pattern (handles both color fill and pattern)
    geom_boxplot_pattern(
      position = position_dodge(preserve = "single"), 
      color = "black", 
      size = 0.2,         # Adjusted for visibility
      outlier.size = 0.1,
      pattern_density = 0.2,
      pattern_size = 0.1,
      pattern_key_scale_factor = 0.3,
      pattern_spacing = 0.05,
      pattern_fill = "white", 
      pattern_color = "white"
    ) +
    
    # Coordinates
    coord_cartesian(ylim = c(0, upper_limit)) +
    
    # Colors and Patterns
    scale_fill_manual(
      values = c("MESuSiE"="#8da0cb", "SuSiE"="#66c2a5", "SuSiE_weighted"="#B2D3A4", 
                 "SuSiE_merged"="#9FBA95", "Paintor"="#fc8d62", "MultiSuSiE"="#e78ac3", 
                 "SuSiEx"="#E89DA0", "XMAP"="#ffd92f", "CARMA-X"="#f2b56e","MFD" = "#377eb8","MFD_MESuSiE"="#5A7EC0","MFD_SuSiEx"="#9467BD","MFD_XMAP"="#2D5E9E"),
      guide = "none" # Replaced deprecated 'guide = FALSE'
    ) +
    
    scale_pattern_manual(values = pattern_input_Power, guide = "none") +
    
    # Facets and Labels
    ylab("Set Size") +
    xlab("") +
    facet_grid(h2 ~ causal_num, labeller = label_parsed, scales = "free", space = "free") +
    
    # Theme
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 7, face = "bold"),
      axis.text.y = element_text(size = 7),  
      axis.title.x = element_text(size = 7, face = "bold"),
      axis.title.y = element_text(size = 7, face = "bold"),
      strip.text.x = element_text(size = 7),
      strip.text.y = element_text(size = 7),
      strip.background = element_blank(),
      legend.position = "none", # Consolidated removal of legends
      plot.title = element_text(size = 9, hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(), 
      axis.line = element_line(color = "black")
    ) + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2)
  
  return(p_size_box)
}


##################################################
#
#         Set Power Plot
#
##################################################
Set_Power_Compare_fun<-function(power_summary,pattern_input_Power){
  p_power_bar = ggplot(data=power_summary,aes(x = Method, y =Power_name,fill=Method,pattern =Pattern))+
    geom_col_pattern(position = "dodge",color = "#212529", size=0.05,aes(fill = Method,pattern = Pattern), pattern_density  = 0.2,pattern_size = 0.1,pattern_key_scale_factor = 0.3,pattern_spacing=0.05,pattern_fill="white", pattern_color="white") +
    scale_fill_manual(values=c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","SuSiE_weighted"="#B2D3A4","SuSiE_merged"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMA-X"="#f2b56e","MFD" = "#377eb8","MFD_MESuSiE"="#5A7EC0","MFD_SuSiEx"="#9467BD","MFD_XMAP"="#2D5E9E"),guide=FALSE)+
    scale_pattern_manual(values = pattern_input_Power)+
    guides(
      fill = guide_legend(
        title = "Method",
        order = 1,                            # Places Method legend first
        override.aes = list(pattern = "none") # Removes texture from the color legend keys
      ),
      pattern = guide_legend(
        override.aes = list(
          fill          = "white",         # white background
          colour        = NA,              # no outline box
          pattern_fill  = "#DEDBD2",
          pattern_color = "#d8c99b"
        )
      )
    ) +
    ylab("Power")+xlab("")+facet_grid(h2~causal_num,labeller=label_parsed)+ylim(c(0,1))+theme_bw()
  
  #p_power_bar = p_power_bar+geom_text(aes(x=Method,y=Power_name,label=round(Power_name,2)),position = position_dodge(width = 1),vjust=-0.5,size = 5*5/14)
  
  p_power_bar = p_power_bar +  theme(
    axis.text.x = element_text(size = 7,face="bold"),
    axis.text.y = element_text(size = 7),  
    axis.title.x = element_text(size = 7, face="bold"),
    axis.title.y = element_text(size = 7, face="bold"),
    strip.text.x = element_text(size = 7),
    strip.text.y= element_text(size = 7),
    # strip.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_blank(),
    legend.text=element_text(size=7),
    legend.title=element_text(size=7, face="bold"),
    plot.title = element_text(size=9, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line = element_line(color = "black")
  )
  
  p_power_bar = p_power_bar + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black" ,size = 0.2)
  return(p_power_bar)
}


Set_Power_Compare_fun_supp <- function(power_summary, pattern_input_Power) {
  
  # Optional: Define order if not already factorized in data
  method_order <- c("MFD","MESuSiE", "SuSiE", "SuSiE_weighted", "SuSiE_merged", 
                    "Paintor", "MultiSuSiE", "SuSiEx", "XMAP", "CARMA-X")
  
  p_power_bar = ggplot(data = power_summary, 
                       aes(x = factor(Method, levels = method_order), # Ensure order
                           y = Power_name, 
                           fill = Method, 
                           pattern = Pattern)) +
    
    # 1. Bar Pattern Setup
    geom_col_pattern(
      position = "dodge",
      color = "#212529", 
      size = 0.05,
      pattern_density = 0.2,
      pattern_size = 0.1,
      pattern_key_scale_factor = 0.3,
      pattern_spacing = 0.05,
      pattern_fill = "white", 
      pattern_color = "white"
    ) +
    
    # 2. Colors without guide=FALSE
    scale_fill_manual(
      values = c("MESuSiE"="#8da0cb", "SuSiE"="#66c2a5", "SuSiE_weighted"="#B2D3A4",
                 "SuSiE_merged"="#9FBA95", "Paintor"="#fc8d62", "MultiSuSiE"="#e78ac3",
                 "SuSiEx"="#E89DA0", "XMAP"="#ffd92f", "CARMA-X"="#f2b56e", "MFD" = "#377eb8","MFD_MESuSiE"="#5A7EC0","MFD_SuSiEx"="#9467BD","MFD_XMAP"="#2D5E9E")
    ) +
    
    scale_pattern_manual(values = pattern_input_Power) +
    
    # 3. Legend configuration
    guides(
      fill = guide_legend(
        title = "Method",
        nrow = 2,               # Display methods in two rows
        order = 1,
        override.aes = list(pattern = "none") # Removes texture from color keys
      ),
      pattern = guide_legend(
        title = "Condition",    # Good practice to have a title
        nrow = 2,               # Matches height of Method legend
        order = 2,
        override.aes = list(
          fill = "white", 
          colour = NA, 
          pattern_fill = "#DEDBD2", 
          pattern_color = "#d8c99b"
        )
      )
    ) +
    
    # 4. Labels & Facets
    ylab("Power") + xlab("") +
    facet_grid(h2 ~ causal_num, labeller = label_parsed) +
    ylim(c(0, 1)) +
    theme_bw() +
    
    # 5. Theme
    theme(
      axis.text.x = element_text(size = 7, face = "bold"),
      axis.text.y = element_text(size = 7),  
      axis.title.x = element_text(size = 7, face = "bold"),
      axis.title.y = element_text(size = 7, face = "bold"),
      strip.text.x = element_text(size = 7),
      strip.text.y = element_text(size = 7),
      strip.background = element_blank(),
      
      # Legend Styling
      legend.position = "bottom", # Ensure it's at the bottom
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7, face = "bold"),
      
      plot.title = element_text(size = 9, hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(), 
      axis.line = element_line(color = "black")
    ) +
    
    # Panel border
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf, colour = "black", size = 0.2)
  
  return(p_power_bar)
}
##############################################################################
#
#
#                       Missing causal SNP
#
#
##############################################################################

####################################################################################
#
#
#			Distance
#
#
#
####################################################################################
###Signal Distance and LD (MESuSiE, SuSiE Either, Paintor)
Signal_Distance_Fun<-function(Signal_distance_data){
  p_signal_distance = ggplot(data=Signal_distance_data,aes(x = Method, y = Distance,fill=Method))+geom_bar(stat = "identity", position = "dodge", aes(fill=Method),alpha = 1.2)+scale_fill_manual(values=c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","SuSiE_weighted"="#B2D3A4","SuSiE_merged"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMA-X"="#f2b56e","MFD" = "#377eb8","MFD_MESuSiE"="#5A7EC0","MFD_SuSiEx"="#9467BD","MFD_XMAP"="#2D5E9E"),guide = FALSE)+ylab("Distance to true signal")+xlab("")+facet_grid(~h2,labeller=label_parsed)+theme_bw()
  p_signal_distance = p_signal_distance +  theme(
    axis.text.x = element_text(size = 5,face="bold"),
    axis.text.y = element_text(size = 7),  
    axis.title.x = element_text(size = 7, face="bold"),
    axis.title.y = element_text(size = 7, face="bold"),
    strip.text.x = element_text(size = 7),
    strip.text.y= element_text(size = 7),
    # strip.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_blank(),
    legend.text=element_text(size=9),
    legend.title=element_text(size=9, face="bold"),
    plot.title = element_text(size=9, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line = element_line(color = "black")
  )
  return(p_signal_distance)
}

####################################################################################
#
#
#			LD
#
#
#
####################################################################################
Signal_LD_Fun<-function(Signal_ld_data){
  p_size_signal_box = ggplot(data=Signal_ld_data,aes(x = Method, y =abs(Cor),fill=Method))+geom_boxplot( aes(fill=Method))+coord_cartesian(ylim = c(0.5, 1))+scale_fill_manual(values=c("MESuSiE"="#8da0cb","SuSiE"="#66c2a5","SuSiE_weighted"="#B2D3A4","SuSiE_merged"="#9FBA95","Paintor"="#fc8d62","MultiSuSiE"="#e78ac3","SuSiEx"="#E89DA0","XMAP"="#ffd92f","CARMA-X"="#f2b56e","MFD" = "#377eb8","MFD_MESuSiE"="#5A7EC0","MFD_SuSiEx"="#9467BD","MFD_XMAP"="#2D5E9E"),guide=FALSE)+ylab("Correlation with true signal")+xlab("")+facet_grid(~h2,labeller=label_parsed)+theme_bw()
  p_size_signal_box = p_size_signal_box +   theme(
    axis.text.x = element_text(size = 5,face="bold"),
    axis.text.y = element_text(size = 7),  
    axis.title.x = element_text(size = 7, face="bold"),
    axis.title.y = element_text(size = 7, face="bold"),
    strip.text.x = element_text(size = 7),
    strip.text.y= element_text(size = 7),
    # strip.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_blank(),
    legend.text=element_text(size=9),
    legend.title=element_text(size=9, face="bold"),
    plot.title = element_text(size=9, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.line = element_line(color = "black")
  )
}
########################################################################################################################
