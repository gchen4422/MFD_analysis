library(dplyr)
library(forcats)
library(ggplot2)
library(ggpattern)
library(patchwork)
library(ggpubr)

## ------------------------------------------------------------------
## 1. Pattern labels for the three LD-control strategies
## ------------------------------------------------------------------

pattern_name_1 <- "Controlled by slalom"
pattern_name_2 <- "Controlled by susie rss"
pattern_name_3 <- "No control"

pattern_input_Power <- c("stripe","circle","none")
names(pattern_input_Power) <- c(pattern_name_1, pattern_name_2, pattern_name_3)

pattern_input_line <- c("dashed","dotted","solid")
names(pattern_input_line) <- c(pattern_name_1, pattern_name_2, pattern_name_3)

## colour palette for methods (same as before)
method_cols <- c(
  "MESuSiE"        = "#8da0cb",
  "SuSiE"          = "#66c2a5",
  "SuSiE_weighted" = "#B2D3A4",
  "SuSiE_merged"   = "#9FBA95",
  "Paintor"        = "#fc8d62",
  "MultiSuSiE"     = "#e78ac3",
  "SuSiEx"         = "#E89DA0",
  "XMAP"           = "#ffd92f",
  "CARMAX"         = "#f2b56e"
)

## ------------------------------------------------------------------
## 2. Bind the three runs together and tag Pattern
##    (replace *_slalom / *_susie / *_nocontrol with the objects)
## ------------------------------------------------------------------

res_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary"

# 12 files: NO CONTROL (this is the existing res_out_name)
res_out_name_nocontrol <- c(
  "flipped_1_rho_1_updated.RData","flipped_1_rho_2_updated.RData","flipped_1_rho_3_updated.RData","flipped_1_rho_4_updated.RData",
  "flipped_2_rho_1_updated.RData","flipped_2_rho_2_updated.RData","flipped_2_rho_3_updated.RData","flipped_2_rho_4_updated.RData",
  "flipped_3_rho_1_updated.RData","flipped_3_rho_2_updated.RData","flipped_3_rho_3_updated.RData","flipped_3_rho_4_updated.RData"
)

# 12 files: controlled by slalom
res_out_name_slalom <- c(
  "flipped_1_rho_1_slalom_updated.RData","flipped_1_rho_2_slalom_updated.RData","flipped_1_rho_3_slalom_updated.RData","flipped_1_rho_4_slalom_updated.RData",
  "flipped_2_rho_1_slalom_updated.RData","flipped_2_rho_2_slalom_updated.RData","flipped_2_rho_3_slalom_updated.RData","flipped_2_rho_4_slalom_updated.RData",
  "flipped_3_rho_1_slalom_updated.RData","flipped_3_rho_2_slalom_updated.RData","flipped_3_rho_3_slalom_updated.RData","flipped_3_rho_4_slalom_updated.RData"
)

# 12 files: controlled by susie_rss
res_out_name_susie_rss <- c(
  "flipped_1_rho_1_susie_rss_updated.RData","flipped_1_rho_2_susie_rss_updated.RData","flipped_1_rho_3_susie_rss_updated.RData","flipped_1_rho_4_susie_rss_updated.RData",
  "flipped_2_rho_1_susie_rss_updated.RData","flipped_2_rho_2_susie_rss_updated.RData","flipped_2_rho_3_susie_rss_updated.RData","flipped_2_rho_4_susie_rss_updated.RData",
  "flipped_3_rho_1_susie_rss_updated.RData","flipped_3_rho_2_susie_rss_updated.RData","flipped_3_rho_3_susie_rss_updated.RData","flipped_3_rho_4_susie_rss_updated.RData"
)

files_nocontrol <- file.path(res_dir, res_out_name_nocontrol)
files_slalom    <- file.path(res_dir, res_out_name_slalom)
files_susie_rss <- file.path(res_dir, res_out_name_susie_rss)


library(dplyr)

load_ld_file <- function(path) {
  env <- new.env()
  load(path, envir = env)
  
  # parse flip and rho from filename, e.g.
  # "flipped_2_rho_3.RData"
  # "flipped_2_rho_3_slalom.RData"
  # "flipped_2_rho_3_susie_rss.RData"
  fname <- basename(path)
  m <- str_match(fname, "flipped_(\\d+)_rho_(\\d+)")
  flip_val <- m[2]
  rho_val  <- m[3]
  
  add_fr <- function(df) {
    df %>%
      mutate(
        flip = flip_val,
        rho  = rho_val
      )
  }
  
  list(
    all_Set_data       = add_fr(env$all_Set_data_dataframe),
    set_power          = add_fr(env$set_power_summary),
    either_ROC         = add_fr(env$either_all_ROC_data_dataframe),
    FDR_Power_either   = add_fr(env$FDR_Power_either),
    shared_ROC         = add_fr(env$shared_all_ROC_data_dataframe),
    FDR_Power_shared   = add_fr(env$FDR_Power_shared),
    ancestry_ROC       = add_fr(env$ancestry_all_ROC_data_dataframe),
    FDR_Power_ancestry = add_fr(env$FDR_Power_ancestry)
  )
}


build_ld_run <- function(files) {
  all_Set_list       <- list()
  set_power_list     <- list()
  either_ROC_list    <- list()
  FDR_either_list    <- list()
  shared_ROC_list    <- list()
  FDR_shared_list    <- list()
  ancestry_ROC_list  <- list()
  FDR_ancestry_list  <- list()
  
  for (i in seq_along(files)) {
    x <- load_ld_file(files[i])
    
    all_Set_list[[i]]      <- x$all_Set_data
    set_power_list[[i]]    <- x$set_power
    either_ROC_list[[i]]   <- x$either_ROC
    FDR_either_list[[i]]   <- x$FDR_Power_either
    shared_ROC_list[[i]]   <- x$shared_ROC
    FDR_shared_list[[i]]   <- x$FDR_Power_shared
    ancestry_ROC_list[[i]] <- x$ancestry_ROC
    FDR_ancestry_list[[i]] <- x$FDR_Power_ancestry
  }
  
  list(
    all_Set_data_all       = bind_rows(all_Set_list),
    set_power_all          = bind_rows(set_power_list),
    either_ROC_all         = bind_rows(either_ROC_list),
    FDR_Power_either_all   = bind_rows(FDR_either_list),
    shared_ROC_all         = bind_rows(shared_ROC_list),
    FDR_Power_shared_all   = bind_rows(FDR_shared_list),
    ancestry_ROC_all       = bind_rows(ancestry_ROC_list),
    FDR_Power_ancestry_all = bind_rows(FDR_ancestry_list)
  )
}

run_nocontrol <- build_ld_run(files_nocontrol)
run_slalom    <- build_ld_run(files_slalom)
run_susie     <- build_ld_run(files_susie_rss)

all_Set_data_all_nocontrol       <- run_nocontrol$all_Set_data_all
set_power_all_nocontrol          <- run_nocontrol$set_power_all
either_ROC_all_nocontrol         <- run_nocontrol$either_ROC_all
FDR_Power_either_all_nocontrol   <- run_nocontrol$FDR_Power_either_all
shared_ROC_all_nocontrol         <- run_nocontrol$shared_ROC_all
FDR_Power_shared_all_nocontrol   <- run_nocontrol$FDR_Power_shared_all
ancestry_ROC_all_nocontrol       <- run_nocontrol$ancestry_ROC_all
FDR_Power_ancestry_all_nocontrol <- run_nocontrol$FDR_Power_ancestry_all

all_Set_data_all_slalom       <- run_slalom$all_Set_data_all
set_power_all_slalom          <- run_slalom$set_power_all
either_ROC_all_slalom         <- run_slalom$either_ROC_all
FDR_Power_either_all_slalom   <- run_slalom$FDR_Power_either_all
shared_ROC_all_slalom         <- run_slalom$shared_ROC_all
FDR_Power_shared_all_slalom   <- run_slalom$FDR_Power_shared_all
ancestry_ROC_all_slalom       <- run_slalom$ancestry_ROC_all
FDR_Power_ancestry_all_slalom <- run_slalom$FDR_Power_ancestry_all

all_Set_data_all_susie       <- run_susie$all_Set_data_all
set_power_all_susie          <- run_susie$set_power_all
either_ROC_all_susie         <- run_susie$either_ROC_all
FDR_Power_either_all_susie   <- run_susie$FDR_Power_either_all
shared_ROC_all_susie         <- run_susie$shared_ROC_all
FDR_Power_shared_all_susie   <- run_susie$FDR_Power_shared_all
ancestry_ROC_all_susie       <- run_susie$ancestry_ROC_all
FDR_Power_ancestry_all_susie <- run_susie$FDR_Power_ancestry_all


pattern_name_1 <- "Controlled by slalom"
pattern_name_2 <- "Controlled by susie rss"
pattern_name_3 <- "No control"

bind_with_pattern <- function(df_sla, df_sus, df_no) {
  bind_rows(
    !!pattern_name_1 := df_sla,
    !!pattern_name_2 := df_sus,
    !!pattern_name_3 := df_no,
    .id = "Pattern"
  )
}

all_Set_data_all       <- bind_with_pattern(all_Set_data_all_slalom,
                                            all_Set_data_all_susie,
                                            all_Set_data_all_nocontrol)
set_power_all          <- bind_with_pattern(set_power_all_slalom,
                                            set_power_all_susie,
                                            set_power_all_nocontrol)
either_ROC_all         <- bind_with_pattern(either_ROC_all_slalom,
                                            either_ROC_all_susie,
                                            either_ROC_all_nocontrol)
FDR_Power_either_all   <- bind_with_pattern(FDR_Power_either_all_slalom,
                                            FDR_Power_either_all_susie,
                                            FDR_Power_either_all_nocontrol)
shared_ROC_all         <- bind_with_pattern(shared_ROC_all_slalom,
                                            shared_ROC_all_susie,
                                            shared_ROC_all_nocontrol)
FDR_Power_shared_all   <- bind_with_pattern(FDR_Power_shared_all_slalom,
                                            FDR_Power_shared_all_susie,
                                            FDR_Power_shared_all_nocontrol)
ancestry_ROC_all       <- bind_with_pattern(ancestry_ROC_all_slalom,
                                            ancestry_ROC_all_susie,
                                            ancestry_ROC_all_nocontrol)
FDR_Power_ancestry_all <- bind_with_pattern(FDR_Power_ancestry_all_slalom,
                                            FDR_Power_ancestry_all_susie,
                                            FDR_Power_ancestry_all_nocontrol)

# make Pattern an ordered factor
for (nm in c("all_Set_data_all","set_power_all",
             "either_ROC_all","FDR_Power_either_all",
             "shared_ROC_all","FDR_Power_shared_all",
             "ancestry_ROC_all","FDR_Power_ancestry_all")) {
  df <- get(nm)
  df$Pattern <- factor(df$Pattern,
                       levels = c(pattern_name_1, pattern_name_2, pattern_name_3))
  assign(nm, df)
}


## ------------------------------------------------------------------
## 3. Facet labeller (same as the LD workflow)
## ------------------------------------------------------------------
flip_map <- c("1" = "0.001",
              "2" = "0.01",
              "3" = "0.1")

rho_map  <- c("1" = "0.5",
              "2" = "0",
              "3" = "-0.5",
              "4" = "-1")

facet_labeller <- labeller(
  flip = function(x) paste0("Flip prop. = ", flip_map[as.character(x)]),
  rho  = function(x) paste0("Rho sig. = ",  rho_map[as.character(x)])
)

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






FDR_Power_either_all <- FDR_Power_either_all %>%
  mutate(
    # Recode 'flip': Match the level "1" to "Flip prop. = 0.001", etc.
    flip = factor(as.character(flip), 
                  levels = names(flip_map), 
                  labels = paste0("Flip prop. = ", flip_map)),
    
    # Recode 'rho': Match the level "1" to "Rho sig. = 0.5", etc.
    rho  = factor(as.character(rho), 
                  levels = names(rho_map), 
                  labels = paste0("Rho sig. = ", rho_map))
  )


FDR_Power_shared_all <- FDR_Power_shared_all %>%
  mutate(
    # Recode 'flip': Match the level "1" to "Flip prop. = 0.001", etc.
    flip = factor(as.character(flip), 
                  levels = names(flip_map), 
                  labels = paste0("Flip prop. = ", flip_map)),
    
    # Recode 'rho': Match the level "1" to "Rho sig. = 0.5", etc.
    rho  = factor(as.character(rho), 
                  levels = names(rho_map), 
                  labels = paste0("Rho sig. = ", rho_map))
  )


FDR_Power_ancestry_all <- FDR_Power_ancestry_all %>%
  mutate(
    # Recode 'flip': Match the level "1" to "Flip prop. = 0.001", etc.
    flip = factor(as.character(flip), 
                  levels = names(flip_map), 
                  labels = paste0("Flip prop. = ", flip_map)),
    
    # Recode 'rho': Match the level "1" to "Rho sig. = 0.5", etc.
    rho  = factor(as.character(rho), 
                  levels = names(rho_map), 
                  labels = paste0("Rho sig. = ", rho_map))
  )



library(openxlsx)
library(stringr)

res_dir_table <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/supp_tables/"
res_out_xlsx <- "LD_outliers_compare.xlsx"
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
  "Either"   = FDR_Power_either_all,
  "Shared"   = FDR_Power_shared_all,
  "Ancestry" = FDR_Power_ancestry_all,
  .id = "Signal_Type" 
) %>%
  select(Method, Signal_Type, Pattern,flip,rho, FDR, Power) %>%
  
  # Apply Factors for ALL sorting columns
  mutate(
    Method = factor(Method, levels = master_method_levels),
    #Pattern = factor(Pattern, levels = pattern_order),
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


## ------------------------------------------------------------------
## 4. Patterned plot functions (Set size, Power, ROC, FDR/Power)
## ------------------------------------------------------------------

Set_Size_fun <- function(all_Set_data_dataframe, upper_limit, pattern_input_Power) {
  
  # Ensure method_cols is available or pass it as an argument
  # method_cols <- ... (defined externally)
  
  ggplot(all_Set_data_dataframe,
         aes(x = reorder(Method, Size, median),
             y = Size,
             fill = Method,
             pattern = Pattern)) +
    
    # Combined Geometry (Handles fill, pattern, and outliers)
    geom_boxplot_pattern(
      position = position_dodge(preserve = "single"),
      color    = "black",
      size     = 0.05,
      outlier.size = 0.01, # Moved outlier size here
      pattern_density          = 0.2,
      pattern_size             = 0.1,
      pattern_key_scale_factor = 0.3,
      pattern_spacing          = 0.05,
      pattern_fill             = "white",
      pattern_color            = "white"
    ) +
    
    coord_cartesian(ylim = c(0, upper_limit)) +
    
    # 1. Enable the fill legend
    scale_fill_manual(values = method_cols) + 
    
    scale_pattern_manual(
      values = pattern_input_Power,
      name   = "LD handling" 
    ) +
    
    # 2. Configure Both Legends
    guides(
      # Method Legend (Fill)
      fill = guide_legend(
        title = "Method",
        order = 1,              # First legend
        nrow = 2,               # Force to 2 rows
        override.aes = list(pattern = "none") # Remove texture from color squares
      ),
      
      # LD Handling Legend (Pattern)
      pattern = guide_legend(
        title = "LD handling",
        order = 2,              # Second legend
        nrow = 2,               # Align height with Method legend
        override.aes = list(
          fill          = "white",
          colour        = "black", # changed to black for visibility in legend
          pattern_fill  = "black", # changed for contrast in legend
          pattern_color = "black"
        )
      )
    ) +
    
    ylab("Set Size") +
    xlab("") +
    theme_bw() +
    theme(
      axis.text.x      = element_text(size = 5, face = "bold"),
      axis.text.y      = element_text(size = 5),
      axis.title.x     = element_text(size = 7, face = "bold"),
      axis.title.y     = element_text(size = 7, face = "bold"),
      strip.text.x     = element_text(size = 7),
      strip.text.y     = element_text(size = 7),
      strip.background = element_blank(),
      legend.text      = element_text(size = 7),
      legend.title     = element_text(size = 7, face = "bold"),
      legend.position  = "bottom",
      plot.title       = element_text(size = 7, hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),
      axis.line        = element_line(color = "black")
    ) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = -Inf, ymax = -Inf,
             colour = "black", size = 0.2)
}


Set_Power_fun <- function(power_summary, pattern_input_Power) {
  
  # Ensure method_cols is available (or pass it as an argument)
  
  p <- ggplot(power_summary,
              aes(x = Method,
                  y = Power_name,
                  fill = Method,
                  pattern = Pattern)) +
    
    # 1. Geometry
    geom_col_pattern(
      position = "dodge",
      color    = "#212529",
      size     = 0.05,
      aes(fill = Method, pattern = Pattern),
      pattern_density          = 0.2,
      pattern_size             = 0.1,
      pattern_key_scale_factor = 0.3,
      pattern_spacing          = 0.05,
      pattern_fill             = "white",
      pattern_color            = "white"
    ) +
    
    # 2. Scales without 'guide = FALSE'
    scale_fill_manual(values = method_cols) + 
    scale_pattern_manual(values = pattern_input_Power) +
    
    # 3. Configure legends with fill
    guides(
      # Method Legend
      fill = guide_legend(
        title = "Method",
        order = 1,              # First legend
        nrow = 2,               # Force to 2 rows
        override.aes = list(pattern = "none") # Remove texture from color keys
      ),
      
      # LD Handling Legend
      pattern = guide_legend(
        title = "LD handling",
        order = 2,              # Second legend
        nrow = 2,               # Match height
        override.aes = list(
          fill          = "white",
          colour        = "black", # Dark border for visibility
          pattern_fill  = "black", # Dark pattern for visibility
          pattern_color = "black"
        )
      )
    ) +
    
    # 4. Facets and labels
    ylab("Power") + xlab("") + 
    ylim(c(0, 1)) + 
    facet_grid(h2 ~ causal_num, labeller = label_parsed) +
    theme_bw()
  
  # 5. Add Text Labels
  #p <- p + geom_text(aes(label = round(Power_name, 2)),
  #                   position = position_dodge(width = 1),
  #                   vjust    = -0.5,
  #                   size     = 5 * 5 / 14)
  
  # 6. Apply Theme
  p + theme(
    axis.text.x      = element_text(size = 5, face = "bold"),
    axis.text.y      = element_text(size = 5),
    axis.title.x     = element_text(size = 7, face = "bold"),
    axis.title.y     = element_text(size = 7, face = "bold"),
    strip.text.x     = element_text(size = 7),
    strip.text.y     = element_text(size = 7),
    strip.background = element_blank(),
    
    # Legend settings
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 7, face = "bold"),
    legend.position  = "bottom",
    
    plot.title       = element_text(size = 7, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black")
  ) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = -Inf, ymax = -Inf,
             colour = "black", size = 0.2)
}

ROC_shared_fun <- function(df, pattern_input_line) {
  ggplot(df, aes(x = Power, y = 1 - FDR)) +
    geom_line(aes(linetype = Pattern, color = Method), size = 0.8) +
    scale_color_manual(values = method_cols) +
    scale_linetype_manual(values = pattern_input_line) +
    xlab("Recall") +
    ylab("Precision") +
    theme_bw() +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
    theme(
      axis.text.x      = element_text(size = 5, face = "bold"),
      axis.text.y      = element_text(size = 5),
      axis.title.x     = element_text(size = 7, face = "bold"),
      axis.title.y     = element_text(size = 7, face = "bold"),
      strip.text.x     = element_text(size = 7),
      strip.text.y     = element_text(size = 7),
      strip.background = element_blank(),
      legend.text      = element_text(size = 7),
      legend.title     = element_text(size = 7, face = "bold"),
      plot.title       = element_text(size = 7, hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),
      axis.line        = element_line(color = "black"),
      legend.position  = "bottom"
    ) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = -Inf, ymax = -Inf,
             colour = "black", size = 0.2)
}

FDR_Power_shared_fun <- function(df, pattern_input_Power) {
  
  ggplot(df, aes(FDR, Power, fill = Method, pattern = Pattern)) +
    
    # 1. Geometry
    geom_col_pattern(
      position = "dodge",
      color    = "#212529",
      size     = 0.05,
      aes(fill = Method, pattern = Pattern),
      pattern_density          = 0.2,
      pattern_size             = 0.1,
      pattern_key_scale_factor = 0.3,
      pattern_spacing          = 0.05,
      pattern_fill             = "white",
      pattern_color            = "white"
    ) +
    
    # 2. Scales without 'guide = FALSE'
    scale_fill_manual(values = method_cols) + 
    scale_pattern_manual(values = pattern_input_Power) +
    
    # 3. Configure legends with fill and pattern
    guides(
      # Method Legend
      fill = guide_legend(
        title = "Method",
        order = 1,              # First legend
        nrow = 2,               # Force to 2 rows
        override.aes = list(pattern = "none") # Solid colors (no stripes on method keys)
      ),
      
      # LD Handling Legend
      pattern = guide_legend(
        title = "LD handling",
        order = 2,              # Second legend
        nrow = 2,               # Match height
        override.aes = list(
          fill          = "white",
          colour        = "black", # Dark border so box is visible
          pattern_fill  = "black", # Dark pattern so stripes are visible
          pattern_color = "black"
        )
      )
    ) +
    
    # 4. Facets and Labels
    ylab("Power") + xlab("FDR") +
    facet_grid(h2 ~ causal_num, labeller = label_parsed) + # Standard faceting
    theme_bw() +
    
    # 5. Theme
    theme(
      legend.position  = "bottom",
      axis.text.x      = element_text(size = 5, face = "bold"),
      axis.text.y      = element_text(size = 5),
      axis.title.x     = element_text(size = 7, face = "bold"),
      axis.title.y     = element_text(size = 7, face = "bold"),
      strip.text.x     = element_text(size = 7),
      strip.text.y     = element_text(size = 7),
      strip.background = element_blank(),
      legend.text      = element_text(size = 7),
      legend.title     = element_text(size = 7, face = "bold"),
      plot.title       = element_text(size = 7, hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),
      axis.line        = element_line(color = "black")
    ) +
    
    # 6. Border
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = -Inf, ymax = -Inf,
             colour = "black", size = 0.2)
}

## Ancestry versions: same idea but the existing splitting of Method/Ancestry
ROC_ancestry_fun <- function(df, pattern_input_line) {
  ggplot(df, aes(x = Power, y = 1 - FDR)) +
    geom_line(aes(linetype = Pattern, color = Method), size = 0.8) +
    scale_color_manual(values = method_cols) +
    scale_linetype_manual(values = pattern_input_line) +
    xlab("Recall") +
    ylab("Precision") +
    theme_bw() +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
    theme(
      axis.text.x      = element_text(size = 5, face = "bold"),
      axis.text.y      = element_text(size = 5),
      axis.title.x     = element_text(size = 7, face = "bold"),
      axis.title.y     = element_text(size = 7, face = "bold"),
      strip.text.x     = element_text(size = 7),
      strip.text.y     = element_text(size = 7),
      strip.background = element_blank(),
      legend.text      = element_text(size = 7),
      legend.title     = element_text(size = 7, face = "bold"),
      plot.title       = element_text(size = 7, hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),
      axis.line        = element_line(color = "black"),
      legend.position  = "bottom"
    ) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = -Inf, ymax = -Inf,
             colour = "black", size = 0.2)
}


FDR_Power_ancestry_fun <- function(df, pattern_input_Power) {
  
  ggplot(df, aes(FDR, Power, fill = Method, pattern = Pattern)) +
    
    # 1. Geometry
    geom_col_pattern(
      position = "dodge",
      color    = "#212529",
      size     = 0.05,
      aes(fill = Method, pattern = Pattern),
      pattern_density          = 0.2,
      pattern_size             = 0.1,
      pattern_key_scale_factor = 0.3,
      pattern_spacing          = 0.05,
      pattern_fill             = "white",
      pattern_color            = "white"
    ) +
    
    # 2. Scales without 'guide = FALSE'
    scale_fill_manual(values = method_cols) + 
    scale_pattern_manual(values = pattern_input_Power) +
    
    # 3. Configure legends with fill
    guides(
      # Method Legend
      fill = guide_legend(
        title = "Method",
        order = 1,              # First legend
        nrow = 2,               # Force to 2 rows
        override.aes = list(pattern = "none") # CRITICAL: Makes color keys solid (no stripes)
      ),
      
      # Ancestry/Pattern Legend
      pattern = guide_legend(
        title = "Ancestry",     # Adjusted title based on function name
        order = 2,              # Second legend
        nrow = 2,               # Match height
        override.aes = list(
          fill          = "white",
          colour        = "black", # Dark border for visibility
          pattern_fill  = "black", # Dark pattern for visibility
          pattern_color = "black"
        )
      )
    ) +
    
    # 4. Facets and Labels
    ylab("Power") + xlab("FDR") +
    facet_grid(h2 ~ causal_num, labeller = label_parsed) + # Included based on context
    theme_bw() +
    
    # 5. Theme
    theme(
      legend.position  = "bottom",
      axis.text.x      = element_text(size = 5, face = "bold"),
      axis.text.y      = element_text(size = 5),
      axis.title.x     = element_text(size = 7, face = "bold"),
      axis.title.y     = element_text(size = 7, face = "bold"),
      strip.text.x     = element_text(size = 7),
      strip.text.y     = element_text(size = 7),
      strip.background = element_blank(),
      legend.text      = element_text(size = 7),
      legend.title     = element_text(size = 7, face = "bold"),
      plot.title       = element_text(size = 7, hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),
      axis.line        = element_line(color = "black")
    ) +
    
    # 6. Border
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = -Inf, ymax = -Inf,
             colour = "black", size = 0.2)
}


## ------------------------------------------------------------------
## 5. Reproduce the LD workflow, but with patterns
## ------------------------------------------------------------------

res_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/"

# global upper limit from Paintor baseline
upper_limit_all <- all_Set_data_all %>%
  filter(Method == "Paintor",
         causal_num == "Num~Causal  == 1 ",
         h2 == "~h^2 == 10^-4") %>%
  summarise(upper = quantile(Size, 0.75)) %>%
  pull(upper) %>%
  round() + 200

p_size_box_all <- Set_Size_fun(
  all_Set_data_all %>% mutate(Size = log2(Size + 1)),
  upper_limit = log2(upper_limit_all),
  pattern_input_Power = pattern_input_Power
) +
  ylab("log2(Set Size + 1)") +
  facet_grid(flip ~ rho, labeller = facet_labeller)+
  # Remove X-axis labels/ticks and suppress legend completely
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) 



ggsave(
  file.path(res_dir, "Figure", "LD_control_Set_size_all_flips_rhos_combined.pdf"),
  p_size_box_all,
  height = 180, width = 210, units = "mm", dpi = 600
)

p_power_bar_all <- Set_Power_fun(
  set_power_all,
  pattern_input_Power = pattern_input_Power
) +
  facet_grid(flip ~ rho, labeller = facet_labeller)+
  # Remove X-axis labels/ticks and suppress legend completely
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) 


ggsave(
  file.path(res_dir, "Figure", "LD_control_power_all_flips_rhos_combined.pdf"),
  p_power_bar_all,
  height = 180, width = 210, units = "mm", dpi = 600
)

## -------- Either ancestry --------
either_ROC_all <- either_ROC_all %>%
  mutate(Method = fct_relevel(Method,
                              "MESuSiE","SuSiE","Paintor",
                              "SuSiEx","XMAP","MultiSuSiE","CARMAX"))

p_ROC_Either_all <- ROC_shared_fun(either_ROC_all, pattern_input_line) +
  facet_grid(flip ~ rho, labeller = facet_labeller)

ggsave(
  file.path(res_dir, "Figure", "LD_control_roc_either_all_flips_rhos_combined.pdf"),
  p_ROC_Either_all,
  height = 180, width = 210, units = "mm", dpi = 600
)

power_upper_limit_either_all <- FDR_Power_either_all %>%
  filter(FDR != 0.5) %>%
  ungroup() %>%
  summarise(upper_limit = min(ceiling(max(Power) * 10) / 10 + 0.1, 1)) %>%
  pull(upper_limit)

p_FDR_Power_Either_all <- FDR_Power_shared_fun(
  FDR_Power_either_all %>% filter(FDR != 0.5),
  pattern_input_Power
) +
  ylim(0, power_upper_limit_either_all) +
  facet_grid(flip ~ rho, labeller = facet_labeller)



ggsave(
  file.path(res_dir, "Figure", "LD_control_fdr_power_either_all_flips_rhos_combined.pdf"),
  p_FDR_Power_Either_all,
  height = 180, width = 210, units = "mm", dpi = 600
)

## -------- Shared signal --------
shared_ROC_all <- shared_ROC_all %>%
  mutate(Method = fct_relevel(Method,
                              "MESuSiE","SuSiE","Paintor",
                              "SuSiEx","XMAP","MultiSuSiE","CARMAX"))

p_ROC_shared_all <- ROC_shared_fun(shared_ROC_all, pattern_input_line) +
  facet_grid(flip ~ rho, labeller = facet_labeller)

power_upper_limit_shared_all <- FDR_Power_shared_all %>%
  filter(FDR != 0.5) %>%
  ungroup() %>%
  summarise(upper_limit = min(ceiling(max(Power) * 10) / 10 + 0.1, 1)) %>%
  pull(upper_limit)

p_FDR_Power_shared_all <- FDR_Power_shared_fun(
  FDR_Power_shared_all %>% filter(FDR != 0.5),
  pattern_input_Power
) +
  ylim(0, power_upper_limit_shared_all) +
  facet_grid(flip ~ rho, labeller = facet_labeller)

ggsave(
  file.path(res_dir, "Figure", "LD_control_roc_shared_all_flips_rhos_combined.pdf"),
  p_ROC_shared_all,
  height = 180, width = 210, units = "mm", dpi = 500
)

ggsave(
  file.path(res_dir, "Figure", "LD_control_fdr_power_shared_all_flips_rhos_combined.pdf"),
  p_FDR_Power_shared_all,
  height = 180, width = 210, units = "mm", dpi = 500
)

## -------- Ancestry-specific --------

# split Method into Method + Ancestry as before
split_list_ROC <- strsplit(as.character(ancestry_ROC_all$Method), " +")

ancestry_ROC_all <- ancestry_ROC_all %>%
  ungroup() %>%
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

p_ROC_ancestry_wb <- ROC_ancestry_fun(ancestry_ROC_all %>% filter(Ancestry == "White British"), pattern_input_line) +
  facet_grid(flip ~ rho, labeller = facet_labeller)

p_ROC_ancestry_bb <- ROC_ancestry_fun(ancestry_ROC_all %>% filter(Ancestry == "Black British"), pattern_input_line) +
  facet_grid(flip ~ rho, labeller = facet_labeller)

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

p_FDR_Power_ancestry_wb <- FDR_Power_ancestry_fun(
  FDR_Power_ancestry_all %>% filter(Ancestry == "White British") ,
  pattern_input_Power
) +
  ylim(0, power_upper_limit_ancestry_all) +
  facet_grid(flip ~ rho, labeller = facet_labeller)


p_FDR_Power_ancestry_bb <- FDR_Power_ancestry_fun(
  FDR_Power_ancestry_all %>% filter(Ancestry == "Black British") ,
  pattern_input_Power
) +
  ylim(0, power_upper_limit_ancestry_all) +
  facet_grid(flip ~ rho, labeller = facet_labeller)



ggsave(
  file.path(res_dir, "Figure", "LD_control_roc_WB_all_flips_rhos_combined.pdf"),
  p_ROC_ancestry_wb,
  height = 180, width = 210, units = "mm", dpi = 600
)

ggsave(
  file.path(res_dir, "Figure", "LD_control_roc_BB_all_flips_rhos_combined.pdf"),
  p_ROC_ancestry_bb,
  height = 180, width = 210, units = "mm", dpi = 600
)

ggsave(
  file.path(res_dir, "Figure", "LD_control_fdr_power_WB_all_flips_rhos_combined.pdf"),
  p_FDR_Power_ancestry_wb,
  height = 180, width = 210, units = "mm", dpi = 600
)


ggsave(
  file.path(res_dir, "Figure", "LD_control_fdr_power_BB_all_flips_rhos_combined.pdf"),
  p_FDR_Power_ancestry_bb,
  height = 180, width = 210, units = "mm", dpi = 600
)
