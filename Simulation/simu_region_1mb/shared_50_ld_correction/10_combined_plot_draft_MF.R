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

custom_theme <- function() {
    theme(
      # 1. Remove all grids
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Configure text sizes and styles
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size =7),
      axis.title.x = element_text(size = 7, face = "bold"),
      axis.title.y = element_text(size = 7, face = "bold"),
      strip.text.x = element_text(size = 7),
      strip.text.y = element_text(size = 7),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 9, face = "bold"),
      plot.title = element_text(size = 9, hjust = 0.5)
    )
}
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
  "CARMA-X"         = "#f2b56e"
)

## ------------------------------------------------------------------
## 2. Bind the three runs together and tag Pattern
##    (replace *_slalom / *_susie / *_nocontrol with the objects)
## ------------------------------------------------------------------

res_dir <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary"

# 12 files: NO CONTROL (this is the existing res_out_name)
res_out_name_nocontrol <- c(
  "flipped_1_rho_1_updated.RData","flipped_1_rho_2_updated.RData","flipped_1_rho_3_updated.RData","flipped_1_rho_4_updated.RData",
  "flipped_2_rho_1_updated.RData","flipped_2_rho_2_updated.RData","flipped_2_rho_3_updated.RData","flipped_2_rho_4_updated.RData"
)

# 12 files: controlled by slalom
res_out_name_slalom <- c(
  "flipped_1_rho_1_slalom_updated.RData","flipped_1_rho_2_slalom_updated.RData","flipped_1_rho_3_slalom_updated.RData","flipped_1_rho_4_slalom_updated.RData",
  "flipped_2_rho_1_slalom_updated.RData","flipped_2_rho_2_slalom_updated.RData","flipped_2_rho_3_slalom_updated.RData","flipped_2_rho_4_slalom_updated.RData"
)

# 12 files: controlled by susie_rss
res_out_name_susie_rss <- c(
  "flipped_1_rho_1_susie_rss_updated.RData","flipped_1_rho_2_susie_rss_updated.RData","flipped_1_rho_3_susie_rss_updated.RData","flipped_1_rho_4_susie_rss_updated.RData",
  "flipped_2_rho_1_susie_rss_updated.RData","flipped_2_rho_2_susie_rss_updated.RData","flipped_2_rho_3_susie_rss_updated.RData","flipped_2_rho_4_susie_rss_updated.RData"
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
              "2" = "0.01")

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
    rho  = factor(rho),
    Method = recode(Method, "CARMAX" = "CARMA-X") 
  )

# Apply changes to the second dataframe
set_power_all <- set_power_all %>%
  mutate(
    flip = factor(flip),
    rho  = factor(rho),
    Method = recode(Method, "CARMAX" = "CARMA-X") 
  )

## ------------------------------------------------------------------
## 4. Patterned plot functions (Set size, Power, ROC, FDR/Power)
## ------------------------------------------------------------------

Set_Size_fun <- function(all_Set_data_dataframe, upper_limit, pattern_input_Power) {
  
  method_order <- c("MESuSiE", "SuSiE", "SuSiE_weighted", "SuSiE_merged", 
                    "Paintor", "MultiSuSiE", "SuSiEx", "XMAP", "CARMA-X","MFD")
  
  ggplot(all_Set_data_dataframe,
         aes(x = factor(Method, levels = method_order), # Apply method order
             y = Size,
             fill = Method,
             pattern = Pattern)) +
    geom_boxplot(aes(fill = Method), outlier.size = 0.01) +
    coord_cartesian(ylim = c(0, upper_limit)) +
    # hide fill legend, keep pattern legend
    scale_fill_manual(values = method_cols, guide = "none") +
    geom_boxplot_pattern(
      position = position_dodge(preserve = "single"),
      color    = "black",
      size     = 0.05,
      aes(fill = Method, pattern = Pattern),
      pattern_density          = 0.2,
      pattern_size             = 0.1,
      pattern_key_scale_factor = 0.3,
      pattern_spacing          = 0.05,
      pattern_fill             = "white",
      pattern_color            = "white"
    ) +
    scale_pattern_manual(
      values = pattern_input_Power,
      name   = "LD handling"        # legend title; change as needed
    ) +
    guides(
      pattern = guide_legend(
        override.aes = list(
          fill          = "white",
          colour        = NA,
          pattern_fill  = "#DEDBD2",
          pattern_color = "#d8c99b"
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
      strip.text.x     = element_text(size = 5),
      strip.text.y     = element_text(size = 5),
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
  
  method_order <- c("MESuSiE", "SuSiE", "SuSiE_weighted", "SuSiE_merged", 
                    "Paintor", "MultiSuSiE", "SuSiEx", "XMAP", "CARMA-X","MFD")
  
  p <- ggplot(power_summary,
              aes(x = factor(Method, levels = method_order), # Apply method order
                  y = Power_name,
                  fill = Method,
                  pattern = Pattern)) +
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
    scale_fill_manual(values = method_cols, guide = FALSE) +
    scale_pattern_manual(values = pattern_input_Power) +
    guides(
      pattern = guide_legend(
        override.aes = list(
          fill          = "white",
          colour        = NA,
          pattern_fill  = "#DEDBD2",
          pattern_color = "#d8c99b"
        )
      )
    ) +
    ylab("Power") + xlab("") + ylim(c(0, 1)) + theme_bw()
  
  p <- p + geom_text(aes(label = round(Power_name, 2)),
                     position = position_dodge(width = 1),
                     vjust    = -0.5,
                     size     = 5 * 5 / 14)
  
  p + theme(
    axis.text.x      = element_text(size = 5, face = "bold"),
    axis.text.y      = element_text(size = 5),
    axis.title.x     = element_text(size = 7, face = "bold"),
    axis.title.y     = element_text(size = 7, face = "bold"),
    strip.text.x     = element_text(size = 5),
    strip.text.y     = element_text(size = 5),
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
  round() + 300

p_size_box_all <- Set_Size_fun(
  all_Set_data_all %>% mutate(Size = log2(Size + 1)),
  upper_limit = log2(upper_limit_all),
  pattern_input_Power = pattern_input_Power
) +
  ylab("log2(Set Size + 1)") +
  facet_grid(flip ~ rho, labeller = facet_labeller)

p_power_bar_all <- Set_Power_fun(
  set_power_all,
  pattern_input_Power = pattern_input_Power
) +
  facet_grid(flip ~ rho, labeller = facet_labeller)

####################################################################################
#
#
#			Combine LD mismatch figures
#
#
####################################################################################
# --- Step 1: Force Plot A to show NOTHING ---
p_a = p_size_box_all+custom_theme()
p_b = p_power_bar_all+custom_theme()


p_a <- p_a + 
  guides(fill = "none", color = "none", pattern = "none") + 
  custom_theme()

# --- Step 2: Force Plot B to show CLEAN Method Legend & Custom Pattern Legend ---
p_b <- p_b +
  guides(
    # Remove pattern from the Method legend boxes
    fill = guide_legend(
      nrow = 1, 
      order = 2,
      override.aes = list(pattern = "none") # Remove pattern strips from legend keys
    ),
    
    # Remove pattern from the Color legend (if applicable)
    color = guide_legend(
      nrow = 1, 
      order = 2,
      override.aes = list(pattern = "none") 
    ),
    
    # Pattern Legend: Keep the custom stripes here
    pattern = guide_legend(
      nrow = 1, 
      order = 1,
      override.aes = list(
        fill          = "white",
        colour        = "black",
        pattern_fill  = "#DEDBD2",
        pattern_color = "#d8c99b"
      )
    ) 
  ) +
  custom_theme()

# --- Step 3: Combine ---
p_out <- (p_a / p_b) +
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box.margin = margin(t = 10),
    legend.spacing.x = unit(1.0, "cm")
  )

p_out


plot_dir<-paste0(res_dir,"Figure/Draft_main_figure/")

ggsave(
  file.path(plot_dir, "LD_control_size_power_flips_rhos_combined_updated.pdf"),
  p_out,
  height = 250, width = 350, units = "mm", dpi = 600)
