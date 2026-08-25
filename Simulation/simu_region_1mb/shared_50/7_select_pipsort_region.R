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
library(tidyverse)

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

## Option A: capture directly from the shell (run in the folder with those files)
x <- system("wc -l /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_2/summary_data/CAUSAL_2_LOCI_*_h2_1", intern = TRUE)

## Option B: when the lines are already stored in R as a character vector `x`, skip the system() call.

## Parse: keep only file lines, then split into columns
df <- read.table(text = x[grepl("_h2_1$", x)],
                 col.names = c("snps", "file"),
                 colClasses = c("integer", "character"),
                 check.names = FALSE)

## Extract the locus index (the number between ..._LOCI_ and _h2_1)
df$locus <- as.integer(sub(".*CAUSAL_2_LOCI_([0-9]+)_h2_1$", "\\1", df$file))

## (Optional) sort by locus
df <- df[order(df$locus), ]

## Peek
df = df %>% arrange(snps)




nbins <- 10 
df$bin <- cut(df$snps, breaks = nbins, include.lowest = TRUE)

library(dplyr)
df_10_alt <- df %>%
  group_by(bin) %>%
  slice(which.min(abs(snps - median(snps)))) %>%
  ungroup() %>%
  arrange(snps)

nrow(df_10_alt)  # may be < 10 if some bins are empty


df_10_alt$snps


df_10_alt$locus


read_first_number <- function(p) {
  if (!file.exists(p)) return(NA_real_)
  ln <- readr::read_lines(p, n_max = 1)
  readr::parse_number(ln) %>% as.numeric()
}

# Modified to include BOTH h2_1 and h2_2 for each row
make_time_cols <- function(df, root_path, n) {
  # tag in df is assumed to be like: "CAUSAL_2_LOCI_37_h2_1"
  df %>%
    mutate(
      rowid = row_number(),
      # Build the tag for the target n value
      tag_n = str_replace(basename(file), "^CAUSAL_\\d+", paste0("CAUSAL_", n))
    ) %>%
    select(rowid, tag_n) %>%
    # Create both h2 targets (1 and 2) per row
    tidyr::crossing(h2_target = c("1", "2")) %>%
    mutate(
      # Replace the h2 part with target
      tag_n_h2 = str_replace(tag_n, "h2_\\d+", paste0("h2_", h2_target)),
      time_path = file.path(root_path, "result/pip_sort_c_1",
                            paste0("time_", tag_n_h2, "_pipsort.txt")),
      time_val  = purrr::map_dbl(time_path, read_first_number),
      # Column name: time_n_x_h_y_c_1, where x=n and y in {1,2}
      colname   = paste0("time_n_", n, "_h_", h2_target, "_c_1")
    ) %>%
    select(rowid, colname, time_val) %>%
    tidyr::pivot_wider(names_from = colname, values_from = time_val)
}

# Roots
root1 <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_1"
root2 <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_2"
root3 <- "/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_3"

# Start with df_10_alt and a stable row id
df_10_alt <- df_10_alt %>%
  mutate(rowid = row_number())

# Build wide time columns for each causal number (now includes h2_1 and h2_2)
wide1 <- make_time_cols(df_10_alt, root1, n = 1)
wide2 <- make_time_cols(df_10_alt, root2, n = 2)
wide3 <- make_time_cols(df_10_alt, root3, n = 3)


df_10_alt <- df_10_alt %>%
  left_join(wide1, by = "rowid") %>%
  left_join(wide2, by = "rowid") %>%
  left_join(wide3, by = "rowid") %>%
  filter(if_all(starts_with("time_n_"), ~ !is.na(.))) %>%  # all time_* not NA
  select(-rowid) %>% mutate(Method = "PIPSORT")


p_time_pipsort<-ggplot(df_10_alt[1:5,], aes(x = snps, y = time_n_1_h_1_c_1, color = Method)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = c("PIPSORT" = "#f7b801")) +
  labs(x = "Number of SNP", y = "Time (minutes)", title = "Runtime of PIPSORT when c=1") +
  theme_bw()+ custom_theme()
p_time_pipsort


res_dir <- paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/res_summary/")
plot_dir_name <- "shared_50_baseline_updated_meta_updated_xmap"
plot_dir <- paste0(res_dir, "Figure/", plot_dir_name, "/")

ggsave(
  paste0(plot_dir, "Runtime_PIPSORT.pdf"),
  p_time_pipsort,
  height = 100, width = 210, units = "mm", dpi = 600
)




