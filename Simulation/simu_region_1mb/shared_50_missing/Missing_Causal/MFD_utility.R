library(dplyr)
library(data.table)
library(tidyr)

# --- Define the Decision Function ---
decide_finemapping_method <- function(sum_eur, sum_afr, ld_eur, ld_afr, p_thresh = 5e-8, r2_thresh = 0.6) {
  
  # 1. Pre-processing
  sum_eur[, SNP := as.character(SNP)]
  sum_afr[, SNP := as.character(SNP)]
  
  check_and_convert_cor <- function(mat) {
    if (any(abs(diag(mat) - 1) > 0.1)) cov2cor(mat) else mat
  }
  R_eur <- check_and_convert_cor(ld_eur)
  R_afr <- check_and_convert_cor(ld_afr)
  
  # 2. Identify Variant Sets
  shared_snps <- intersect(sum_eur$SNP, sum_afr$SNP)
  asv_eur <- setdiff(sum_eur$SNP, sum_afr$SNP)
  asv_afr <- setdiff(sum_afr$SNP, sum_eur$SNP)
  
  # 3. Check Condition 1: Are AS-Vs Significant?
  sig_asv_eur <- sum_eur[SNP %in% asv_eur & PVAL < p_thresh, SNP]
  sig_asv_afr <- sum_afr[SNP %in% asv_afr & PVAL < p_thresh, SNP]
  has_sig_asv <- length(sig_asv_eur) > 0 || length(sig_asv_afr) > 0
  
  if (!has_sig_asv) {
    return(list(method = "MESuSiE", reason_code = 1,reason = "AS-V not significant"))
  }
  
  # 4. Prepare for LD Checks
  sig_shared_eur <- sum_eur[SNP %in% shared_snps & PVAL < p_thresh, SNP]
  sig_shared_afr <- sum_afr[SNP %in% shared_snps & PVAL < p_thresh, SNP]
  all_sig_shared <- unique(c(sig_shared_eur, sig_shared_afr))
  
  if (length(all_sig_shared) == 0) { # This case is that no significant shared signals even present
    return(list(method = "SuSiE post-hoc",reason_code = 2,reason = "AS-V significant, but not in high LD with shared signals"))
  }
  
  # 5. Check High LD Logic
  high_ld_found <- FALSE
  high_ld_partners <- character()
  
  # EUR Check
  if (length(sig_asv_eur) > 0) {
    valid_asv <- intersect(sig_asv_eur, rownames(R_eur))
    valid_shared <- intersect(all_sig_shared, colnames(R_eur))
    if (length(valid_asv) > 0 && length(valid_shared) > 0) {
      ld_sub <- R_eur[valid_asv, valid_shared, drop=FALSE]
      if (max(ld_sub^2, na.rm = TRUE) > r2_thresh) {
        high_ld_found <- TRUE
        high_ld_partners <- c(high_ld_partners, colnames(ld_sub)[apply(ld_sub^2, 2, max) > r2_thresh])
      }
    }
  }
  
  # AFR Check
  if (length(sig_asv_afr) > 0) {
    valid_asv <- intersect(sig_asv_afr, rownames(R_afr))
    valid_shared <- intersect(all_sig_shared, colnames(R_afr))
    if (length(valid_asv) > 0 && length(valid_shared) > 0) {
      ld_sub <- R_afr[valid_asv, valid_shared, drop=FALSE]
      if (max(ld_sub^2, na.rm = TRUE) > r2_thresh) {
        high_ld_found <- TRUE
        high_ld_partners <- c(high_ld_partners, colnames(ld_sub)[apply(ld_sub^2, 2, max) > r2_thresh])
      }
    }
  }
  
  high_ld_partners <- unique(high_ld_partners)
  
  if (!high_ld_found) {
    return(list(method = "SuSiE post-hoc",reason_code = 2,reason = "AS-V significant, but not in high LD with shared signals"))
  }
  
  # 6. Check Condition 3: Nature of Shared SNPs
  p_eur_check <- sum_eur[SNP %in% high_ld_partners, PVAL]
  p_afr_check <- sum_afr[SNP %in% high_ld_partners, PVAL]
  is_sig_both <- (p_eur_check < p_thresh) & (p_afr_check < p_thresh)
  
  if (any(is_sig_both)) {
    return(list(method = "MESuSiE", reason_code = 3,reason = "AS-V significant and in high LD with a Shared SNP significant in ALL ancestries"))
  } else {
    return(list(method = "SuSiE post-hoc",reason_code = 4,reason = "AS-V significant and in high LD with a Shared SNP significant in ONLY specific ancestry"))
  }
}

in_CS <- function (res, coverage = 0.95) {
  
  # Apply the in_CS_x function to each row of the 'res' matrix
  result_matrix = t(apply(res, 1, function(x) in_CS_x(x, coverage)))
  
  return(result_matrix)
}

get_purity <- function(pos, Xcorr) {
  
  # If there's only one position, return a vector of ones
  if (length(pos) == 1) {
    return(c(1, 1, 1))
  } else {
    # Extract absolute correlations for specified positions
    value_list = lapply(Xcorr, function(x) c(abs(x[pos, pos])))
    
    # Combine the lists column-wise
    value_matrix = Reduce(cbind, value_list)
    
    # Compute the maximum for each row
    value_max = do.call(pmax, data.frame(value_matrix))
    
    # Return the desired statistics
    return(c(min(value_max, na.rm = TRUE),
             mean(value_max, na.rm = TRUE),
             median(value_max, na.rm = TRUE)))
  }
}

meSuSie_get_cs_specific <- function(res, Xcorr, target_idx, coverage = 0.95, prior_tol = 1e-9, cor_method = "min.abs.corr", cor_threshold = 0.5){
  # target_idx: 1 = WB-specific, 2 = BB-specific, 3 = Shared (configuration labels)
  
  # 1. Identify effects that exist (total sum > tolerance)
  include_idx = unlist(lapply(res$V, function(x) max(diag(x)) > prior_tol))
  
  # 2. Extract Alpha for the SPECIFIC component and NORMALIZE (Conditioning)
  # Original: alpha_either = t(Reduce(cbind,lapply(res$alpha,function(x)apply(x,1,sum))))
  alpha_specific = t(Reduce(cbind, lapply(res$alpha, function(x) {
    # Extract the specific column (WB, BB, or Shared)
    target_vec <- x[, target_idx]
    
    # Check if the component has enough mass to matter, otherwise return 0s
    if(sum(target_vec) < 1e-10) return(rep(0, length(target_vec)))
    
    # NORMALIZE: "Given this is a Shared effect, what is the CS?"
    # This forces the sum to 1 so in_CS can find a 95% set.
    return(target_vec / sum(target_vec)) 
  })))
  
  # 3. Use internal SuSiE function to find CS based on this specific alpha
  # Note: ensure 'in_CS' is loaded from susieR namespace or defined in the calling script
  status = in_CS(alpha_specific, coverage = coverage)
  
  cs = lapply(1:nrow(status), function(i) which(status[i,] != 0))
  
  # 4. Filter and Purity Checks (Same as original)
  include_idx = include_idx * (lapply(cs, length) > 0)
  include_idx = include_idx * (!(duplicated(cs)))
  include_idx = as.logical(include_idx)
  
  if (sum(include_idx) == 0)
    return(list(cs = NULL, cs_index = NULL))
  
  cs = cs[include_idx]
  
  # Calculate Purity based on Xcorr
  purity = data.frame(do.call(rbind, lapply(1:length(cs), function(i){
    get_purity(cs[[i]], Xcorr) 
  })))
  
  colnames(purity) = c("min.abs.corr", "mean.abs.corr", "median.abs.corr")
  
  # Filter by correlation threshold
  is_pure = which(purity[, colnames(purity) %in% cor_method] >= cor_threshold)
  
  if (length(is_pure) > 0) {
    cs = cs[is_pure]
    cs_index = which(include_idx)[is_pure]
    
    # Assign names (L1, L2, etc)
    row_names = paste0("L", cs_index)
    names(cs) = row_names
    
    return(list(cs = cs, cs_index = cs_index))
  } else {
    return(list(cs = NULL, cs_index = NULL))
  }
}


#' Map Credible Set indices to a full SNP vector
#'
#' @param cs_res The result object from meSuSie_get_cs_specific (containing $cs and $cs_index)
#' @param n_snps The total number of SNPs (usually length(MESuSiE_res$pip))
#' @param renumber Logical. If TRUE (default), renumbers sets 1, 2, 3... 
#'                 If FALSE, keeps original master indices (e.g., 4, 5...)
#'
#' @return A vector of length n_snps containing integer labels
get_cs_index_vector <- function(cs_res, n_snps, renumber = TRUE) {
  
  # 1. Initialize the vector with 0s
  cs_vec <- rep(0, n_snps)
  
  # 2. Check if CS results exist
  if (!is.null(cs_res) && !is.null(cs_res$cs)) {
    
    # Iterate through each Credible Set found
    for (j in seq_along(cs_res$cs)) {
      
      # Determine which ID to assign
      if (renumber) {
        # NEW WAY: Sequential numbering (1, 2, 3...)
        # Useful for individual population plots
        assign_id <- j 
      } else {
        # OLD WAY: Master indices (e.g., 4, 5...)
        # Useful for comparing identity across populations
        assign_id <- cs_res$cs_index[j]
      }
      
      # Extract SNP indices for this set
      snp_indices <- cs_res$cs[[j]]
      
      # Assign the ID to the vector
      cs_vec[snp_indices] <- assign_id
    }
  }
  
  return(cs_vec)
}




#' Generalized Multi-Ancestry Fine-Mapping Function
#' 
#' @param gwas_1 Dataframe for Pop 1 (SNP, CHR, POS, Z, Beta, Se, N)
#' @param gwas_2 Dataframe for Pop 2 (SNP, CHR, POS, Z, Beta, Se, N)
#' @param ld_1 Matrix or file path for Pop 1 LD (matched to gwas_1)
#' @param ld_2 Matrix or file path for Pop 2 LD (matched to gwas_2)
#' @param pop_names Vector of two strings, e.g., c("EUR", "AFR")

run_mf_decision_fm <- function(gwas_1, gwas_2, ld_1, ld_2, pop_names = c("Pop1", "Pop2")) {
  
  require(dplyr)
  require(data.table)
  require(susieR)
  require(MESuSiE)
  
  
  colnames(ld_1)<-gwas_1$SNP
  rownames(ld_1)<-gwas_1$SNP
  
  colnames(ld_2)<-gwas_2$SNP
  rownames(ld_2)<-gwas_2$SNP
  
  
  # 1. Decision Logic (Using the custom decision function)
  # Ensure input is setDT for consistency
  setDT(gwas_1); setDT(gwas_2)
  decision <- decide_finemapping_method(gwas_1, gwas_2, ld_1, ld_2)
  
  # --- CASE 1: SuSiE post-hoc (Independent signals) ---
  if(decision$method == "SuSiE post-hoc") {
    
    susie_1 <- susie_rss(gwas_1$Z, ld_1, check_prior = FALSE)
    susie_2 <- susie_rss(gwas_2$Z, ld_2, check_prior = FALSE)
    
    # Map CS indices
    cs_1 <- get_cs_index_vector(susie_1$sets, nrow(gwas_1), renumber = TRUE)
    cs_2 <- get_cs_index_vector(susie_2$sets, nrow(gwas_2), renumber = TRUE)
    
    res1_df <- data.frame(SNP = gwas_1$SNP, PIP1 = susie_1$pip, CS1 = cs_1)
    res2_df <- data.frame(SNP = gwas_2$SNP, PIP2 = susie_2$pip, CS2 = cs_2)
    
    mf_result <- gwas_1 %>% 
      select(SNP, CHR, POS) %>% 
      full_join(gwas_2 %>% select(SNP, CHR, POS), by = "SNP") %>%
      mutate(CHR = coalesce(CHR.x, CHR.y), POS = coalesce(POS.x, POS.y)) %>%
      select(SNP, CHR, POS) %>%
      left_join(res1_df, by = "SNP") %>%
      left_join(res2_df, by = "SNP") %>%
      mutate(across(everything(), ~replace_na(., 0))) %>%
      mutate(
        PIP_Either = pmax(PIP1, PIP2),
        PIP_Shared = pmin(PIP1, PIP2),
        CS = ifelse(CS1 + CS2 == 0, 0, 1)
      ) %>%
      select(SNP, CHR, POS, PIP_Either, PIP_Shared, PIP1, PIP2, CS, CS1, CS2)
    
  } else {
    
    # --- CASE 2: MESuSiE (Shared/Common variants) ---
    
    # MESuSiE requires variants present in both datasets
    common_snps <- intersect(gwas_1$SNP, gwas_2$SNP)
    
    g1_sub <- gwas_1[SNP %in% common_snps][order(match(SNP, common_snps))]
    g2_sub <- gwas_2[SNP %in% common_snps][order(match(SNP, common_snps))]
    
    # Subset LD matrices based on SNP IDs (assuming LD has row/col names)
    ld1_sub <- ld_1[common_snps, common_snps]
    ld2_sub <- ld_2[common_snps, common_snps]
    
    summary_stat_sd_list = list(g1_sub, g2_sub)
    names(summary_stat_sd_list) <- pop_names
    R_mat_list = list(ld1_sub, ld2_sub)
    names(R_mat_list) <- pop_names
    
    mesusie_res <- meSuSie_core(R_mat_list, summary_stat_sd_list, L = 10)
    
    # Helper to count CS categories dynamically
    get_pop_cs_vec <- function(mesusie_obj, pop_name) {
      vec <- rep(0, length(mesusie_obj$pip))
      # Extract indices where category matches user-provided pop_name
      pop_sets <- mesusie_obj$cs$cs[mesusie_obj$cs$cs_category == pop_name]
      if (length(pop_sets) > 0) {
        for (j in seq_along(pop_sets)) {
          vec[pop_sets[[j]]] <- j
        }
      }
      return(vec)
    }
    
    cs_either <- get_cs_index_vector(mesusie_res$cs, length(mesusie_res$pip), renumber = TRUE)
    cs_p1 <- get_pop_cs_vec(mesusie_res, pop_names[1])
    cs_p2 <- get_pop_cs_vec(mesusie_res, pop_names[2])
    
    mesusie_df <- data.frame(
      SNP = common_snps,
      PIP_Either = mesusie_res$pip,
      PIP1 = mesusie_res$pip_config[,1],
      PIP2 = mesusie_res$pip_config[,2],
      PIP_Shared = mesusie_res$pip_config[,3],
      CS = cs_either,
      CS1 = cs_p1,
      CS2 = cs_p2
    )
    
    mf_result <- gwas_1 %>% 
      select(SNP, CHR, POS) %>% 
      full_join(gwas_2 %>% select(SNP, CHR, POS), by = "SNP") %>%
      mutate(CHR = coalesce(CHR.x, CHR.y), POS = coalesce(POS.x, POS.y)) %>%
      select(SNP, CHR, POS) %>%
      left_join(mesusie_df, by = "SNP") %>%
      mutate(across(everything(), ~replace_na(., 0))) %>%
      select(SNP, CHR, POS, PIP_Either, PIP_Shared, PIP1, PIP2, CS, CS1, CS2)
  }
  
  # 3. Final Column Formatting
  colnames(mf_result) <- c("SNP", "CHR", "POS", "PIP_Either", "PIP_Shared", 
                           "PIP_Ancestry_1", "PIP_Ancestry_2", "CS", 
                           "CS_Ancestry_1", "CS_Ancestry_2")
  
  return(mf_result %>% arrange(CHR, POS))
}





#' Generalized Multi-Ancestry Fine-Mapping Function
#' 
#' @param gwas_1 Dataframe for Pop 1 (SNP, CHR, POS, Z, Beta, Se, N)
#' @param gwas_2 Dataframe for Pop 2 (SNP, CHR, POS, Z, Beta, Se, N)
#' @param ld_1 Matrix or file path for Pop 1 LD (matched to gwas_1)
#' @param ld_2 Matrix or file path for Pop 2 LD (matched to gwas_2)
#' @param pop_names Vector of two strings, e.g., c("EUR", "AFR")


run_mf_decision <- function(gwas_1, gwas_2, ld_1, ld_2, 
                                             pop_names = c("Pop1", "Pop2"),
                                             # --- A. Model Configuration ---
                                             L = 10,
                                             prior_weights = NULL,
                                             ancestry_weight = NULL,
                                             # --- B. Decision Logic Thresholds ---
                                             p_thresh = 5e-8,
                                             r2_thresh = 0.6) {
  
  require(dplyr)
  require(data.table)
  require(susieR)
  require(MESuSiE)
  
  
  colnames(ld_1)<-gwas_1$SNP
  rownames(ld_1)<-gwas_1$SNP
  
  colnames(ld_2)<-gwas_2$SNP
  rownames(ld_2)<-gwas_2$SNP
  
  
  # 1. Decision Logic (Using the custom decision function)
  # Ensure input is setDT for consistency
  setDT(gwas_1); setDT(gwas_2)
  decision <- decide_finemapping_method(gwas_1, gwas_2, ld_1, ld_2, 
                                        p_thresh = p_thresh, 
                                        r2_thresh = r2_thresh)
  raw_objects <- list()
  
  # --- CASE 1: SuSiE post-hoc (Independent signals) ---
  if(decision$method == "SuSiE post-hoc") {
    
    susie_1 <- susie_rss(gwas_1$Z, ld_1, L = L, prior_weights = prior_weights, check_prior = FALSE)
    susie_2 <- susie_rss(gwas_2$Z, ld_2, L = L, prior_weights = prior_weights, check_prior = FALSE)
    
    # Store raw objects
    raw_objects$susie_1 <- susie_1
    raw_objects$susie_2 <- susie_2
    
    # 1. Standardize CS indices locally

    cs_1 <- get_cs_index_vector(susie_1$sets, nrow(gwas_1))
    cs_2 <- get_cs_index_vector(susie_2$sets, nrow(gwas_2))
    
    table(cs_1)
    table(cs_2)
    
    # 2. Prepare Data for Merging
    res1_df <- data.frame(SNP = gwas_1$SNP, PIP1 = susie_1$pip, CS1 = cs_1, stringsAsFactors = FALSE)
    res2_df <- data.frame(SNP = gwas_2$SNP, PIP2 = susie_2$pip, CS2 = cs_2, stringsAsFactors = FALSE)
    
    # Merge into a single frame
    temp_res <- gwas_1 %>% 
      select(SNP, CHR, POS) %>% 
      full_join(gwas_2 %>% select(SNP, CHR, POS), by = "SNP") %>%
      mutate(CHR = coalesce(CHR.x, CHR.y), POS = coalesce(POS.x, POS.y)) %>%
      select(SNP, CHR, POS) %>%
      left_join(res1_df, by = "SNP") %>%
      left_join(res2_df, by = "SNP") %>%
      mutate(across(c(PIP1, PIP2, CS1, CS2), ~replace_na(., 0))) %>%
      mutate(PIP_Either = pmax(PIP1, PIP2))
    
    # 3. GENERATE 'EITHER' CS INDEX using pmax(LD)
    # Identify all SNPs that are in a CS in *either* pop
    cs_snps_idx <- which(temp_res$CS1 > 0 | temp_res$CS2 > 0)
    
    if(length(cs_snps_idx) > 0) {
      # Extract these SNPs
      active_snps <- temp_res$SNP[cs_snps_idx]
      
      # A. Create the consensus LD matrix
      #    Only active CS SNPs are required.
      idx1 <- match(active_snps, colnames(ld_1))
      idx2 <- match(active_snps, colnames(ld_2))
      
      # Handle cases where SNP might be missing in one LD ref (set NA to 0 for pmax)
      sub_ld1 <- ld_1[idx1, idx1]; sub_ld1[is.na(sub_ld1)] <- 0
      sub_ld2 <- ld_2[idx2, idx2]; sub_ld2[is.na(sub_ld2)] <- 0
      
      # The magic step: Max correlation across ancestries
      ld_either <- pmax(abs(sub_ld1), abs(sub_ld2)) 
      
      # B. Cluster these SNPs into 'Either' Credible Sets
      #    Simple logic: If SNPs are highly correlated in ld_either (>0.6), 
      #    they belong to the same merged signal.
      
      #    Use distance-based clustering:
      #    Distance = 1 - correlation. Cut tree at 0.4 (corr=0.6).
      if (length(active_snps) > 1) {
        h_dynamic <- 1 - sqrt(r2_thresh)
        dist_mat <- as.dist(1 - ld_either)
        
        # Use 'complete' linkage to avoid chaining disparate signals
        hc <- hclust(dist_mat, method = "complete") 
        clusters <- cutree(hc, h = h_dynamic) 
      } else {
        clusters <- 1
      }
      
      # Map clusters back to the main dataframe
      temp_res$CS <- 0
      temp_res$CS[cs_snps_idx] <- clusters
      
    } else {
      temp_res$CS <- 0
    }
    
    mf_result <- temp_res %>%
      mutate(PIP_Shared = pmin(PIP1, PIP2)) %>%
      select(SNP, CHR, POS, PIP_Either, PIP_Shared, PIP1, PIP2, CS, CS1, CS2)
  } else {
    
    # --- CASE 2: MESuSiE (Shared/Common variants) ---
    
    # MESuSiE requires variants present in both datasets
    common_snps <- intersect(gwas_1$SNP, gwas_2$SNP)
    
    g1_sub <- gwas_1[SNP %in% common_snps][order(match(SNP, common_snps))]
    g2_sub <- gwas_2[SNP %in% common_snps][order(match(SNP, common_snps))]
    
    # Subset LD matrices based on SNP IDs (assuming LD has row/col names)
    ld1_sub <- ld_1[common_snps, common_snps]
    ld2_sub <- ld_2[common_snps, common_snps]
    
    summary_stat_sd_list = list(g1_sub, g2_sub)
    names(summary_stat_sd_list) <- pop_names
    R_mat_list = list(ld1_sub, ld2_sub)
    names(R_mat_list) <- pop_names
    
    mesusie_res <- meSuSie_core(R_mat_list, summary_stat_sd_list, 
                                L = L, 
                                prior_weights = prior_weights,
                                ancestry_weight = ancestry_weight)
    
    # Store raw object
    raw_objects$mesusie_res <- mesusie_res
    
    
    # Helper to count CS categories dynamically
    get_pop_cs_vec <- function(mesusie_obj, pop_name) {
      vec <- rep(0, length(mesusie_obj$pip))
      # Extract indices where category matches user-provided pop_name
      pop_sets <- mesusie_obj$cs$cs[mesusie_obj$cs$cs_category == pop_name]
      if (length(pop_sets) > 0) {
        for (j in seq_along(pop_sets)) {
          vec[pop_sets[[j]]] <- j
        }
      }
      return(vec)
    }
    
    cs_either <- get_cs_index_vector(mesusie_res$cs, length(mesusie_res$pip), renumber = TRUE)
    cs_p1 <- get_pop_cs_vec(mesusie_res, pop_names[1])
    cs_p2 <- get_pop_cs_vec(mesusie_res, pop_names[2])
    
    mesusie_df <- data.frame(
      SNP = common_snps,
      PIP_Either = mesusie_res$pip,
      PIP1 = mesusie_res$pip_config[,1],
      PIP2 = mesusie_res$pip_config[,2],
      PIP_Shared = mesusie_res$pip_config[,3],
      CS = cs_either,
      CS1 = cs_p1,
      CS2 = cs_p2
    )
    
    mf_result <- gwas_1 %>% 
      select(SNP, CHR, POS) %>% 
      full_join(gwas_2 %>% select(SNP, CHR, POS), by = "SNP") %>%
      mutate(CHR = coalesce(CHR.x, CHR.y), POS = coalesce(POS.x, POS.y)) %>%
      select(SNP, CHR, POS) %>%
      left_join(mesusie_df, by = "SNP") %>%
      mutate(across(everything(), ~replace_na(., 0))) %>%
      select(SNP, CHR, POS, PIP_Either, PIP_Shared, PIP1, PIP2, CS, CS1, CS2)
  }
  
  # 3. Final Column Formatting
  colnames(mf_result) <- c("SNP", "CHR", "POS", "PIP_Either", "PIP_Shared", 
                           "PIP_Ancestry_1", "PIP_Ancestry_2", "CS", 
                           "CS_Ancestry_1", "CS_Ancestry_2")
  mf_result <- mf_result %>% arrange(CHR, POS)
  
  return(list(
    decision = decision,
    results = mf_result,
    raw_objects = raw_objects
  ))
  
}
