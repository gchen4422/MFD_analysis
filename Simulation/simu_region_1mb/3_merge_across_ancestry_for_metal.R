library(data.table)
library(dplyr)
library(snpStats)
library(MESuSiE)


preprocess_data <- function(data_matrix) {
  # Replace missing values with column means and scale columns
  data_matrix <- apply(data_matrix, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x <- scale(x)
    return(x)
  })
  
  return(data_matrix)
}

# Function to compute covariance matrix and write it to a file
compute_cov <- function(data_matrix, file_name) {
  # Compute the covariance matrix
  cov_matrix <- cov2cor(crossprod(data_matrix))
  
  # Ensure the covariance matrix is symmetric
  cov_matrix <- (cov_matrix + t(cov_matrix)) / 2
  
  # Write the covariance matrix to a file
  fwrite(cov_matrix, file_name,sep =" ")
}





for(i in 1:100){
  
  snplist = paste("loci",i, sep = "_")
  eur_ldname = paste("./risk_loci_ld_eur/loci", i,sep = "_")
  afr_ldname = paste("./risk_loci_ld_afr/loci", i,sep = "_")
  out_ldname = paste("./risk_loci_ld_merge/loci", i,sep = "_")
  plink_bin <- "~/gwas_software/plink2/plink"   # plink v1.9 binary
  ref_filename = paste0("./risk_loci_ld_eur/loci_", i,".bim")
  
  
  command <- sprintf(
    '%s --bfile "%s" --bmerge "%s.bed" "%s.bim" "%s.fam" --make-bed --allow-no-sex --a1-allele %s 5 2 --out "%s"',
    plink_bin, eur_ldname, afr_ldname, afr_ldname, afr_ldname,ref_filename, out_ldname
  )
  t1 <- system(command, intern = TRUE)
  
  
  #Check LD and summary stats allele consistency
  info = fread(paste0("risk_loci_ld_merge/loci_",i,".bim"))
  merged_plink <- read.plink(paste0(out_ldname,".bed"))
  merged_plink_geno <- as(merged_plink$genotypes, "numeric")
  merged_mat = as.matrix(merged_plink_geno)
  merged_mat = preprocess_data(merged_mat)
  
  compute_cov(merged_mat, paste0(out_ldname,".ld"))
  #save(merged_mat, file = paste0("risk_loci_ld/loci_",i,".Rdata"))
  
  
  
}


library(pheatmap)


merged_ld = as.matrix(fread("risk_loci_ld_merge/loci_1.ld"))


pheatmap(
  merged_ld,
  color = colorRampPalette(c("white","orange","red"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  legend = TRUE,
  display_numbers = FALSE,
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  breaks = seq(0, 1, length.out = 101)  # use seq(-1,1,101) if it's r not r^2
)

eur_ld = as.matrix(fread("risk_loci_ld_eur/loci_1.ld"))

pheatmap(
  eur_ld,
  color = colorRampPalette(c("white","orange","red"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  legend = TRUE,
  display_numbers = FALSE,
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  breaks = seq(0, 1, length.out = 101)  # use seq(-1,1,101) if it's r not r^2
)


afr_ld = as.matrix(fread("risk_loci_ld_afr/loci_1.ld"))

pheatmap(
  afr_ld,
  color = colorRampPalette(c("white","orange","red"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  legend = TRUE,
  display_numbers = FALSE,
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  breaks = seq(0, 1, length.out = 101)  # use seq(-1,1,101) if it's r not r^2
)

n_eur = 50000
n_afr = 6374

weighted_ld <- (n_eur * eur_ld + n_afr * afr_ld) / (n_eur + n_afr)


pheatmap(
  weighted_ld,
  color = colorRampPalette(c("white","orange","red"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  legend = TRUE,
  display_numbers = FALSE,
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  breaks = seq(0, 1, length.out = 101)  # use seq(-1,1,101) if it's r not r^2
)

n_eur = 1
n_afr = 1

weighted_ld <- (n_eur * eur_ld + n_afr * afr_ld) / (n_eur + n_afr)


pheatmap(
  weighted_ld,
  color = colorRampPalette(c("white","orange","red"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  legend = TRUE,
  display_numbers = FALSE,
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE,
  breaks = seq(0, 1, length.out = 101)  # use seq(-1,1,101) if it's r not r^2
)


quantify_ld_diff <- function(merged_ld, eur_ld) {
  stopifnot(all(dim(merged_ld) == dim(eur_ld)))
  # use only one triangle, exclude diagonal
  tri <- lower.tri(merged_ld, diag = FALSE)
  x <- as.numeric(merged_ld[tri])
  y <- as.numeric(eur_ld[tri])
  
  # basic diffs
  d  <- x - y
  mae  <- mean(abs(d), na.rm = TRUE)
  mse  <- mean(d^2, na.rm = TRUE)
  rmse <- sqrt(mse)
  max_abs_diff <- max(abs(d), na.rm = TRUE)
  
  # norms (Frobenius on triangles) + relative error
  frob_diff <- sqrt(sum(d^2, na.rm = TRUE))
  frob_ref  <- sqrt(sum(y^2, na.rm = TRUE))  # relative to eur_ld
  rel_frob  <- frob_diff / (frob_ref + 1e-12)
  
  # agreement/correlation of entries
  r_pearson  <- suppressWarnings(cor(x, y, method = "pearson",  use = "complete.obs"))
  r_spearman <- suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
  r2         <- r_pearson^2
  
  # Bland–Altman style summary
  mean_diff <- mean(d, na.rm = TRUE)
  sd_diff   <- stats::sd(d, na.rm = TRUE)
  
  # return a tidy list
  list(
    n_pairs            = sum(tri, na.rm = TRUE),
    mae                = mae,
    rmse               = rmse,
    max_abs_diff       = max_abs_diff,
    frobenius_diff     = frob_diff,
    relative_frobenius = rel_frob,
    pearson_r          = r_pearson,
    spearman_r         = r_spearman,
    r_squared          = r2,
    mean_diff          = mean_diff,
    sd_diff            = sd_diff
  )
}

# Example:
res <- quantify_ld_diff(merged_ld, eur_ld)
str(res)

res <- quantify_ld_diff(merged_ld, afr_ld)
str(res)

res <- quantify_ld_diff(merged_ld, weighted_ld)
str(res)

res <- quantify_ld_diff(eur_ld, afr_ld)
str(res)


