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
  eur_ldname = paste("./risk_loci_ld_eur/loci_1kg", i,sep = "_")
  afr_ldname = paste("./risk_loci_ld_afr/loci_1kg", i,sep = "_")
  out_ldname = paste("./risk_loci_ld_merge/loci_1kg", i,sep = "_")
  plink_bin <- "~/gwas_software/plink2/plink"   # plink v1.9 binary
  ref_filename = paste0("./risk_loci_ld_eur/loci_1kg_", i,".bim")
  
  
  command <- sprintf(
    '%s --bfile "%s" --bmerge "%s.bed" "%s.bim" "%s.fam" --make-bed --allow-no-sex --a1-allele %s 5 2 --out "%s"',
    plink_bin, eur_ldname, afr_ldname, afr_ldname, afr_ldname,ref_filename, out_ldname
  )
  t1 <- system(command, intern = TRUE)
  
  
  #Check LD and summary stats allele consistency
  info = fread(paste0("risk_loci_ld_merge/loci_1kg_",i,".bim"))
  merged_plink <- read.plink(paste0(out_ldname,".bed"))
  merged_plink_geno <- as(merged_plink$genotypes, "numeric")
  merged_mat = as.matrix(merged_plink_geno)
  merged_mat = preprocess_data(merged_mat)
  
  compute_cov(merged_mat, paste0(out_ldname,".ld"))
  #save(merged_mat, file = paste0("risk_loci_ld/loci_",i,".Rdata"))
  
  
  
}

test_eur = fread("./risk_loci_ld_eur/loci_1kg_1.bim")
test_afr = fread("./risk_loci_ld_afr/loci_1kg_1.bim")
test_merge = fread("./risk_loci_ld_merge/loci_1kg_1.bim")

  
  
  