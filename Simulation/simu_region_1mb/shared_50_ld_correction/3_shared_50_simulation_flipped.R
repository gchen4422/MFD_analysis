# Loading necessary libraries
library(data.table)
library(snpStats)
library(mvtnorm)
library(dplyr)

# Setting constants and initializing parameters
#args <- as.numeric(commandArgs(TRUE))
#LD_BLOCK <- args[1]
num_EU <- 300000
num_BB <- 300000



for (LD_BLOCK in 1:100){
  
  
  for(num_causal in 1:3){
    
    
    data_dir_summ<-paste0(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_",num_causal,"/"),"summary_data/")
    wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_ld_correction/causal_num_",num_causal,"/")
    system(paste0("mkdir -p ",wrk_dir))
    data_dir<-paste0(wrk_dir,"summary_data/")
    result_dir<-paste0(wrk_dir,"result/")
    out_dir<-paste0(wrk_dir,"out/")
    system(paste0("mkdir -p ",data_dir))
    system(paste0("mkdir -p ",result_dir))
    system(paste0("mkdir -p ",out_dir))

    
    

    print(LD_BLOCK)
    for(h2_num in 1:2){
      
      
      zfile<-read.table(paste0(data_dir_summ,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),header=T)
      EU_cov<-as.matrix(fread(paste0(data_dir_summ,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")))
      BB_cov<-as.matrix(fread(paste0(data_dir_summ,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")))
      
      
      
      for (flip_prop_num in 1:3){
      
        
        flip_prop = c(0.001,0.01,0.1)[flip_prop_num]
        k_eu <- ceiling(length(zfile$zscore_1) * flip_prop)           ###Think about this as rsprasepro only select non-causal snps to be mismatched
        k_bb <- ceiling(length(zfile$zscore_2) * flip_prop)
        
      
        for (rho_num in 1:4){
          
          z_EU_copy = zfile$zscore_1
          z_BB_copy = zfile$zscore_2
          
          rho = c(0.5,0,-0.5,-1)[rho_num]
          idx_eu <- sample(seq_along(1:length(zfile$zscore_1)), k_eu)
          idx_bb <- sample(seq_along(1:length(zfile$zscore_2)), k_bb)
          z_eu_samp <- zfile$zscore_1[idx_eu]
          z_bb_samp <- zfile$zscore_2[idx_bb]
          
          
          z_EU_copy[idx_eu] = z_eu_samp*rho
          z_BB_copy[idx_bb] = z_bb_samp*rho
          
          
          flip_eu_ind <- seq_len(nrow(zfile)) %in% idx_eu
          flip_bb_ind <- seq_len(nrow(zfile)) %in% idx_bb
          
          
          z_file_flipped <- zfile %>%
            mutate(
              zscore_1   = z_EU_copy,
              zscore_2   = z_BB_copy,
              flipped_EU = as.integer(flip_eu_ind),  # 1 = flipped in EU, 0 = not
              flipped_BB = as.integer(flip_bb_ind)   # 1 = flipped in BB, 0 = not
            )
          
          
          data_dir_flipped = paste0(data_dir,"flipped_",flip_prop_num,"_rho_",rho_num,"/")
          system(paste0("mkdir -p ",data_dir_flipped))
          
          
          z_file_out <- paste0(
            data_dir_flipped,
            "CAUSAL_", num_causal,
            "_LOCI_",  LD_BLOCK,
            "_h2_",    h2_num,
            "_flipped_",flip_prop_num,
            "_rho_",   rho_num
          )
          
          write.table(z_file_flipped,z_file_out,col.names = T,row.names = F,quote=F,sep=" ")  
        }
      }
      
      
      ld_EU_paintor_name<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")
      ld_BB_paintor_name<-paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")
      fwrite(EU_cov,ld_EU_paintor_name,sep=" ",col.names=F)
      fwrite(BB_cov,ld_BB_paintor_name,sep=" ",col.names=F)
    }
  }
  
}
