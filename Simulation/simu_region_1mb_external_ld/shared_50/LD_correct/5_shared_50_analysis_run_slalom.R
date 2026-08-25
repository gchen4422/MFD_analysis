library(Rcpp)
#sourceCpp("/net/fantasia/home/borang/MESuSiE_inf/software/develop_mode/MESuSiE_inf.cpp")
#source("/net/fantasia/home/borang/MESuSiE_inf/software/develop_mode/MESuSiE_inf.R")
library(dplyr)
library(data.table)
#source("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/MESuSiE_original/modified_kriging_rss.R")


LD_BLOCK_range <- 1:100
h2_num_range   <- 1:2
num_causal_range <- 1:3



for (num_causal in num_causal_range) {
      for (LD_BLOCK in LD_BLOCK_range) {
        for (h2_num in h2_num_range) {
          
          wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb_external_ld/shared_50/causal_num_",num_causal,"/")
          system(paste0("mkdir -p ",wrk_dir))
          data_dir<-paste0(wrk_dir,"summary_data/")
          output_dir<-paste0(data_dir,"slalom/")
          system(paste0("mkdir -p ",output_dir))
          
          
          
          
          summary_stat_1_susiex = fread(paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_eur"))
          summary_stat_2_susiex = fread(paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_afr"))
          
          summary_stat_1_susiex <- summary_stat_1_susiex %>%
            dplyr::rename(
              chromosome = CHR,
              position   = POS,
              rsid       = RSID,   # not required but useful / consistent with slalom
              allele1    = A1,
              allele2    = A2,
              beta       = Beta,
              se         = Se,
              n_nfe      = N,
              p          = P
            )
          
          summary_stat_2_susiex <- summary_stat_2_susiex %>%
            dplyr::rename(
              chromosome = CHR,
              position   = POS,
              rsid       = RSID,
              allele1    = A1,
              allele2    = A2,
              beta       = Beta,
              se         = Se,
              n_afr      = N,
              p          = P
            )
          
          
          
          write.table(summary_stat_1_susiex,file = paste0(output_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_eur.snp"),col.names = T,row.names = F,quote=F, sep ="\t")
          write.table(summary_stat_2_susiex,file = paste0(output_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_afr.snp"),col.names = T,row.names = F,quote=F, sep ="\t")
          #write.table(summary_stat_2_bim,file = paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_",LD_BLOCK,".bim"),col.names = F,row.names = F,quote=F, sep ="\t")
          
          
          
          
          
        }
      }
    }

