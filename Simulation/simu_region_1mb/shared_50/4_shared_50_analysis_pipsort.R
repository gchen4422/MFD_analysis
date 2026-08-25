args<-as.numeric(commandArgs(TRUE))
LD_BLOCK = args[1]
h2_num = args[2]
num_causal = args[3]

library(data.table)

wrk_dir<-paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50/causal_num_",num_causal,"/")
system(paste0("mkdir -p ",wrk_dir))
data_dir<-paste0(wrk_dir,"summary_data/")
system(paste0("mkdir -p ",data_dir))
result_dir<-paste0(wrk_dir,"result/")
system(paste0("mkdir -p ",result_dir))
paintor_dir<-data_dir


##############################################
#
#		Format MESuSiE Input
#
##############################################
zfile<-read.table(paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num),header=T)
summary_stat_1 = data.frame("SNP" = zfile$RSID, "Beta"=zfile$zscore_1/sqrt(zfile$N_1),"Se"=1/sqrt(zfile$N_1), "Z" =zfile$zscore_1 ,  "N" = zfile$N_1 )
summary_stat_2 = data.frame("SNP" = zfile$RSID, "Beta"=zfile$zscore_2/sqrt(zfile$N_2),"Se"=1/sqrt(zfile$N_2), "Z" =zfile$zscore_2 ,  "N" = zfile$N_2 )

summary_stat_1_bim = read.table(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_eur/loci_",LD_BLOCK,".bim"))
summary_stat_2_bim = read.table(paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_",LD_BLOCK,".bim"))
summary_stat_2_bim$V2 = summary_stat_1_bim$V2


summary_stat_1_allele = data.frame("A1" = summary_stat_1_bim$V5, "A2" = summary_stat_1_bim$V6)
summary_stat_2_allele = data.frame("A1" = summary_stat_2_bim$V5, "A2" = summary_stat_2_bim$V6)

summary_stat_1_susiex = cbind(zfile[,c(1,2,3)],summary_stat_1_allele,summary_stat_1[,c(2:5)])
summary_stat_2_susiex = cbind(zfile[,c(1,2,3)],summary_stat_2_allele,summary_stat_2[,c(2:5)])

summary_stat_1_susiex$P = 2*pnorm(abs(summary_stat_1_susiex$Z),lower.tail = F)
summary_stat_2_susiex$P = 2*pnorm(abs(summary_stat_2_susiex$Z),lower.tail = F)

# colnames(summary_stat_1_susiex)[c(3,4,5)] = c("SNP","REF","ALT")
# colnames(summary_stat_2_susiex)[c(3,4,5)] = c("SNP","REF","ALT")



write.table(summary_stat_1_susiex,file = paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_eur"),col.names = T,row.names = F,quote=F, sep ="\t")
write.table(summary_stat_2_susiex,file = paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_afr"),col.names = T,row.names = F,quote=F, sep ="\t")
write.table(summary_stat_2_bim,file = paste0("/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/risk_loci_ld_afr/loci_",LD_BLOCK,".bim"),col.names = F,row.names = F,quote=F, sep ="\t")



summary_stat_sd_list = list("EUR" = summary_stat_1,"AFR"=summary_stat_2 ) 

EU_cov<-as.matrix(fread(paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")))
colnames(EU_cov)<-zfile$RSID
BB_cov<-as.matrix(fread(paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")))
colnames(BB_cov)<-zfile$RSID
R_mat_list=list("EUR" = EU_cov,"AFR" = BB_cov)

library(dplyr)

paintor_out_name<-paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num)
input_file_name<-paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".input")
ld_EU_paintor_name<-paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD1")
ld_BB_paintor_name<-paste0(paintor_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,".LD2")
paintor_suffix<-paste0("CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num)
input_dir<-paintor_dir


z_file_eur = summary_stat_1_susiex %>% dplyr::select(RSID,Z) 
z_file_afr = summary_stat_2_susiex %>% dplyr::select(RSID,Z)


write.table(z_file_eur,file = paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_z_eur"),col.names = F,row.names = F,quote=F, sep ="\t")
write.table(z_file_afr,file = paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_z_afr"),col.names = F,row.names = F,quote=F, sep ="\t")

line1 <- paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_z_eur")
line2 <- paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_z_afr")
newfile <- data.frame(z_file = c(line1, line2), stringsAsFactors = FALSE)
out_path <- paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_zfiles.txt")
write.table(newfile, file = out_path, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


line1_ld <- ld_EU_paintor_name
line2_ld <- ld_BB_paintor_name
newfile_ld <- data.frame(z_file = c(line1_ld, line2_ld), stringsAsFactors = FALSE)
out_path_ld <- paste0(data_dir, "CAUSAL_", num_causal, "_LOCI_", LD_BLOCK, "_h2_", h2_num, "_ldfiles.txt")
write.table(newfile_ld, file = out_path_ld, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")




snp_map = summary_stat_1_susiex %>% dplyr::select(RSID) %>% mutate(SNP1 = 0:(nrow(summary_stat_1_susiex)-1), SNP2 =0:(nrow(summary_stat_1_susiex)-1) )
snp_map_name = paste0(data_dir,"CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_snp_map")
write.table(snp_map,file =snp_map_name ,col.names = F,row.names = F,quote=F, sep =",")


pipsort_outname = paste0(result_dir,paintor_suffix,"_pipsort")


#limit the maximum number of causal variants at 1 if either study has more than 1700 variants and at 2 if either study has more than 325 variants

n_max <- max(nrow(summary_stat_1), nrow(summary_stat_2))
#causal_number_pipsort <- if (n_max > 1700) 1 else if (n_max > 325) 2 else causal_number_pipsort
causal_number_pipsort = 2

cmd <- sprintf(
  '/home/chen4422/gwas_software/pipsort/PIPSORT -c %d -l %s -z %s -m %s -n 300000,300000 -p 0.25 -o %s',
  causal_number_pipsort, shQuote(out_path_ld), shQuote(out_path), shQuote(snp_map_name), shQuote(pipsort_outname)
)



start.time<-Sys.time()

system(cmd)

pipsort.time<-Sys.time()-start.time
pipsort.time<-as.numeric(pipsort.time,units = "mins")






py <- "/apps/spack/negishi/apps/anaconda/2022.10-py39-gcc-8.5.0-sjvibry/bin/python"
script <- "/home/chen4422/gwas_software/pipsort/utils/get_global_pips.py"
args <- c(
  paste0(pipsort_outname, "_study0_post.txt"),
  paste0(pipsort_outname, "_study1_post.txt"),
  paste0(pipsort_outname, "_shared_pips.txt"),
  paste0(pipsort_outname, "_global_pips.txt")
)

system2(py, args = c(script, args))

script <- "/home/chen4422/gwas_software/pipsort/utils/get_not_shared_pips.py"
system2(
  command = py,
  args = c(
    script,
    paste0(pipsort_outname, "_shared_pips.txt"),
    paste0(pipsort_outname, "_global_pips.txt"),
    paste0(pipsort_outname, "_not_shared_pips.txt")
  )
)

time_out_name<-paste0(result_dir,"time_CAUSAL_",num_causal,"_LOCI_",LD_BLOCK,"_h2_",h2_num,"_pipsort.txt")
write.table(pipsort.time,time_out_name,col.names = F,row.names = F,quote=F)
