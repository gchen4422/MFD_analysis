args<-as.numeric(commandArgs(TRUE))
LD_BLOCK = args[1]
External_missing_index = args[2]

library(data.table)
library(snpStats)
library(mvtnorm)
library(dplyr)

LD_dir<-paste0("/net/fantasia/home/borang/Susie_Mult/Revision_Round_1/Simulation/091223/data/LD_data/")
ld_EU_name<-paste0(LD_dir,"Region_",LD_BLOCK,".LD1")
EU_cov<-as.matrix(fread(ld_EU_name))

ld_BB_name<-paste0(LD_dir,"Region_",LD_BLOCK,".LD2")
BB_cov<-as.matrix(fread(ld_BB_name))
if(External_missing_index ==1 ){
	num_snp_missing = round(nrow(BB_cov)*0.3)
}else if(External_missing_index ==2){
	num_snp_missing = round(nrow(BB_cov)*0.3+nrow(BB_cov)*(1-0.3)*0.3)
}

for(num_causal in 1:3){
  
  input_dir<-paste0("/net/fantasia/home/borang/Susie_Mult/Revision_Round_1/Simulation/091223/data/shared_50/causal_num_",num_causal,"/summary_data/")
	
	for(h2_num in 1:2){ 
		z_file<-read.table(paste0(input_dir,"CAUSAL_",num_causal,"_GENE_",LD_BLOCK,"_h2_",h2_num),header=T)
		non_causal_SNP = z_file%>%filter(Signal==0)%>%pull(RSID)
		EU_missing_SNP = sample(non_causal_SNP,num_snp_missing)
		BB_missing_SNP = EU_missing_SNP
		z_file<-z_file%>%mutate(EU_missing = ifelse(RSID%in%EU_missing_SNP,1,0),BB_missing = ifelse(RSID%in%BB_missing_SNP,1,0))
	 
	  
	  for(causal_index in c("Both","One")){
		if(External_missing_index ==1){
			wrk_dir<-paste0("/net/fantasia/home/borang/Susie_Mult/Revision_Round_1/Simulation/091223/data/Missing_Non_Causal/",causal_index,"/causal_num_",num_causal,"/")
		}else if(External_missing_index ==2){
			wrk_dir<-paste0("/net/fantasia/home/borang/Susie_Mult/Revision_Round_1/Simulation/091223/data/Missing_Non_Causal/External_",causal_index,"/causal_num_",num_causal,"/")
		}
	  system(paste0("mkdir -p ",wrk_dir))
	  data_dir<-paste0(wrk_dir,"summary_data/")
	  result_dir<-paste0(wrk_dir,"result/")
	  out_dir<-paste0(wrk_dir,"out/")
	  system(paste0("mkdir -p ",data_dir))
	  system(paste0("mkdir -p ",result_dir))
	  system(paste0("mkdir -p ",out_dir))
	  

     
    z_file_out<-paste0(data_dir,"CAUSAL_",num_causal,"_GENE_",LD_BLOCK,"_h2_",h2_num)
    write.table(z_file,z_file_out,col.names = T,row.names = F,quote=F,sep=" ")  
    
    ld_EU_paintor_name<-paste0(data_dir,"CAUSAL_",num_causal,"_GENE_",LD_BLOCK,"_h2_",h2_num,".LD1")
    ld_BB_paintor_name<-paste0(data_dir,"CAUSAL_",num_causal,"_GENE_",LD_BLOCK,"_h2_",h2_num,".LD2")
    fwrite(EU_cov,ld_EU_paintor_name,sep=" ",col.names=F)
    fwrite(BB_cov,ld_BB_paintor_name,sep=" ",col.names=F)
    
    annotation_file<-matrix(rep(1,nrow(z_file)),ncol=1)
    colnames(annotation_file)<-"coding"
    annotation_paintor_name<-paste0(data_dir,"CAUSAL_",num_causal,"_GENE_",LD_BLOCK,"_h2_",h2_num,".annotations")
    write.table(annotation_file,annotation_paintor_name,col.names = T,row.names = F,quote=F)
    
    input_file_name<-paste0(data_dir,"CAUSAL_",num_causal,"_GENE_",LD_BLOCK,"_h2_",h2_num,".input")
    write.table(paste0("CAUSAL_",num_causal,"_GENE_",LD_BLOCK,"_h2_",h2_num),input_file_name,col.names = F,row.names = F,quote=F)
    }	
  }
}

