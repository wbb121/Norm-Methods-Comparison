########################################################################################################################
### validate the performance of different normalization methods on binary phenotype prediction in CRC/IBD datasets
### 8 available CRC datasets: "FengQ_2015","GuptaA_2019","ThomasAM_2018","VogtmannE_2016",
#                             "WirbelJ_2018","YachidaS_2019","YuJ_2015","ZellerG_2014"
### 5 available IBD datasets: "HallAB_2017","HMP_2019_ibdmdb","IjazUZ_2017","VilaAV_2018","NielsenHB_2014"
### repeat 30 times
### 2023/06/05
########################################################################################################################

setwd("/home/wangbb/normalization_comparison")

### packages
all(sapply(c("foreach","doParallel"), require, character.only=TRUE))

### load the helper function
source("helper.R")


#======================================================================================================================#
### arguments
command_args=commandArgs(trailingOnly=T)
#command_args=c("data/crc_meta.rds","data/crc_count.rds","TSS",10,"ranger")

### parameters
meta=readRDS(command_args[1])
count=readRDS(command_args[2])
norm_method=command_args[3]
pred_cluster=as.numeric(command_args[4])
pred_method=command_args[5]

studies <- unique(meta$study_name)
meta$phenotype <- meta$study_condition
meta$phenotype <- sub(unique(meta$phenotype)[unique(meta$phenotype)!="control"],"case",meta$phenotype)


#======================================================================================================================#
### normalization
for(study1 in studies){
  for(study2 in studies){
    if(study1!=study2){
      # data to be normalized
      meta1 <- meta[meta$study_name==study1,]
      meta2 <- meta[meta$study_name==study2,]
      data1 <- count[,colnames(count)%in%meta1$sample_id]
      data2 <- count[,colnames(count)%in%meta2$sample_id]
      # normalization
      norm_data <- norm.func(p1=data1,p2=data2,norm_method=norm_method)
      # save the results
      if(!dir.exists("real_data_norm")) dir.create("real_data_norm")
      if(!dir.exists(paste0("./real_data_norm/",norm_method))) dir.create(paste0("./real_data_norm/",norm_method))
      saved_file <- paste0("./real_data_norm/",norm_method,"/",norm_method,"_trn_",study1,"_tst_",study2,".rds")
      saveRDS(norm_data, saved_file)
      print(paste0("trn ",study1,", tst ",study2,", ",norm_method))
    }
  }
}


#======================================================================================================================#
### prediction ###
for(study1 in crc_study){
  for(study2 in crc_study){
    if(study1!=study2){
      
      # normalized data
      norm_data <- readRDS(paste0("./real_data_norm/",norm_method,"/",norm_method,"_trn_",study1,"_tst_",study2,".rds"))
      
      # prediction
      cl<- makeCluster(pred_cluster,setup_strategy="sequential")   
      registerDoParallel(cl)
      auc_values <- foreach(x=1:30,.packages=c("caret","pROC"),.combine=c,.errorhandling="pass") %dopar% {
        pred.func(trn=norm_data[[1]],tst=norm_data[[2]],meta=meta,pred_method=pred_method)
      }
      stopCluster(cl)
      
      # summarize the results
      df <- data.frame(study1=study1,study2=study2,norm_method=norm_method,pred_method=pred_method,auc_values=auc_values)
      
      # save the results
      if(!dir.exists("./real_data_pred")) dir.create("./real_data_pred")
      saved_file <- paste0("./real_data_pred/",pred_method,"_",norm_method,"_trn_",study1,"_tst_",study2,".rds")
      saveRDS(df, saved_file)
      print(paste0("trn ",study1,", tst ", study2,", ",norm_method,", ",pred_method))

    }
  }
}





