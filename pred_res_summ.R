######################################################################################################################
### summarize the results
### 2023/06/05
######################################################################################################################

setwd("/home/wangbb/normalization_comparison")

library(reshape2)


#=====================================================================================================================#
### information for normalization methods ###
#=====================================================================================================================#
norm_methods <- data.frame(
  norm_method=c("TSS","UQ","MED","CSS","TMM","RLE_poscounts","GMPR","CLR",
                "LOG","AST","STD","rank","blom","NPN","logcpm","VST",
                "QN","FSQN","BMC","limma","combat","conqur"),
  annotated_method=c("TSS","UQ","MED","CSS","TMM","RLE","GMPR","CLR",
                     "LOG","AST","STD","Rank","Blom","NPN","logCPM","VST",
                     "QN","FSQN","BMC","Limma","ComBat","ConQuR"),
  class=c(rep("Scaling",7),rep("CoDA",1),rep("Transformation",8),rep("Batch Correction",6))
)
norm_methods$annotated_method <- factor(norm_methods$annotated_method,levels=norm_methods$annotated_method)
norm_methods$class <- factor(norm_methods$class,levels=c("Scaling","CoDA","Transformation","Batch Correction"))
saveRDS(norm_methods,"data/norm_methods.rds")


#=====================================================================================================================#
### prediction results for simulated data ###
#=====================================================================================================================#
### simulation parameters
ep <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.4, 0.6, 0.8, 1)  
ed <- c(1.02,1.04,1.06)
nd <- 10
pred_method <- "ranger"

### combine the prediction results in "E:/USC/normalization/20230419_sim_crc_ranger_gpu/"
sim_pred <- data.frame(alpha=NA,ed=NA,nd=NA,norm_method=NA,pred_method=NA,auc_value=NA)
sim_pred <- sim_pred[-1,]
# results with normal ep
for(i in ed){
  for(j in norm_methods$norm_method){
    pred_res_file <- paste0("./sim_data_pred/auc_nd",nd,"_ed",i,"_",j,"_",pred_method,".rds")
    pred_res <- readRDS(pred_res_file)
    sim_pred <- rbind(sim_pred,pred_res)
  }
}
# save the results
saveRDS(sim_pred,"data/sim_pred.rds")

### compute the average auc value for different combinations
sim_pred_summ <- data.frame(alpha=NA,ed=NA,nd=NA,pred_method=NA,norm_method=NA,average_auc=NA)
sim_pred_summ <- sim_pred_summ[-1,]
for(i in unique(sim_pred$alpha)){
  for(j in unique(sim_pred$ed)){
    for(k in unique(sim_pred$norm_method)){
      data <- sim_pred[sim_pred$alpha==i & sim_pred$ed==j & sim_pred$norm_method==k, ]
      df <- data.frame(alpha=i,ed=j,nd=nd,pred_method=pred_method,norm_method=k,
                       average_auc=round(mean(as.numeric(data$auc_value)),5))
      sim_pred_summ <- rbind(sim_pred_summ,df)
    }
  }
}

### add the orders of auc_values
final_data <- data.frame(alpha=NA,ed=NA,nd=NA,pred_method=NA,norm_method=NA,average_auc=NA,order_average_auc=NA)
final_data <- final_data[-1,]
for(i in unique(sim_pred_summ$alpha)){
  for(j in unique(sim_pred_summ$ed)){
    # i=0;j=1.025
    data <- sim_pred_summ[sim_pred_summ$alpha==i & sim_pred_summ$ed==j,]
    # orders for average auc
    data <- data[order(data$average_auc,decreasing=T),]
    data$order_average_auc <- 1:nrow(data)
    # final data
    final_data <- rbind(final_data,data)
  }
}
saveRDS(final_data,"data/sim_pred_summ.rds")


#=====================================================================================================================#
### prediction results for CRC data ###
#=====================================================================================================================#
### related information for CRC datasets 
crc_meta <- readRDS("data/crc_meta.rds")
crc_studies <- unique(crc_meta$study_name)
norm_methods <- readRDS("data/norm_methods.rds")

### results of cross-study prediction using ranger
crc_pred <- data.frame(study1=NA,study2=NA,norm_method=NA,pred_method=NA,auc_values=NA)
crc_pred <- crc_pred[-1,]
for(i in norm_methods$norm_method){
  for(j in crc_studies){
    for(k in crc_studies){
      if(j!=k){
        pred_res_file <- paste0("./real_data_pred/",pred_method,"_",i,"_trn_",j,"_tst_",k,".rds")
        pred_res <- readRDS(pred_res_file)
        crc_pred <- rbind(crc_pred,pred_res)
      }
    }
  }
}
# save the results
saveRDS(crc_pred,"data/crc_pred.rds")

### compute the mean auc value for different combinations
crc_pred_summ <- data.frame(study1=NA,study2=NA,pred_method=NA,norm_method=NA,average_auc=NA)
crc_pred_summ <- crc_pred_summ[-1,]
for(i in unique(crc_pred$study1)){
  for(j in unique(crc_pred$study2)){
    for(k in unique(crc_pred$norm_method)){
      data <- crc_pred[crc_pred$study1==i & crc_pred$study2==j & crc_pred$norm_method==k, ]
      df <- data.frame(study1=i,study2=j,pred_method="crc",norm_method=k,
                       average_auc=round(mean(as.numeric(data$auc_values)),5))
      crc_pred_summ <- rbind(crc_pred_summ,df)
    }
  }
}

### add the orders of auc_values
final_data <- data.frame(study1=NA,study2=NA,pred_method=NA,norm_method=NA,average_auc=NA,order_average_auc=NA)
final_data <- final_data[-1,]
for(i in unique(crc_pred_summ$study1)){
  for(j in unique(crc_pred_summ$study2)){
    # i="FengQ_2015";j="ZellerG_2014"
    data <- crc_pred_summ[crc_pred_summ$study1==i & crc_pred_summ$study2==j,]
    # orders for average auc
    data <- data[order(data$average_auc,decreasing=T),]
    data$order_average_auc <- 1:nrow(data)
    # final data
    final_data <- rbind(final_data,data)
  }
}
# save the results
saveRDS(final_data,"data/crc_pred_summ.rds")


#=====================================================================================================================#
### prediction results for ibd data ###
#=====================================================================================================================#
### related information for ibd datasets 
ibd_meta <- readRDS("data/ibd_meta.rds")
ibd_studies <- unique(ibd_meta$study_name)
norm_methods <- readRDS("data/norm_methods.rds")

### results of cross-study prediction using ranger
ibd_pred <- data.frame(study1=NA,study2=NA,norm_method=NA,pred_method=NA,auc_values=NA)
ibd_pred <- ibd_pred[-1,]
for(i in norm_methods$norm_method){
  for(j in ibd_studies){
    for(k in ibd_studies){
      if(j!=k){
        pred_res_file <- paste0("./real_data_pred/",pred_method,"_",i,"_trn_",j,"_tst_",k,".rds")
        pred_res <- readRDS(pred_res_file)
        crc_pred <- rbind(crc_pred,pred_res)
      }
    }
  }
}
# save the results
saveRDS(ibd_pred,"data/ibd_pred.rds")

### compute the mean auc value for different combinations
ibd_pred_summ <- data.frame(study1=NA,study2=NA,pred_method=NA,norm_method=NA,average_auc=NA)
ibd_pred_summ <- ibd_pred_summ[-1,]
for(i in unique(ibd_pred$study1)){
  for(j in unique(ibd_pred$study2)){
    for(k in unique(ibd_pred$norm_method)){
      data <- ibd_pred[ibd_pred$study1==i & ibd_pred$study2==j & ibd_pred$norm_method==k, ]
      df <- data.frame(study1=i,study2=j,pred_method="ibd",norm_method=k,
                       average_auc=round(mean(as.numeric(data$auc_values)),5))
      ibd_pred_summ <- rbind(ibd_pred_summ,df)
    }
  }
}

### add the orders of auc_values
final_data <- data.frame(study1=NA,study2=NA,pred_method=NA,norm_method=NA,average_auc=NA,order_average_auc=NA)
final_data <- final_data[-1,]
for(i in unique(ibd_pred_summ$study1)){
  for(j in unique(ibd_pred_summ$study2)){
    # i="FengQ_2015";j="ZellerG_2014"
    data <- ibd_pred_summ[ibd_pred_summ$study1==i & ibd_pred_summ$study2==j,]
    # orders for average auc
    data <- data[order(data$average_auc,decreasing=T),]
    data$order_average_auc <- 1:nrow(data)
    # final data
    final_data <- rbind(final_data,data)
  }
}
# save the results
saveRDS(final_data,"data/ibd_pred_summ.rds")




