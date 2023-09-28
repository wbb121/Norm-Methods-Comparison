######################################################################################################################
### summarize the results
### 2023/06/05
######################################################################################################################

setwd("/home/wangbb/normalization_comparison")

library(reshape2)

#=====================================================================================================================#
### normalization methods ###
#=====================================================================================================================#
### normalization methods 
norm_methods <- data.frame(
  norm_method=c("TSS","UQ","MED","CSS","TMM","RLE_poscounts","GMPR","CLR+",
                "logcpm","LOG","AST","STD","rank","blom","VST","NPN",
                "QN","FSQN","BMC","limma","combat","conqur"),
  annotated_method=c("TSS","UQ","MED","CSS","TMM","RLE","GMPR","CLR",
                     "logCPM","LOG","AST","STD","Rank","Blom","VST","NPN",
                     "QN","FSQN","BMC","Limma","ComBat","ConQuR"),
  class=c(rep("Scaling",7),rep("CoDA",1),rep("Transformation",8),rep("Batch Correction",6))
)
norm_methods$annotated_method <- factor(norm_methods$annotated_method,
                                        levels=c("TSS","UQ","MED","CSS","TMM","RLE","GMPR","CLR",
                                                 "LOG","AST","STD","Rank","Blom","NPN","logCPM","VST",
                                                 "QN","FSQN","BMC","Limma","ComBat","ConQuR"))
norm_methods$class <- factor(norm_methods$class,levels=c("Scaling","CoDA","Transformation","Batch Correction"))
saveRDS(norm_methods,"res_summ/norm_methods.rds")


#=====================================================================================================================#
### simulation scenario1 ###
#=====================================================================================================================#
### parameters
population_effects <- c(0,0.05,0.1,0.15,0.2,0.25,0.4,0.6,0.8,1)  # population effect
disease_effects <- c(1.02,1.04,1.06)                             # disease effect
norm_methods <- readRDS("res_summ/norm_methods.rds")
pred_method <- "rfr"

### prediction results
scenario1_df <- data.frame(ep=numeric(0),ed=numeric(0),norm_method=character(0),pred_method=character(0),auc=numeric(0))
scenario1_summ_df <- data.frame(ep=numeric(0),ed=numeric(0),norm_method=character(0),pred_method=character(0),average_auc=numeric(0))
for(alpha in population_effects){
  for(disease_effect in disease_effects){
    for(norm_method in norm_methods$norm_method){
      #alpha=0;disease_effect=1.02;norm_method="TSS"
      df <- readRDS(paste0("scenario1/pred_results/",norm_method,"/norm_alpha",alpha,"_ed",disease_effect,"_",norm_method,"_",pred_method,".rds"))
      summ_df <- data.frame(ep=alpha,ed=disease_effect,norm_method=norm_method,pred_method=pred_method,
                            average_auc=round(mean(as.numeric(df$auc_value)),3))
      scenario1_df <- rbind(scenario1_df,df)
      scenario1_summ_df <- rbind(scenario1_summ_df,summ_df)
    }
  }
}
saveRDS(scenario1_df,"res_summ/scenario1_res.rds")
saveRDS(scenario1_summ_df,"res_summ/scenario1_res_summ.rds")


#=====================================================================================================================#
### simulation scenario2 ###
#=====================================================================================================================#
### parameters
parameters <- list(c(0,1),c(500,1),c(1000,1),c(0,2),c(0,4))
disease_effects <- c(1.02,1.04,1.06)
norm_methods <- readRDS("res_summ/norm_methods.rds")
pred_method <- "rfr"

### prediction results
scenario2_df <- data.frame(ed=numeric(0),batch_mean=numeric(0),batch_var=numeric(0),
                           norm_method=character(0),pred_method=character(0),auc=numeric(0))
scenario2_summ_df <- data.frame(ed=numeric(0),batch_mean=numeric(0),batch_var=numeric(0),
                                norm_method=character(0),pred_method=character(0),average_auc=numeric(0))
for(i in 1:length(parameters)){
  for(disease_effect in disease_effects){
    for(norm_method in norm_methods$norm_method){
      #alpha=0;disease_effect=1.02;norm_method="TSS"
      df <- readRDS(paste0("scenario2/pred_results/",norm_method,"/pred_ed",disease_effect,"_mean",parameters[[i]][1],"_var",parameters[[i]][2],"_",norm_method,"_",pred_method,".rds"))
      summ_df <- data.frame(ed=disease_effect,batch_mean=parameters[[i]][1],batch_var=parameters[[i]][2],
                            norm_method=norm_method,pred_method=pred_method,average_auc=round(mean(as.numeric(df$auc_value)),3))
      scenario2_df <- rbind(scenario2_df,df)
      scenario2_summ_df <- rbind(scenario2_summ_df,summ_df)
    }
  }
}
saveRDS(scenario2_df,"res_summ/scenario2_res.rds")
saveRDS(scenario2_summ_df,"res_summ/scenario2_res_summ.rds")


#=====================================================================================================================#
### simulation scenario3 ###
#=====================================================================================================================#
### parameters
overlaps <- c(2,4,6,8,10)
disease_effects <- c(1.02,1.04,1.06)
norm_methods <- readRDS("res_summ/norm_methods.rds")
pred_method <- "rfr"

### prediction results
scenario3_df <- data.frame(ed=numeric(0),overlap=character(0),norm_method=character(0),pred_method=character(0),auc=numeric(0))
scenario3_summ_df <- data.frame(ed=numeric(0),overlap=character(0),norm_method=character(0),pred_method=character(0),average_auc=numeric(0))
for(overlap in overlaps){
  for(disease_effect in disease_effects){
    for(norm_method in norm_methods$norm_method){
      df <- readRDS(paste0("scenario3/pred_results/",norm_method,"/pred_ed",disease_effect,"_overlap",overlap,"_",norm_method,"_",pred_method,".rds"))
      summ_df <- data.frame(ed=disease_effect,overlap=overlap,norm_method=norm_method,pred_method=pred_method,
                            average_auc=round(mean(as.numeric(df$auc_value)),3))
      scenario3_df <- rbind(scenario3_df,df)
      scenario3_summ_df <- rbind(scenario3_summ_df,summ_df)
    }
  }
}
saveRDS(scenario3_df,"res_summ/scenario3_res.rds")
saveRDS(scenario3_summ_df,"res_summ/scenario3_res_summ.rds")


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




