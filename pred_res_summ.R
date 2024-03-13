######################################################################################################################
### summarize the results
### 2024/03/12
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
saveRDS(norm_methods,"data/norm_methods.rds")


#=====================================================================================================================#
### simulation scenario1 ###
#=====================================================================================================================#
### parameters
population_effects <- c(0,0.05,0.1,0.15,0.2,0.25,0.4,0.6,0.8,1)  # population effect
disease_effects <- c(1.02,1.04,1.06)                             # disease effect
norm_methods <- readRDS("data/norm_methods.rds")
pred_method <- "rfr"

### prediction results
scenario1_df <- data.frame(ep=numeric(0),ed=numeric(0),norm_method=character(0),pred_method=character(0),
                           auc=numeric(0),accuracy=numeric(0),specificity=numeric(0),sensitivity=numeric(0))
scenario1_summ_df <- data.frame(ep=numeric(0),ed=numeric(0),norm_method=character(0),pred_method=character(0),
                                average_auc=numeric(0),average_accuracy=numeric(0),average_specificity=numeric(0),average_sensitivity=numeric(0))
for(alpha in population_effects){
  for(disease_effect in disease_effects){
    for(norm_method in norm_methods$norm_method){
      #alpha=0;disease_effect=1.02;norm_method="TSS"
      df <- readRDS(paste0("scenario1/pred_results/",norm_method,"/norm_alpha",alpha,"_ed",disease_effect,"_",norm_method,"_",pred_method,".rds"))
      summ_df <- data.frame(ep=alpha,ed=disease_effect,norm_method=norm_method,pred_method=pred_method,
                            average_auc=round(mean(as.numeric(df$auc_value)),3),
                            average_accuracy=round(mean(as.numeric(df$accuracy)),3),
                            average_specificity=round(mean(as.numeric(df$specificity)),3),
                            average_sensitivity=round(mean(as.numeric(df$sensitivity)),3))
      scenario1_df <- rbind(scenario1_df,df)
      scenario1_summ_df <- rbind(scenario1_summ_df,summ_df)
    }
  }
}
saveRDS(scenario1_df,"data/scenario1_res.rds")
saveRDS(scenario1_summ_df,"data/scenario1_res_summ.rds")


#=====================================================================================================================#
### simulation scenario2 ###
#=====================================================================================================================#
### parameters
parameters <- list(c(0,1),c(500,1),c(1000,1),c(0,2),c(0,4))
disease_effects <- c(1.02,1.04,1.06)
norm_methods <- readRDS("data/norm_methods.rds")
pred_method <- "rfr"

### prediction results
scenario2_df <- data.frame(ed=numeric(0),batch_mean=numeric(0),batch_var=numeric(0),norm_method=character(0),pred_method=character(0),
                           auc=numeric(0),accuracy=numeric(0),specificity=numeric(0),sensitivity=numeric(0))
scenario2_summ_df <- data.frame(ed=numeric(0),batch_mean=numeric(0),batch_var=numeric(0),norm_method=character(0),pred_method=character(0),
                                average_auc=numeric(0),average_accuracy=numeric(0),average_specificity=numeric(0),average_sensitivity=numeric(0))
for(i in 1:length(parameters)){
  for(disease_effect in disease_effects){
    for(norm_method in norm_methods$norm_method){
      #alpha=0;disease_effect=1.02;norm_method="TSS"
      df <- readRDS(paste0("scenario2/pred_results/",norm_method,"/pred_ed",disease_effect,"_mean",parameters[[i]][1],"_var",parameters[[i]][2],"_",norm_method,"_",pred_method,".rds"))
      summ_df <- data.frame(ed=disease_effect,batch_mean=parameters[[i]][1],batch_var=parameters[[i]][2],
                            norm_method=norm_method,pred_method=pred_method,
                            average_auc=round(mean(as.numeric(df$auc_value)),3),
                            average_accuracy=round(mean(as.numeric(df$accuracy)),3),
                            average_specificity=round(mean(as.numeric(df$specificity)),3),
                            average_sensitivity=round(mean(as.numeric(df$sensitivity)),3))
      scenario2_df <- rbind(scenario2_df,df)
      scenario2_summ_df <- rbind(scenario2_summ_df,summ_df)
    }
  }
}
saveRDS(scenario2_df,"data/scenario2_res.rds")
saveRDS(scenario2_summ_df,"data/scenario2_res_summ.rds")


#=====================================================================================================================#
### simulation scenario3 ###
#=====================================================================================================================#
### parameters
overlaps <- c(2,4,6,8,10)
disease_effects <- c(1.02,1.04,1.06)
norm_methods <- readRDS("data/norm_methods.rds")
pred_method <- "rfr"

### prediction results
scenario3_df <- data.frame(ed=numeric(0),overlap=character(0),norm_method=character(0),pred_method=character(0),
                           auc=numeric(0),accuracy=numeric(0),specificity=numeric(0),sensitivity=numeric(0))
scenario3_summ_df <- data.frame(ed=numeric(0),overlap=character(0),norm_method=character(0),pred_method=character(0),
                                average_auc=numeric(0),average_accuracy=numeric(0),average_specificity=numeric(0),average_sensitivity=numeric(0))
for(overlap in overlaps){
  for(disease_effect in disease_effects){
    for(norm_method in norm_methods$norm_method){
      df <- readRDS(paste0("scenario3/pred_results/",norm_method,"/pred_ed",disease_effect,"_overlap",overlap,"_",norm_method,"_",pred_method,".rds"))
      summ_df <- data.frame(ed=disease_effect,overlap=overlap,norm_method=norm_method,pred_method=pred_method,
                            average_auc=round(mean(as.numeric(df$auc_value)),3),
                            average_accuracy=round(mean(as.numeric(df$accuracy)),3),
                            average_specificity=round(mean(as.numeric(df$specificity)),3),
                            average_sensitivity=round(mean(as.numeric(df$sensitivity)),3))
      scenario3_df <- rbind(scenario3_df,df)
      scenario3_summ_df <- rbind(scenario3_summ_df,summ_df)
    }
  }
}
saveRDS(scenario3_df,"data/scenario3_res.rds")
saveRDS(scenario3_summ_df,"data/scenario3_res_summ.rds")


#=====================================================================================================================#
### prediction results for CRC data ###
#=====================================================================================================================#
### related information for CRC datasets 
crc_meta <- readRDS("data/crc_meta.rds")
crc_studies <- unique(crc_meta$study_name)
norm_methods <- readRDS("data/norm_methods.rds")

### prediction results of cross-study prediction using ranger
crc_df <- data.frame(study1=character(0),study2=character(0),norm_method=character(0),pred_method=character(0),
                     auc=numeric(0),accuracy=numeric(0),specificity=numeric(0),sensitivity=numeric(0))
crc_summ_df <- data.frame(study1=character(0),study2=character(0),norm_method=character(0),pred_method=character(0),
                          average_auc=numeric(0),average_accuracy=numeric(0),average_specificity=numeric(0),average_sensitivity=numeric(0))
for(study1 in crc_studies){
  for(study2 in crc_studies){
    for(norm_method in norm_methods$norm_method){
      if(study1 != study2){
        # study1="FengQ_2015";study2="GuptaA_2019";norm_method="TSS";pred_method="ranger"
		df <- readRDS(paste0("./real_data_pred/",pred_method,"_",norm_method,"_trn_",study1,"_tst_",study2,".rds"))
        df$accuracy <- 1-df$misclassification_rate
        summ_df <- data.frame(study1=study1,study2=study2,norm_method=norm_method,pred_method=pred_method,
                              average_auc=round(mean(as.numeric(df$auc_value)),3),
                              average_accuracy=round(mean(as.numeric(df$accuracy)),3),
                              average_specificity=round(mean(as.numeric(df$specificity)),3),
                              average_sensitivity=round(mean(as.numeric(df$sensitivity)),3))
        crc_df <- rbind(crc_df,df)
        crc_summ_df <- rbind(crc_summ_df,summ_df)
      }
    }
  }
}
saveRDS(crc_df,"data/crc_pred.rds")
saveRDS(crc_summ_df,"data/crc_pred_summ.rds")

### add the ranks according to the average auc, accuracy, specificity, sensitivity
crc_summ_final <- data.frame(study1=character(0),study2=character(0),norm_method=character(0),pred_method=character(0),
                             average_auc=numeric(0),average_accuracy=numeric(0),average_specificity=numeric(0),average_sensitivity=numeric(0),
                             order_auc=numeric(0),order_accuracy=numeric(0),order_specificity=numeric(0),order_sensitivity=numeric(0))
for(study1 in crc_studies){
  for(study2 in crc_studies){
    if(study1 != study2){
      df <- crc_summ_df[crc_summ_df$study1==study1 & crc_summ_df$study2==study2,]
      # orders according to the descending order of average auc
      df <- df[order(df$average_auc,decreasing=TRUE),]
      df$order_auc <- 1:nrow(df)
      # orders according to the descending order of average accuracy
      df <- df[order(df$average_accuracy,decreasing=TRUE),]
      df$order_accuracy <- 1:nrow(df)
      # orders according to the descending order of average specificity
      df <- df[order(df$average_specificity,decreasing=TRUE),]
      df$order_specificity <- 1:nrow(df)
      # orders according to the descending order of average sensitivity
      df <- df[order(df$average_sensitivity,decreasing=TRUE),]
      df$order_sensitivity <- 1:nrow(df)
      # combine the results
      crc_summ_final <- rbind(crc_summ_final,df)
    }
  }
}
saveRDS(crc_summ_final,"data/crc_pred_summ.rds")


#=====================================================================================================================#
### prediction results for ibd data ###
#=====================================================================================================================#
### related information for ibd datasets 
ibd_meta <- readRDS("data/ibd_meta.rds")
ibd_studies <- unique(ibd_meta$study_name)
norm_methods <- readRDS("data/norm_methods.rds")

### results of cross-study prediction using ranger
ibd_df <- data.frame(study1=character(0),study2=character(0),norm_method=character(0),pred_method=character(0),
                     auc=numeric(0),accuracy=numeric(0),specificity=numeric(0),sensitivity=numeric(0))
ibd_summ_df <- data.frame(study1=character(0),study2=character(0),norm_method=character(0),pred_method=character(0),
                          average_auc=numeric(0),average_accuracy=numeric(0), average_specificity=numeric(0),average_sensitivity=numeric(0))
for(study1 in ibd_studies){
  for(study2 in ibd_studies){
    for(norm_method in norm_methods$norm_method){
      if(study1 != study2){
        # study1="FengQ_2015";study2="GuptaA_2019";norm_method="TSS";pred_method="ranger"
		df <- readRDS(paste0("./real_data_pred/",pred_method,"_",norm_method,"_trn_",study1,"_tst_",study2,".rds"))
        summ_df <- data.frame(study1=study1,study2=study2,norm_method=norm_method,pred_method=pred_method,
                              average_auc=round(mean(as.numeric(df$auc_value)),3),
                              average_accuracy=round(mean(as.numeric(df$accuracy)),3),
                              average_specificity=round(mean(as.numeric(df$specificity)),3),
                              average_sensitivity=round(mean(as.numeric(df$sensitivity)),3))
        ibd_df <- rbind(ibd_df,df)
        ibd_summ_df <- rbind(ibd_summ_df,summ_df)
      }
    }
  }
}
saveRDS(ibd_df,"data/ibd_pred.rds")
saveRDS(ibd_summ_df,"data/ibd_pred_summ.rds")


### add the ranks according to the average auc, misclassification rate, specificity, sensitivity
ibd_summ_final <- data.frame(study1=character(0),study2=character(0),norm_method=character(0),pred_method=character(0),
                             average_auc=numeric(0),average_accuracy=numeric(0),average_specificity=numeric(0),average_sensitivity=numeric(0),
                             order_auc=numeric(0),order_accuracy=numeric(0),order_specificity=numeric(0),order_sensitivity=numeric(0))
for(study1 in ibd_studies){
  for(study2 in ibd_studies){
    if(study1 != study2){
      df <- ibd_summ_df[ibd_summ_df$study1==study1 & ibd_summ_df$study2==study2,]
      # orders according to the descending order of average auc
      df <- df[order(df$average_auc,decreasing=TRUE),]
      df$order_auc <- 1:nrow(df)
      # orders according to the descending order of average accuracy
      df <- df[order(df$average_accuracy,decreasing=TRUE),]
      df$order_accuracy <- 1:nrow(df)
      # orders according to the descending order of average specificity
      df <- df[order(df$average_specificity,decreasing=TRUE),]
      df$order_specificity <- 1:nrow(df)
      # orders according to the descending order of average sensitivity
      df <- df[order(df$average_sensitivity,decreasing=TRUE),]
      df$order_sensitivity <- 1:nrow(df)
      # combine the results
      ibd_summ_final <- rbind(ibd_summ_final,df)
    }
  }
}
saveRDS(ibd_summ_final,"data/ibd_pred_summ.rds")




