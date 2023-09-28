######################################################################################################################
### figures
### 2023/09/28
######################################################################################################################

setwd("/home/wangbb/normalization_comparison")

### packages
all(sapply(c("vegan","ape","reshape2","ggplot2","ggtree","viridis","RColorBrewer","ggsci","aplot","ggpubr",
             "grid","ggplotify","patchwork"),require, character.only=TRUE))

#======================================================================================================================#
### function for merging two count table
merge.func <- function(table1,table2){
  table <- merge(table1,table2,by="row.names",all=T)
  rownames(table) <- table$Row.names
  table <- table[,-grep("Row.names",colnames(table))]
  table[is.na(table)] <- 0
  return(table)
}


#=====================================================================================================================#
### basic statistical analysis for CRC data ###
#=====================================================================================================================#
### data
crc_meta <- readRDS("data/crc_meta.rds")
crc_count <- readRDS("data/crc_count.rds")
crc_meta$study_name_simplified <- factor(crc_meta$study_name_simplified,levels=c("Feng","Gupta","Thomas","Vogtmann","Wirbel","Yachida","Yu","Zeller"))

### PCoA plot
crc_count_norm <- as.data.frame(apply(crc_count,2,function(x) x/sum(x)))
crc_distance <- vegdist(t(crc_count_norm),method="bray")
crc_pcoa <- pcoa(crc_distance)
crc_pcoa_df <- as.data.frame(crc_pcoa$vectors[,1:2])
crc_pcoa_df$sample_id<- rownames(crc_pcoa_df)
crc_pcoa_df <- merge(crc_pcoa_df,crc_meta,by="sample_id")
# axis
crc_pro1 <- as.numeric(sprintf("%.3f",crc_pcoa$values[,"Relative_eig"][1]))*100  
crc_pro2 <- as.numeric(sprintf("%.3f",crc_pcoa$values[,"Relative_eig"][2]))*100
# plot
crc_pcoa_scatter <- ggscatter(crc_pcoa_df,x="Axis.1",y="Axis.2",color="study_name_simplified",fill="study_name_simplified",
                              ellipse=TRUE,mean.point=FALSE,star.plot=TRUE)+
  labs(x=paste0("PCoA1(",crc_pro1,"%)"), y=paste0("PCoA2(",crc_pro2,"%)"),title="PCoA")+
  scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")+
  theme(legend.justification=c(0.02,0.02), legend.position=c(0.02,0.02),#legend.position="right",
        legend.background=element_rect(fill=NA),legend.title=element_blank())+
  theme(plot.margin=margin(t=0, r=0, b=0, l=0, unit="lines"))
crc_pcoa_scatter

# permanova analysis
permanova.func <- function(distance,meta_df,variable){
  meta_df <- meta_df[!is.na(meta_df[,variable]),]
  distance <- as.dist(as.matrix(distance)[meta_df$sample_id,meta_df$sample_id])
  formulas <- as.formula(paste0("distance~",variable))
  permanova <- adonis2(formulas,data=meta_df,permutations=1000,parallel=1)
  c(variable=variable,r2=round(permanova[variable,"R2"],3),pvalue=round(permanova[variable,"Pr(>F)"],3))
}
crc_study_permanova <- permanova.func(distance=crc_distance,meta_df=crc_meta,variable="study_name")
crc_study_r2 <- sprintf("italic(R^2) == %.3f",as.numeric(crc_study_permanova["r2"]))
crc_study_pval <- sprintf("italic(p) == %.3f",as.numeric(crc_study_permanova["pvalue"]))
crc_study_labels <- data.frame(r2=crc_study_r2,pvalue=crc_study_pval,stringsAsFactors = FALSE)
crc_pcoa_scatter <- crc_pcoa_scatter+
  geom_text(data=crc_study_labels,mapping=aes(x=0.35,y=-0.4,label=r2),parse=TRUE,inherit.aes=FALSE,hjust=0,size=4)+
  geom_text(data=crc_study_labels,mapping=aes(x=0.35,y=-0.45,label=pvalue),parse=TRUE,inherit.aes=FALSE,hjust=0,size=4)
crc_pcoa_scatter

### heatmaps of bray curtis distance between different studies
crc_studies <- unique(crc_meta$study_name_simplified)   # crc studies
crc_distance <- as.matrix(crc_distance)
crc_distance_table <- matrix(NA,nrow=length(crc_studies),ncol=length(crc_studies),dimnames=list(crc_studies,crc_studies))
for(i in 1:length(crc_studies)){
  for(j in 1:i){
    distance_df <- crc_distance[crc_meta[crc_meta$study_name_simplified==crc_studies[i],"sample_id"],
                                crc_meta[crc_meta$study_name_simplified==crc_studies[j],"sample_id"]]
    crc_distance_table[i,j] <- round(mean(as.vector(distance_df)),3)
  }
}
crc_distance_table <- as.data.frame(crc_distance_table)
crc_distance_table$study1 <- rownames(crc_distance_table)
crc_distance_table_long <- melt(crc_distance_table,id.vars="study1",variable.name="study2",value.name="distance")
crc_distance_heatmap <- ggplot(crc_distance_table_long, aes(x=study2, y=study1, fill=distance))+ 
  geom_tile() +
  scale_fill_distiller(palette="RdBu", direction=-1, na.value="transparent")+
  geom_text(aes(study2, study1, label=distance), color="black", size=4)+
  labs(x="Dataset 1", y="Dataset 2", fill="Average\nBray-Curtis\ndistance")+
  labs(title="Average Bray-Curtis distance")+
  theme_bw()+
  theme(panel.grid.minor=element_blank())+
  theme(panel.border=element_blank())+
  theme(plot.margin=margin(t=0, r=1, b=0, l=0, unit="lines"))
crc_distance_heatmap

### arrange the plots
ggarrange(crc_pcoa_scatter,crc_distance_heatmap,nrow=1,ncol=2,labels=letters[1:2],widths=c(1,1.2))
ggsave("crc_analysis.jpg",width=12,height=5,dpi=300)


#=====================================================================================================================#
### basic statistical analysis for IBD data ###
#=====================================================================================================================#
### data
ibd_meta <- readRDS("data/ibd_meta.rds")
ibd_count <- readRDS("data/ibd_count.rds")
ibd_meta$study_name_simplified <- factor(ibd_meta$study_name_simplified,levels=c("Hall","HMP","Ijaz","Nielsen","Vila"))

### PCoA plot
ibd_count_norm <- as.data.frame(apply(ibd_count,2,function(x) x/sum(x)))
ibd_distance <- vegdist(t(ibd_count_norm),method="bray")
ibd_pcoa <- pcoa(ibd_distance)
ibd_pcoa_df <- as.data.frame(ibd_pcoa$vectors[,1:2])
ibd_pcoa_df$sample_id<- rownames(ibd_pcoa_df)
ibd_pcoa_df <- merge(ibd_pcoa_df,ibd_meta,by="sample_id")
# axis
ibd_pro1 <- as.numeric(sprintf("%.3f",ibd_pcoa$values[,"Relative_eig"][1]))*100  
ibd_pro2 <- as.numeric(sprintf("%.3f",ibd_pcoa$values[,"Relative_eig"][2]))*100
# plot
ibd_pcoa_df$Axis.2 <- -ibd_pcoa_df$Axis.2
ibd_pcoa_scatter <- ggscatter(ibd_pcoa_df,x="Axis.1",y="Axis.2",color="study_name_simplified",fill="study_name_simplified",
                              ellipse=TRUE,mean.point=FALSE,star.plot=TRUE)+
  labs(x=paste0("PCoA1(",ibd_pro1,"%)"), y=paste0("PCoA2(",ibd_pro2,"%)"),title="PCoA")+
  scale_color_brewer(palette="Set2")+scale_fill_brewer(palette="Set2")+
  theme(legend.justification=c(0.02,0.02), legend.position=c(0.02,0.02),#legend.position="right",
        legend.background=element_rect(fill=NA),legend.title=element_blank())+
  theme(plot.margin=margin(t=0, r=0, b=0, l=0, unit="lines"))
ibd_pcoa_scatter

# permanova analysis
ibd_study_permanova <- permanova.func(distance=ibd_distance,meta_df=ibd_meta,variable="study_name")
ibd_study_r2 <- sprintf("italic(R^2) == %.3f",as.numeric(ibd_study_permanova["r2"]))
ibd_study_pval <- sprintf("italic(p) == %.3f",as.numeric(ibd_study_permanova["pvalue"]))
ibd_study_labels <- data.frame(r2=ibd_study_r2,pvalue=ibd_study_pval,stringsAsFactors = FALSE)
ibd_pcoa_scatter <- ibd_pcoa_scatter+
  geom_text(data=ibd_study_labels,mapping=aes(x=0.35,y=-0.5,label=r2),parse=TRUE,inherit.aes=FALSE,hjust=0,size=4)+
  geom_text(data=ibd_study_labels,mapping=aes(x=0.35,y=-0.55,label=pvalue),parse=TRUE,inherit.aes=FALSE,hjust=0,size=4)
ibd_pcoa_scatter

### heatmaps of bray curtis distance between different studies
ibd_studies <- levels(ibd_meta$study_name_simplified)   # ibd studies
ibd_distance <- as.matrix(ibd_distance)
ibd_distance_table <- matrix(NA,nrow=length(ibd_studies),ncol=length(ibd_studies),dimnames=list(ibd_studies,ibd_studies))
for(i in 1:length(ibd_studies)){
  for(j in 1:i){
    distance_df <- ibd_distance[ibd_meta[ibd_meta$study_name_simplified==ibd_studies[i],"sample_id"],
                                ibd_meta[ibd_meta$study_name_simplified==ibd_studies[j],"sample_id"]]
    ibd_distance_table[i,j] <- round(mean(as.vector(distance_df)),3)
  }
}
ibd_distance_table <- as.data.frame(ibd_distance_table)
ibd_distance_table$study1 <- rownames(ibd_distance_table)
ibd_distance_table_long <- melt(ibd_distance_table,id.vars="study1",variable.name="study2",value.name="distance")
ibd_distance_heatmap <- ggplot(ibd_distance_table_long, aes(x=study2, y=study1, fill=distance))+ 
  geom_tile() +
  scale_fill_distiller(palette="RdBu", direction=-1, na.value="transparent")+
  geom_text(aes(study2, study1, label=distance), color="black", size=4)+
  labs(x="Dataset 1", y="Dataset 2", fill="Average\nBray-Curtis\ndistance")+
  labs(title="Average Bray-Curtis distance")+
  theme_bw()+
  theme(panel.grid.minor=element_blank())+
  theme(panel.border=element_blank())+
  theme(plot.margin=margin(t=0, r=1, b=0, l=0, unit="lines"))
ibd_distance_heatmap

### arrange the plots
ggarrange(ibd_pcoa_scatter,ibd_distance_heatmap,nrow=1,ncol=2,labels=letters[1:2],widths=c(1,1.2))
ggsave("ibd_analysis.jpg",width=12,height=5,dpi=300)


#=====================================================================================================================#
### simulation scenario1 ###
#=====================================================================================================================#
### data 
norm_methods <- readRDS("data/norm_methods.rds")
scenario1_summ <- readRDS("data/scenario1_res_summ.rds")
scenario1_summ <- merge(scenario1_summ,norm_methods,by="norm_method")
scenario1_summ$ep <- factor(scenario1_summ$ep,levels=c(0,0.05,0.1,0.15,0.2,0.25,0.4,0.6,0.8,1))

### margins of normalization methods
margins <- ggplot(norm_methods,aes("",annotated_method,fill=class))+
  geom_tile() + 
  scale_fill_jco()+
  labs(fill="Group of\nnormalization\nmethods")+
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.y=element_blank())
margins

### heatmap
heatmap.scenario1.func <- function(df_data,ed,title){
  # heatmaps
  df_data <- df_data[df_data$ed==ed,]
  df_data <- df_data[df_data$norm_method%in%norm_methods$norm_method,]
  df_data$average_auc <- round(df_data$average_auc,3)
  heatmaps <- ggplot(df_data, aes(x=ep, y=annotated_method, fill=average_auc))+ 
    geom_tile() +
    geom_text(aes(x=ep, y=annotated_method, label=average_auc), color="black", size=4)+
    labs(x="population effect",y="",fill="Average\nAUC",title=title)+
    scale_fill_distiller(palette="Reds",limits=c(0.4,1),direction=1)+
    theme_bw()+
    theme(panel.grid.minor=element_blank(), panel.border=element_blank())
  heatmaps
}
p1 <- heatmap.scenario1.func(df_data=scenario1_summ,ed=1.02,title=expression(paste(bold("a")," ed = 1.02")))
p2 <- heatmap.scenario1.func(df_data=scenario1_summ,ed=1.04,title=expression(paste(bold("b")," ed = 1.04")))+theme(axis.text.y=element_blank())
p3 <- heatmap.scenario1.func(df_data=scenario1_summ,ed=1.06,title=expression(paste(bold("c")," ed = 1.06")))+theme(axis.text.y=element_blank())
p1 %>% insert_left(margins,width=0.05) %>% 
  insert_right(p2, width=1) %>%
  insert_right(p3,width=1)
ggsave("sim_scenario1.jpg",width=18,height=6,dpi=300)


#=====================================================================================================================#
### simulation scenario2 ###
#=====================================================================================================================#
### data 
norm_methods <- readRDS("data/norm_methods.rds")
scenario2_summ <- readRDS("data/scenario2_res_summ.rds")
scenario2_summ <- merge(scenario2_summ,norm_methods,by="norm_method")
scenario2_summ$factor <- factor(paste0("m",scenario2_summ$batch_mean,"v",scenario2_summ$batch_var),
                                levels=c("m0v1","m0v2","m0v4","m500v1","m1000v1"))

### margins of normalization methods
margins <- ggplot(norm_methods,aes("",annotated_method,fill=class))+
  geom_tile() + 
  scale_fill_jco()+
  labs(fill="Group of\nnormalization\nmethods")+
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.y=element_blank())
margins

### heatmap
heatmap.scenario2.func <- function(df_data,ed,title){
  # heatmaps
  df_data <- df_data[df_data$ed==ed,]
  df_data <- df_data[df_data$norm_method%in%norm_methods$norm_method,]
  df_data$average_auc <- round(df_data$average_auc,3)
  heatmaps <- ggplot(df_data, aes(x=factor, y=annotated_method, fill=average_auc))+ 
    geom_tile() +
    geom_text(aes(x=factor, y=annotated_method, label=average_auc), color="black", size=4)+
    labs(x="factor",y="",fill="Average\nAUC",title=title)+
    scale_fill_distiller(palette="Reds",limits=c(0.4,1),direction=1)+
    theme_bw()+
    theme(panel.grid.minor=element_blank(), panel.border=element_blank())
  heatmaps
}
p1 <- heatmap.scenario2.func(df_data=scenario2_summ,ed=1.02,title=expression(paste(bold("a")," ed = 1.02")))
p2 <- heatmap.scenario2.func(df_data=scenario2_summ,ed=1.04,title=expression(paste(bold("b")," ed = 1.04")))+theme(axis.text.y=element_blank())
p3 <- heatmap.scenario2.func(df_data=scenario2_summ,ed=1.06,title=expression(paste(bold("c")," ed = 1.06")))+theme(axis.text.y=element_blank())
p1 %>% insert_left(margins,width=0.05) %>% 
  insert_right(p2, width=1) %>%
  insert_right(p3,width=1)
ggsave("sim_scenario2.jpg",width=18,height=6,dpi=300)


#=====================================================================================================================#
### simulation scenario3 ###
#=====================================================================================================================#
### data 
norm_methods <- readRDS("data/norm_methods.rds")
scenario3_summ <- readRDS("data/scenario3_res_summ.rds")
scenario3_summ <- merge(scenario3_summ,norm_methods,by="norm_method")
scenario3_summ$overlap <- factor(scenario3_summ$overlap,levels=c(2,4,6,8,10))

### margins of normalization methods
margins <- ggplot(norm_methods,aes("",annotated_method,fill=class))+
  geom_tile() + 
  scale_fill_jco()+
  labs(fill="Group of\nnormalization\nmethods")+
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.y=element_blank())
margins

### heatmap
heatmap.scenario3.func <- function(df_data,ed,title){
  # heatmaps
  df_data <- df_data[df_data$ed==ed,]
  df_data <- df_data[df_data$norm_method%in%norm_methods$norm_method,]
  df_data$average_auc <- round(df_data$average_auc,3)
  heatmaps <- ggplot(df_data, aes(x=overlap, y=annotated_method, fill=average_auc))+ 
    geom_tile() +
    geom_text(aes(x=overlap, y=annotated_method, label=average_auc), color="black", size=4)+
    labs(x="overlap",y="",fill="Average\nAUC",title=title)+
    scale_fill_distiller(palette="Reds",limits=c(0.56,1),direction=1)+
    theme_bw()+
    theme(panel.grid.minor=element_blank(), panel.border=element_blank())
  heatmaps
}
p1 <- heatmap.scenario3.func(df_data=scenario3_summ,ed=1.02,title=expression(paste(bold("a")," ed = 1.02")))
p2 <- heatmap.scenario3.func(df_data=scenario3_summ,ed=1.04,title=expression(paste(bold("b")," ed = 1.04")))+theme(axis.text.y=element_blank())
p3 <- heatmap.scenario3.func(df_data=scenario3_summ,ed=1.06,title=expression(paste(bold("c")," ed = 1.06")))+theme(axis.text.y=element_blank())
p1 %>% insert_left(margins,width=0.05) %>% 
  insert_right(p2, width=1) %>%
  insert_right(p3,width=1)
ggsave("sim_scenario3.jpg",width=18,height=6,dpi=300)


#=====================================================================================================================#
### CRC data ###
#=====================================================================================================================#
### data 
crc_pred <- readRDS("data/crc_pred.rds")
crc_pred_summ <- readRDS("data/crc_pred_summ.rds")
norm_methods <- readRDS("data/norm_methods.rds")

### change the study name to the simplified version
crc_pred$study1 <- gsub("([A-Z][a-z]+)[A-Z]+_[0-9]+", "\\1", crc_pred$study1)
crc_pred$study2 <- gsub("([A-Z][a-z]+)[A-Z]+_[0-9]+", "\\1", crc_pred$study2)
crc_pred_summ$study1 <- gsub("([A-Z][a-z]+)[A-Z]+_[0-9]+", "\\1", crc_pred_summ$study1)
crc_pred_summ$study2 <- gsub("([A-Z][a-z]+)[A-Z]+_[0-9]+", "\\1", crc_pred_summ$study2)

### combine the pred results and norm methods
crc_pred <- merge(crc_pred,norm_methods,by="norm_method")
crc_pred_summ <- merge(crc_pred_summ,norm_methods,by="norm_method")

### box plot for auc values using different normalization method of different training and testing set
# add the trn and tst information
crc_pred$study1 <- paste("trn:",crc_pred$study1)
crc_pred$study2 <- paste("tst:",crc_pred$study2)
# plot 
ggboxplot(crc_pred,x="annotated_method",y="auc_values",color="class")+
  geom_hline(aes(yintercept=0.4), colour="grey", linetype="dashed")+
  geom_hline(aes(yintercept=0.6), colour="grey", linetype="dashed")+
  geom_hline(aes(yintercept=0.8), colour="grey", linetype="dashed")+
  scale_color_jco()+
  scale_y_continuous(limits=c(0.2,1), breaks=seq(0.4,1,0.2))+
  labs(x="Normalization methods",y="AUC", color="Group of normalization methods")+
  theme(panel.border=element_rect(fill=NA,color="black"))+
  theme(axis.text.y=element_text(size=8))+
  facet_grid(study1~study2,scales="free")+
  coord_flip()
ggsave("crc_pred_box.jpg",width=16,height=18,dpi=300)


### figure in the main text, 
#   including the boxplot for auc values using the simulation template data, 
#   and the boxplot for the rank of different normalization methods for CRC datasets
# boxplot for auc values using trn: Gupta and tst: Feng
box_auc_template <- ggboxplot(crc_pred[crc_pred$study1=="trn: Gupta" & crc_pred$study2=="tst: Feng",],
                              x="annotated_method",y="auc_values",color="class")+
  geom_hline(aes(yintercept=0.5), colour="grey", linetype="dashed")+
  scale_color_jco()+
  labs(x="Normalization methods",y="AUC",color="Group of normalization methods",
       title="AUC values of model trained on Gupta and tested on Feng")+
  theme(axis.text.x=element_text(vjust=0.5,hjust=1,angle=90))+
  theme(legend.position="bottom")
box_auc_template

# boxplot for the rank of different normalization methods
box_crc_rank <- ggboxplot(crc_pred_summ,x="annotated_method",y="order_average_auc",color="class")+
  scale_color_jco()+
  labs(x="Normalization methods",y="Ranks",color="Group of normalization methods",
       title="Ranks of normalizaton methods for CRC prediction")+
  theme(axis.text.x=element_text(vjust=0.5,hjust=1,angle=90))+
  theme(legend.position="bottom")
box_crc_rank

# combine box_auc_template and box_crc_rank
ggarrange(box_auc_template,box_crc_rank,ncol=1,labels=letters[1:2],common.legend=T)
ggsave("crc_pred_res.jpg",width=10,height=5,dpi=300)


#=====================================================================================================================#
### prediction results for IBD datasets ###
#=====================================================================================================================#
### data
ibd_pred <- readRDS("data/ibd_pred.rds")
ibd_pred_summ <- readRDS("data/ibd_pred_summ.rds")
norm_methods <- readRDS("data/norm_methods.rds")

### change the study name to the simplified version
ibd_pred$study1 <- gsub("([A-Z][a-z]+)[A-Z]+_[0-9]+", "\\1", ibd_pred$study1)
ibd_pred$study2 <- gsub("([A-Z][a-z]+)[A-Z]+_[0-9]+", "\\1", ibd_pred$study2)
ibd_pred_summ$study1 <- gsub("([A-Z][a-z]+)[A-Z]+_[0-9]+", "\\1", ibd_pred_summ$study1)
ibd_pred_summ$study2 <- gsub("([A-Z][a-z]+)[A-Z]+_[0-9]+", "\\1", ibd_pred_summ$study2)
# modify HMP that in a different pattern separately
ibd_pred$study1 <- sub("HMP_2019_ibdmdb","HMP",ibd_pred$study1)
ibd_pred$study2 <- sub("HMP_2019_ibdmdb","HMP",ibd_pred$study2)
ibd_pred_summ$study1 <- sub("HMP_2019_ibdmdb","HMP",ibd_pred_summ$study1)
ibd_pred_summ$study2 <- sub("HMP_2019_ibdmdb","HMP",ibd_pred_summ$study2)

### combine the pred results and norm methods
ibd_pred <- merge(ibd_pred,norm_methods,by="norm_method")
ibd_pred_summ <- merge(ibd_pred_summ,norm_methods,by="norm_method")

### box plot for auc values using different normalization method of different training and testing set
# add the trn and tst information
ibd_pred$study1 <- paste("trn:",ibd_pred$study1)
ibd_pred$study2 <- paste("tst:",ibd_pred$study2)
# plot 
ggboxplot(ibd_pred,x="annotated_method",y="auc_values",color="class")+
  geom_hline(aes(yintercept=0.5), colour="grey", linetype="dashed")+
  scale_color_jco()+
  scale_y_continuous(limits=c(0.2,1), breaks=seq(0.4,1,0.2))+
  labs(x="Normalization methods",y="AUC", color="Group of normalization methods")+
  theme(panel.border=element_rect(fill=NA,color="black"))+
  theme(axis.text.y=element_text(size=8))+
  facet_grid(study1~study2,scales="free")+
  coord_flip()
ggsave("ibd_pred_box.jpg",width=10,height=12,dpi=300)

### boxplot for the rank of different normalization methods
ggboxplot(ibd_pred_summ,x="annotated_method",y="order_average_auc",color="class")+
  scale_color_jco()+
  labs(x="Normalization methods",y="Ranks",color="Group of normalization methods")+
  #title="Ranks of normalizaton methods for IBD prediction")+
  theme(axis.text.x=element_text(vjust=0.5,hjust=1,angle=90))+
  theme(legend.position="top")
ggsave("ibd_pred_rank.jpg",width=10,height=3,dpi=300)

















