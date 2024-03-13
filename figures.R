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


### heatmap for auc 
heatmap.auc.scenario1.func <- function(df_data,ed,title){
  # heatmaps
  df_data <- df_data[df_data$ed==ed,]
  df_data <- df_data[df_data$norm_method%in%norm_methods$norm_method,]
  df_data$average_auc <- round(df_data$average_auc,3)
  heatmaps <- ggplot(df_data, aes(x=ep, y=annotated_method, fill=average_auc))+ 
    geom_tile() +
    geom_text(aes(x=ep, y=annotated_method, label=average_auc), color="black", size=4)+
    labs(x="population effect",y="",fill="Average\nAUC",title=title)+
    scale_fill_distiller(palette="Reds",limits=c(0.48,1),direction=1)+
    theme_bw()+
    theme(panel.grid.minor=element_blank(), panel.border=element_blank())
  heatmaps
}
p1 <- heatmap.auc.scenario1.func(df_data=scenario1_summ,ed=1.02,title=expression(paste(bold("a")," ed = 1.02")))
p2 <- heatmap.auc.scenario1.func(df_data=scenario1_summ,ed=1.04,title=expression(paste(bold("b")," ed = 1.04")))+theme(axis.text.y=element_blank())
p3 <- heatmap.auc.scenario1.func(df_data=scenario1_summ,ed=1.06,title=expression(paste(bold("c")," ed = 1.06")))+theme(axis.text.y=element_blank())
p1 %>% insert_left(margins,width=0.05) %>% 
  insert_right(p2, width=1) %>%
  insert_right(p3,width=1)
ggsave("figures/sim_scenario1_auc.jpg",width=18,height=6,dpi=300)


### heatmap for accuracy
heatmap.accuracy.scenario1.func <- function(df_data,ed,title){
  # heatmaps
  df_data <- df_data[df_data$ed==ed,]
  df_data <- df_data[df_data$norm_method%in%norm_methods$norm_method,]
  df_data$average_accuracy <- round(df_data$average_accuracy,3)
  heatmaps <- ggplot(df_data, aes(x=ep, y=annotated_method, fill=average_accuracy))+ 
    geom_tile() +
    geom_text(aes(x=ep, y=annotated_method, label=average_accuracy), color="black", size=4)+
    labs(x="population effect",y="",fill="Average\naccuracy",title=title)+
    scale_fill_distiller(palette="Purples",limits=c(0.49,1),direction=1)+
    theme_bw()+
    theme(panel.grid.minor=element_blank(), panel.border=element_blank())
  heatmaps
}
p1 <- heatmap.accuracy.scenario1.func(df_data=scenario1_summ,ed=1.02,title=expression(paste(bold("a")," ed = 1.02")))
p2 <- heatmap.accuracy.scenario1.func(df_data=scenario1_summ,ed=1.04,title=expression(paste(bold("b")," ed = 1.04")))+theme(axis.text.y=element_blank())
p3 <- heatmap.accuracy.scenario1.func(df_data=scenario1_summ,ed=1.06,title=expression(paste(bold("c")," ed = 1.06")))+theme(axis.text.y=element_blank())
p1 %>% insert_left(margins,width=0.05) %>% 
  insert_right(p2, width=1) %>%
  insert_right(p3,width=1)
ggsave("figures/sim_scenario1_accuracy.jpg",width=18,height=6,dpi=300)


### heatmap for specificity
heatmap.specificity.scenario1.func <- function(df_data,ed,title){
  # heatmaps
  df_data <- df_data[df_data$ed==ed,]
  df_data <- df_data[df_data$norm_method%in%norm_methods$norm_method,]
  df_data$average_specificity <- round(df_data$average_specificity,3)
  heatmaps <- ggplot(df_data, aes(x=ep, y=annotated_method, fill=average_specificity))+ 
    geom_tile() +
    geom_text(aes(x=ep, y=annotated_method, label=average_specificity), color="black", size=4)+
    labs(x="population effect",y="",fill="Average\nspecificity",title=title)+
    scale_fill_distiller(palette="Greens",limits=c(0,1),direction=1)+
    theme_bw()+
    theme(panel.grid.minor=element_blank(), panel.border=element_blank())
  heatmaps
}
p1 <- heatmap.specificity.scenario1.func(df_data=scenario1_summ,ed=1.02,title=expression(paste(bold("a")," ed = 1.02")))
p2 <- heatmap.specificity.scenario1.func(df_data=scenario1_summ,ed=1.04,title=expression(paste(bold("b")," ed = 1.04")))+theme(axis.text.y=element_blank())
p3 <- heatmap.specificity.scenario1.func(df_data=scenario1_summ,ed=1.06,title=expression(paste(bold("c")," ed = 1.06")))+theme(axis.text.y=element_blank())
p1 %>% insert_left(margins,width=0.05) %>% 
  insert_right(p2, width=1) %>%
  insert_right(p3,width=1)
ggsave("figures/sim_scenario1_specificity.jpg",width=18,height=6,dpi=300)


### heatmap for sensitivity
heatmap.sensitivity.scenario1.func <- function(df_data,ed,title){
  # heatmaps
  df_data <- df_data[df_data$ed==ed,]
  df_data <- df_data[df_data$norm_method%in%norm_methods$norm_method,]
  df_data$average_sensitivity <- round(df_data$average_sensitivity,3)
  heatmaps <- ggplot(df_data, aes(x=ep, y=annotated_method, fill=average_sensitivity))+ 
    geom_tile() +
    geom_text(aes(x=ep, y=annotated_method, label=average_sensitivity), color="black", size=4)+
    labs(x="population effect",y="",fill="Average\nsensitivity",title=title)+
    scale_fill_distiller(palette="Blues",limits=c(0,1),direction=1)+
    theme_bw()+
    theme(panel.grid.minor=element_blank(), panel.border=element_blank())
  heatmaps
}
p1 <- heatmap.sensitivity.scenario1.func(df_data=scenario1_summ,ed=1.02,title=expression(paste(bold("a")," ed = 1.02")))
p2 <- heatmap.sensitivity.scenario1.func(df_data=scenario1_summ,ed=1.04,title=expression(paste(bold("b")," ed = 1.04")))+theme(axis.text.y=element_blank())
p3 <- heatmap.sensitivity.scenario1.func(df_data=scenario1_summ,ed=1.06,title=expression(paste(bold("c")," ed = 1.06")))+theme(axis.text.y=element_blank())
p1 %>% insert_left(margins,width=0.05) %>% 
  insert_right(p2, width=1) %>%
  insert_right(p3,width=1)
ggsave("figures/sim_scenario1_sensitivity.jpg",width=18,height=6,dpi=300)
                                   



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

                                      
### heatmap for auc
heatmap.auc.scenario2.func <- function(df_data,ed,title){
  # heatmaps
  df_data <- df_data[df_data$ed==ed,]
  df_data <- df_data[df_data$norm_method%in%norm_methods$norm_method,]
  df_data$average_auc <- round(df_data$average_auc,3)
  heatmaps <- ggplot(df_data, aes(x=factor, y=annotated_method, fill=average_auc))+ 
    geom_tile() +
    geom_text(aes(x=factor, y=annotated_method, label=average_auc), color="black", size=4)+
    labs(x="factor",y="",fill="Average\nAUC",title=title)+
    scale_fill_distiller(palette="Reds",limits=c(0.47,1),direction=1)+
    theme_bw()+
    theme(panel.grid.minor=element_blank(), panel.border=element_blank())
  heatmaps
}
p1 <- heatmap.auc.scenario2.func(df_data=scenario2_summ,ed=1.02,title=expression(paste(bold("a")," ed = 1.02")))
p2 <- heatmap.auc.scenario2.func(df_data=scenario2_summ,ed=1.04,title=expression(paste(bold("b")," ed = 1.04")))+theme(axis.text.y=element_blank())
p3 <- heatmap.auc.scenario2.func(df_data=scenario2_summ,ed=1.06,title=expression(paste(bold("c")," ed = 1.06")))+theme(axis.text.y=element_blank())
p1 %>% insert_left(margins,width=0.05) %>% 
  insert_right(p2, width=1) %>%
  insert_right(p3,width=1)
ggsave("figures/sim_scenario2_auc.jpg",width=18,height=6,dpi=300)


### heatmap for accuracy
heatmap.accuracy.scenario2.func <- function(df_data,ed,title){
  # heatmaps
  df_data <- df_data[df_data$ed==ed,]
  df_data <- df_data[df_data$norm_method%in%norm_methods$norm_method,]
  df_data$average_accuracy <- round(df_data$average_accuracy,3)
  heatmaps <- ggplot(df_data, aes(x=factor, y=annotated_method, fill=average_accuracy))+ 
    geom_tile() +
    geom_text(aes(x=factor, y=annotated_method, label=average_accuracy), color="black", size=4)+
    labs(x="factor",y="",fill="Average\naccuracy",title=title)+
    scale_fill_distiller(palette="Purples",limits=c(0.47,0.88),direction=1)+
    theme_bw()+
    theme(panel.grid.minor=element_blank(), panel.border=element_blank())
  heatmaps
}
p1 <- heatmap.accuracy.scenario2.func(df_data=scenario2_summ,ed=1.02,title=expression(paste(bold("a")," ed = 1.02")))
p2 <- heatmap.accuracy.scenario2.func(df_data=scenario2_summ,ed=1.04,title=expression(paste(bold("b")," ed = 1.04")))+theme(axis.text.y=element_blank())
p3 <- heatmap.accuracy.scenario2.func(df_data=scenario2_summ,ed=1.06,title=expression(paste(bold("c")," ed = 1.06")))+theme(axis.text.y=element_blank())
p1 %>% insert_left(margins,width=0.05) %>% 
  insert_right(p2, width=1) %>%
  insert_right(p3,width=1)
ggsave("figures/sim_scenario2_accuracy.jpg",width=18,height=6,dpi=300)


### heatmap for specificity
heatmap.specificity.scenario2.func <- function(df_data,ed,title){
  # heatmaps
  df_data <- df_data[df_data$ed==ed,]
  df_data <- df_data[df_data$norm_method%in%norm_methods$norm_method,]
  df_data$average_specificity <- round(df_data$average_specificity,3)
  heatmaps <- ggplot(df_data, aes(x=factor, y=annotated_method, fill=average_specificity))+ 
    geom_tile() +
    geom_text(aes(x=factor, y=annotated_method, label=average_specificity), color="black", size=4)+
    labs(x="factor",y="",fill="Average\nspecificity",title=title)+
    scale_fill_distiller(palette="Greens",limits=c(0,1),direction=1)+
    theme_bw()+
    theme(panel.grid.minor=element_blank(), panel.border=element_blank())
  heatmaps
}
p1 <- heatmap.specificity.scenario2.func(df_data=scenario2_summ,ed=1.02,title=expression(paste(bold("a")," ed = 1.02")))
p2 <- heatmap.specificity.scenario2.func(df_data=scenario2_summ,ed=1.04,title=expression(paste(bold("b")," ed = 1.04")))+theme(axis.text.y=element_blank())
p3 <- heatmap.specificity.scenario2.func(df_data=scenario2_summ,ed=1.06,title=expression(paste(bold("c")," ed = 1.06")))+theme(axis.text.y=element_blank())
p1 %>% insert_left(margins,width=0.05) %>% 
  insert_right(p2, width=1) %>%
  insert_right(p3,width=1)
ggsave("figures/sim_scenario2_specificity.jpg",width=18,height=6,dpi=300)


### heatmap for sensitivity
heatmap.sensitivity.scenario2.func <- function(df_data,ed,title){
  # heatmaps
  df_data <- df_data[df_data$ed==ed,]
  df_data <- df_data[df_data$norm_method%in%norm_methods$norm_method,]
  df_data$average_sensitivity <- round(df_data$average_sensitivity,3)
  heatmaps <- ggplot(df_data, aes(x=factor, y=annotated_method, fill=average_sensitivity))+ 
    geom_tile() +
    geom_text(aes(x=factor, y=annotated_method, label=average_sensitivity), color="black", size=4)+
    labs(x="factor",y="",fill="Average\nsensitivity",title=title)+
    scale_fill_distiller(palette="Blues",limits=c(0,1),direction=1)+
    theme_bw()+
    theme(panel.grid.minor=element_blank(), panel.border=element_blank())
  heatmaps
}
p1 <- heatmap.sensitivity.scenario2.func(df_data=scenario2_summ,ed=1.02,title=expression(paste(bold("a")," ed = 1.02")))
p2 <- heatmap.sensitivity.scenario2.func(df_data=scenario2_summ,ed=1.04,title=expression(paste(bold("b")," ed = 1.04")))+theme(axis.text.y=element_blank())
p3 <- heatmap.sensitivity.scenario2.func(df_data=scenario2_summ,ed=1.06,title=expression(paste(bold("c")," ed = 1.06")))+theme(axis.text.y=element_blank())
p1 %>% insert_left(margins,width=0.05) %>% 
  insert_right(p2, width=1) %>%
  insert_right(p3,width=1)
ggsave("figures/sim_scenario2_sensitivity.jpg",width=18,height=6,dpi=300)



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


### heatmap for auc
heatmap.auc.scenario3.func <- function(df_data,ed,title){
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
p1 <- heatmap.auc.scenario3.func(df_data=scenario3_summ,ed=1.02,title=expression(paste(bold("a")," ed = 1.02")))
p2 <- heatmap.auc.scenario3.func(df_data=scenario3_summ,ed=1.04,title=expression(paste(bold("b")," ed = 1.04")))+theme(axis.text.y=element_blank())
p3 <- heatmap.auc.scenario3.func(df_data=scenario3_summ,ed=1.06,title=expression(paste(bold("c")," ed = 1.06")))+theme(axis.text.y=element_blank())
p1 %>% insert_left(margins,width=0.05) %>% 
  insert_right(p2, width=1) %>%
  insert_right(p3,width=1)
ggsave("figures/sim_scenario3_auc.jpg",width=18,height=6,dpi=300)


### heatmap for accuracy
heatmap.accuracy.scenario3.func <- function(df_data,ed,title){
  # heatmaps
  df_data <- df_data[df_data$ed==ed,]
  df_data <- df_data[df_data$norm_method%in%norm_methods$norm_method,]
  df_data$average_accuracy <- round(df_data$average_accuracy,3)
  heatmaps <- ggplot(df_data, aes(x=overlap, y=annotated_method, fill=average_accuracy))+ 
    geom_tile() +
    geom_text(aes(x=overlap, y=annotated_method, label=average_accuracy), color="black", size=4)+
    labs(x="overlap",y="",fill="Average\naccuracy",title=title)+
    scale_fill_distiller(palette="Purples",limits=c(0.54,1),direction=1)+
    theme_bw()+
    theme(panel.grid.minor=element_blank(), panel.border=element_blank())
  heatmaps
}
p1 <- heatmap.accuracy.scenario3.func(df_data=scenario3_summ,ed=1.02,title=expression(paste(bold("a")," ed = 1.02")))
p2 <- heatmap.accuracy.scenario3.func(df_data=scenario3_summ,ed=1.04,title=expression(paste(bold("b")," ed = 1.04")))+theme(axis.text.y=element_blank())
p3 <- heatmap.accuracy.scenario3.func(df_data=scenario3_summ,ed=1.06,title=expression(paste(bold("c")," ed = 1.06")))+theme(axis.text.y=element_blank())
p1 %>% insert_left(margins,width=0.05) %>% 
  insert_right(p2, width=1) %>%
  insert_right(p3,width=1)
ggsave("figures/sim_scenario3_accuracy.jpg",width=18,height=6,dpi=300)


### heatmap for specificity
heatmap.specificity.scenario3.func <- function(df_data,ed,title){
  # heatmaps
  df_data <- df_data[df_data$ed==ed,]
  df_data <- df_data[df_data$norm_method%in%norm_methods$norm_method,]
  df_data$average_accuracy <- round(df_data$average_specificity,3)
  heatmaps <- ggplot(df_data, aes(x=overlap, y=annotated_method, fill=average_specificity))+ 
    geom_tile() +
    geom_text(aes(x=overlap, y=annotated_method, label=average_specificity), color="black", size=4)+
    labs(x="overlap",y="",fill="Average\nspecificity",title=title)+
    scale_fill_distiller(palette="Greens",limits=c(0.41,1),direction=1)+
    theme_bw()+
    theme(panel.grid.minor=element_blank(), panel.border=element_blank())
  heatmaps
}
p1 <- heatmap.specificity.scenario3.func(df_data=scenario3_summ,ed=1.02,title=expression(paste(bold("a")," ed = 1.02")))
p2 <- heatmap.specificity.scenario3.func(df_data=scenario3_summ,ed=1.04,title=expression(paste(bold("b")," ed = 1.04")))+theme(axis.text.y=element_blank())
p3 <- heatmap.specificity.scenario3.func(df_data=scenario3_summ,ed=1.06,title=expression(paste(bold("c")," ed = 1.06")))+theme(axis.text.y=element_blank())
p1 %>% insert_left(margins,width=0.05) %>% 
  insert_right(p2, width=1) %>%
  insert_right(p3,width=1)
ggsave("figures/sim_scenario3_specificity.jpg",width=18,height=6,dpi=300)


### heatmap for sensitivity
heatmap.sensitivity.scenario3.func <- function(df_data,ed,title){
  # heatmaps
  df_data <- df_data[df_data$ed==ed,]
  df_data <- df_data[df_data$norm_method%in%norm_methods$norm_method,]
  df_data$average_accuracy <- round(df_data$average_sensitivity,3)
  heatmaps <- ggplot(df_data, aes(x=overlap, y=annotated_method, fill=average_sensitivity))+ 
    geom_tile() +
    geom_text(aes(x=overlap, y=annotated_method, label=average_sensitivity), color="black", size=4)+
    labs(x="overlap",y="",fill="Average\nsensitivity",title=title)+
    scale_fill_distiller(palette="Blues",limits=c(0.41,1),direction=1)+
    theme_bw()+
    theme(panel.grid.minor=element_blank(), panel.border=element_blank())
  heatmaps
}
p1 <- heatmap.sensitivity.scenario3.func(df_data=scenario3_summ,ed=1.02,title=expression(paste(bold("a")," ed = 1.02")))
p2 <- heatmap.sensitivity.scenario3.func(df_data=scenario3_summ,ed=1.04,title=expression(paste(bold("b")," ed = 1.04")))+theme(axis.text.y=element_blank())
p3 <- heatmap.sensitivity.scenario3.func(df_data=scenario3_summ,ed=1.06,title=expression(paste(bold("c")," ed = 1.06")))+theme(axis.text.y=element_blank())
p1 %>% insert_left(margins,width=0.05) %>% 
  insert_right(p2, width=1) %>%
  insert_right(p3,width=1)
ggsave("figures/sim_scenario3_sensitivity.jpg",width=18,height=6,dpi=300)



                                      
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
# add the trn and tst information
crc_pred$group <- paste0("trn:",crc_pred$study1,", tst:",crc_pred$study2)
crc_pred_summ$group <- paste0("trn:",crc_pred_summ$study1,", tst:",crc_pred_summ$study2)


#=====================================================================================================================#
### box plot for auc values
crc_pred_auc <- list()
crc_groups <- unique(crc_pred$group)[order(unique(crc_pred$group))]
for(i in 1:length(crc_groups)){
  df <- crc_pred[crc_pred$group==crc_groups[i],]
  crc_pred_auc[[i]] <- ggplot(df,aes(x=annotated_method,y=auc_value,color=class))+
    geom_boxplot()+
    scale_color_jco()+
    labs(x="Normalization methods",y="AUC",title=crc_groups[i],color="Group of normalization methods")+
    theme_bw()+
    theme(plot.title=element_text(size=9))+
    theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))+
    theme(axis.title.x=element_text(size=7),axis.title.y=element_text(size=7))+
    theme(panel.grid.minor=element_blank())+
    theme(panel.border=element_rect(fill=NA,color="black"))+
    coord_flip()
}
ggarrange(plotlist=crc_pred_auc,nrow=8,ncol=7,labels=paste0(rep(letters[1:8],each=7),rep(1:7,8)),common.legend=T)
ggsave("figures/crc_pred_auc.jpg",width=15,height=22,dpi=300)


### box plot accuracy
crc_pred_accuracy <- list()
crc_groups <- unique(crc_pred$group)[order(unique(crc_pred$group))]
for(i in 1:length(crc_groups)){
  df <- crc_pred[crc_pred$group==crc_groups[i],]
  crc_pred_accuracy[[i]] <- ggplot(df,aes(x=annotated_method,y=accuracy,color=class))+
    geom_boxplot()+
    scale_color_jco()+
    labs(x="Normalization methods",y="accuracy",title=crc_groups[i],color="Group of normalization methods")+
    theme_bw()+
    theme(plot.title=element_text(size=9))+
    theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))+
    theme(axis.title.x=element_text(size=7),axis.title.y=element_text(size=7))+
    theme(panel.grid.minor=element_blank())+
    theme(panel.border=element_rect(fill=NA,color="black"))+
    coord_flip()
}
ggarrange(plotlist=crc_pred_accuracy,nrow=8,ncol=7,labels=paste0(rep(letters[1:8],each=7),rep(1:7,8)),common.legend=T)
ggsave("figures/crc_pred_accuracy.jpg",width=15,height=22,dpi=300)


### box plot specificity
crc_pred_specificity <- list()
crc_groups <- unique(crc_pred$group)[order(unique(crc_pred$group))]
for(i in 1:length(crc_groups)){
  df <- crc_pred[crc_pred$group==crc_groups[i],]
  crc_pred_specificity[[i]] <- ggplot(df,aes(x=annotated_method,y=specificity,color=class))+
    geom_boxplot()+
    scale_color_jco()+
    labs(x="Normalization methods",y="specificity",title=crc_groups[i],color="Group of normalization methods")+
    theme_bw()+
    theme(plot.title=element_text(size=9))+
    theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))+
    theme(axis.title.x=element_text(size=7),axis.title.y=element_text(size=7))+
    theme(panel.grid.minor=element_blank())+
    theme(panel.border=element_rect(fill=NA,color="black"))+
    coord_flip()
}
ggarrange(plotlist=crc_pred_specificity,nrow=8,ncol=7,labels=paste0(rep(letters[1:8],each=7),rep(1:7,8)),common.legend=T)
ggsave("figures/crc_pred_specificity.jpg",width=15,height=22,dpi=300)


### box plot sensitivity
crc_pred_sensitivity <- list()
crc_groups <- unique(crc_pred$group)[order(unique(crc_pred$group))]
for(i in 1:length(crc_groups)){
  df <- crc_pred[crc_pred$group==crc_groups[i],]
  crc_pred_sensitivity[[i]] <- ggplot(df,aes(x=annotated_method,y=sensitivity,color=class))+
    geom_boxplot()+
    scale_color_jco()+
    labs(x="Normalization methods",y="sensitivity",title=crc_groups[i],color="Group of normalization methods")+
    theme_bw()+
    theme(plot.title=element_text(size=9))+
    theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))+
    theme(axis.title.x=element_text(size=7),axis.title.y=element_text(size=7))+
    theme(panel.grid.minor=element_blank())+
    theme(panel.border=element_rect(fill=NA,color="black"))+
    coord_flip()
}
ggarrange(plotlist=crc_pred_sensitivity,nrow=8,ncol=7,labels=paste0(rep(letters[1:8],each=7),rep(1:7,8)),common.legend=T)
ggsave("figures/crc_pred_sensitivity.jpg",width=15,height=22,dpi=300)
                                    
                                      
#=====================================================================================================================#
### box plot for the ranks of auc
crc_ranks_auc <- ggboxplot(crc_pred_summ,x="annotated_method",y="order_auc",color="class")+
  scale_color_jco()+
  labs(x="Normalization methods",y="Ranks of AUC",color="Group of normalization methods",
       title="CRC Datasets: AUC Ranking Across Normalization Methods")+
  theme(axis.text.x=element_text(vjust=0.5,hjust=1,angle=90))+
  theme(legend.position="bottom")


### box plot for the ranks of accuracy
crc_ranks_accuracy <- ggboxplot(crc_pred_summ,x="annotated_method",y="order_accuracy",color="class")+
  scale_color_jco()+
  labs(x="Normalization methods",y="Ranks of accuracy",color="Group of normalization methods",
       title="CRC Datasets: Accuracy Ranking Across Normalization Methods")+
  theme(axis.text.x=element_text(vjust=0.5,hjust=1,angle=90))+
  theme(legend.position="bottom")

                                      
### box plot for the ranks of specificity
crc_ranks_specificity <- ggboxplot(crc_pred_summ,x="annotated_method",y="order_specificity",color="class")+
  scale_color_jco()+
  labs(x="Normalization methods",y="Ranks of specificity",color="Group of normalization methods",
       title="CRC Datasets: Specificity Ranking Across Normalization Methods")+
  theme(axis.text.x=element_text(vjust=0.5,hjust=1,angle=90))+
  theme(legend.position="bottom")

                                      
### box plot for the ranks of sensitivity
crc_ranks_sensitivity <- ggboxplot(crc_pred_summ,x="annotated_method",y="order_sensitivity",color="class")+
  scale_color_jco()+
  labs(x="Normalization methods",y="Ranks of sensitivity",color="Group of normalization methods",
       title="CRC Datasets: Sensitivity Ranking Across Normalization Methods")+
  theme(axis.text.x=element_text(vjust=0.5,hjust=1,angle=90))+
  theme(legend.position="bottom")


ggarrange(crc_ranks_auc,crc_ranks_accuracy,crc_ranks_sensitivity,crc_ranks_specificity,
          nrow=4,ncol=1,labels=letters[1:4],common.legend=T)
ggsave("figures/crc_pred_ranks.jpg",width=10,height=12,dpi=300)  


                                    

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
# add the trn and tst information
ibd_pred$group <- paste0("trn:",ibd_pred$study1,", tst:",ibd_pred$study2)
ibd_pred_summ$group <- paste0("trn:",ibd_pred_summ$study1,", tst:",ibd_pred_summ$study2)

                                      
#=====================================================================================================================#
### box plot for auc values
ibd_pred_auc <- list()
ibd_groups <- unique(ibd_pred$group)[order(unique(ibd_pred$group))]
for(i in 1:length(ibd_groups)){
  df <- ibd_pred[ibd_pred$group==ibd_groups[i],]
  ibd_pred_auc[[i]] <- ggplot(df,aes(x=annotated_method,y=auc_value,color=class))+
    geom_boxplot()+
    scale_color_jco()+
    labs(x="Normalization methods",y="AUC",title=ibd_groups[i],color="Group of normalization methods")+
    theme_bw()+
    theme(plot.title=element_text(size=9))+
    theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))+
    theme(axis.title.x=element_text(size=7),axis.title.y=element_text(size=7))+
    theme(panel.grid.minor=element_blank())+
    theme(panel.border=element_rect(fill=NA,color="black"))+
    coord_flip()
}
ggarrange(plotlist=ibd_pred_auc,nrow=5,ncol=4,labels=paste0(rep(letters[1:5],each=4),rep(1:4,5)),common.legend=T)
ggsave("figures/ibd_pred_auc.jpg",width=10,height=14,dpi=300)


### box plot accuracy
ibd_pred_accuracy <- list()
ibd_groups <- unique(ibd_pred$group)[order(unique(ibd_pred$group))]
for(i in 1:length(ibd_groups)){
  df <- ibd_pred[ibd_pred$group==ibd_groups[i],]
  ibd_pred_accuracy[[i]] <- ggplot(df,aes(x=annotated_method,y=accuracy,color=class))+
    geom_boxplot()+
    scale_color_jco()+
    labs(x="Normalization methods",y="accuracy",title=ibd_groups[i],color="Group of normalization methods")+
    theme_bw()+
    theme(plot.title=element_text(size=9))+
    theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))+
    theme(axis.title.x=element_text(size=7),axis.title.y=element_text(size=7))+
    theme(panel.grid.minor=element_blank())+
    theme(panel.border=element_rect(fill=NA,color="black"))+
    coord_flip()
}
ggarrange(plotlist=ibd_pred_accuracy,nrow=5,ncol=4,labels=paste0(rep(letters[1:5],each=4),rep(1:4,5)),common.legend=T)
ggsave("figures/ibd_pred_accuracy.jpg",width=10,height=14,dpi=300)


### box plot specificity
ibd_pred_specificity <- list()
ibd_groups <- unique(ibd_pred$group)[order(unique(ibd_pred$group))]
for(i in 1:length(ibd_groups)){
  df <- ibd_pred[ibd_pred$group==ibd_groups[i],]
  ibd_pred_specificity[[i]] <- ggplot(df,aes(x=annotated_method,y=specificity,color=class))+
    geom_boxplot()+
    scale_color_jco()+
    labs(x="Normalization methods",y="specificity",title=ibd_groups[i],color="Group of normalization methods")+
    theme_bw()+
    theme(plot.title=element_text(size=9))+
    theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))+
    theme(axis.title.x=element_text(size=7),axis.title.y=element_text(size=7))+
    theme(panel.grid.minor=element_blank())+
    theme(panel.border=element_rect(fill=NA,color="black"))+
    coord_flip()
}
ggarrange(plotlist=ibd_pred_specificity,nrow=5,ncol=4,labels=paste0(rep(letters[1:5],each=4),rep(1:4,5)),common.legend=T)
ggsave("figures/ibd_pred_specificity.jpg",width=10,height=14,dpi=300)

                                      

### box plot sensitivity
ibd_pred_sensitivity <- list()
ibd_groups <- unique(ibd_pred$group)[order(unique(ibd_pred$group))]
for(i in 1:length(ibd_groups)){
  df <- ibd_pred[ibd_pred$group==ibd_groups[i],]
  ibd_pred_sensitivity[[i]] <- ggplot(df,aes(x=annotated_method,y=sensitivity,color=class))+
    geom_boxplot()+
    scale_color_jco()+
    labs(x="Normalization methods",y="sensitivity",title=ibd_groups[i],color="Group of normalization methods")+
    theme_bw()+
    theme(plot.title=element_text(size=9))+
    theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))+
    theme(axis.title.x=element_text(size=7),axis.title.y=element_text(size=7))+
    theme(panel.grid.minor=element_blank())+
    theme(panel.border=element_rect(fill=NA,color="black"))+
    coord_flip()
}
ggarrange(plotlist=ibd_pred_sensitivity,nrow=5,ncol=4,labels=paste0(rep(letters[1:5],each=4),rep(1:4,5)),common.legend=T)
ggsave("figures/ibd_pred_sensitivity.jpg",width=10,height=14,dpi=300)
                                      


#=====================================================================================================================#
### box plot for the ranks of auc
ibd_ranks_auc <- ggboxplot(ibd_pred_summ,x="annotated_method",y="order_auc",color="class")+
  scale_color_jco()+
  labs(x="Normalization methods",y="Ranks of AUC",color="Group of normalization methods",
       title="IBD Datasets: AUC Ranking Across Normalization Methods")+
  theme(axis.text.x=element_text(vjust=0.5,hjust=1,angle=90))+
  theme(legend.position="bottom")


### box plot for the ranks of accuracy
ibd_ranks_accuracy <- ggboxplot(ibd_pred_summ,x="annotated_method",y="order_accuracy",color="class")+
  scale_color_jco()+
  labs(x="Normalization methods",y="Ranks of accuracy",color="Group of normalization methods",
       title="IBD Datasets: Accuracy Ranking Across Normalization Methods")+
  theme(axis.text.x=element_text(vjust=0.5,hjust=1,angle=90))+
  theme(legend.position="bottom")


### box plot for the ranks of specificity
ibd_ranks_specificity <- ggboxplot(ibd_pred_summ,x="annotated_method",y="order_specificity",color="class")+
  scale_color_jco()+
  labs(x="Normalization methods",y="Ranks of specificity",color="Group of normalization methods",
       title="IBD Datasets: Specificity Ranking Across Normalization Methods")+
  theme(axis.text.x=element_text(vjust=0.5,hjust=1,angle=90))+
  theme(legend.position="bottom")


### box plot for the ranks of sensitivity
ibd_ranks_sensitivity <- ggboxplot(ibd_pred_summ,x="annotated_method",y="order_sensitivity",color="class")+
  scale_color_jco()+
  labs(x="Normalization methods",y="Ranks of sensitivity",color="Group of normalization methods",
       title="IBD Datasets: Sensitivity Ranking Across Normalization Methods")+
  theme(axis.text.x=element_text(vjust=0.5,hjust=1,angle=90))+
  theme(legend.position="bottom")


ggarrange(ibd_ranks_auc,ibd_ranks_accuracy,ibd_ranks_sensitivity,ibd_ranks_specificity,
          nrow=4,ncol=1,labels=letters[1:4],common.legend=T)
ggsave("figures/ibd_pred_ranks.jpg",width=10,height=12,dpi=300)
















