########################################################################################################################
### validate the performance of different normalization methods on binary phenotype prediction in simulated datasets
### simulations
#   template data: controls in GuptaA_2019 and FengQ_2015
#   parameters: population effect ep=0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.4, 0.6, 0.8, 1
#               disease effect ed=1.02, 1.04,1.06
#               library size ls=1,000,000
#               sample size n=100
#               disease related gene nd=10, with half enriched, half depleted
### 2023/06/05
########################################################################################################################

setwd("/home/wangbb/normalization_comparison")

### packages
all(sapply(c("foreach","doParallel"), require, character.only=TRUE))

### load the helper function
source("helper.R")


#======================================================================================================================#
### function for simulate count table
sim.func <- function(population1,population2,alpha,ed,ns,ls,enriched_gene,depleted_gene){
  
  # alpha-population effect; ed-disease effect; ns-sample size; ls-library size; 
  # enriched_gene; depleted_gene
  
  # merge two populations
  merged_tab <- merge.func(population1,population2)
  tab1 <- merged_tab[,colnames(population1)]
  tab2 <- merged_tab[,colnames(population2)]
  genes <- rownames(merged_tab)
  
  # let v1 and v2 be the probability vectors for the two populations.
  v1 <- rowSums(tab1)/sum(rowSums(tab1))
  v2 <- rowSums(tab2)/sum(rowSums(tab2))
  # create a pseudo-population with relative abundance vector v1(alpha) = alpha*v1 +(1-alpha)*v2.
  # the differences between the two simulated populations (alpha*(v1-v2)) increases with alpha.
  v1 <- alpha*v1+(1-alpha)*v2
  
  # generate the case and control vectors
  # take half of the samples as controls without changing their probability
  v1_control <- v1
  v2_control <- v2
  # probability for half of the genes are multiplied by a disease effect ed and the other half by 1/ed
  v1_case <- v1
  v1_case[enriched_gene] <- v1_case[enriched_gene]*ed
  v1_case[depleted_gene] <- v1_case[depleted_gene]/ed
  v2_case <- v2
  v2_case[enriched_gene] <- v2_case[enriched_gene]*ed
  v2_case[depleted_gene] <- v2_case[depleted_gene]/ed
  # normalize the probability vector
  v1_control <- v1_control/sum(v1_control)
  v2_control <- v2_control/sum(v2_control)
  v1_case <- v1_case/sum(v1_case)
  v2_case <- v2_case/sum(v2_case)
  
  # get the counts table using multinomial distribution MN(library size, v)
  sim1_control <- as.data.frame(matrix(NA,nrow=length(genes),ncol=ns/2,dimnames=list(genes,paste0("p1_control_",1:(ns/2)))))
  sim2_control <- as.data.frame(matrix(NA,nrow=length(genes),ncol=ns/2,dimnames=list(genes,paste0("p2_control_",1:(ns/2)))))
  sim1_case <- as.data.frame(matrix(NA,nrow=length(genes),ncol=ns/2,dimnames=list(genes,paste0("p1_case_",1:(ns/2)))))
  sim2_case <- as.data.frame(matrix(NA,nrow=length(genes),ncol=ns/2,dimnames=list(genes,paste0("p2_case_",1:(ns/2)))))
  for(i in 1:(ns/2)){
    sim1_control[,i] <- table(sample(genes, size=ls, prob=v1_control, replace=TRUE))[genes]
    sim2_control[,i] <- table(sample(genes, size=ls, prob=v2_control, replace=TRUE))[genes]
    sim1_case[,i] <- table(sample(genes, size=ls, prob=v1_case, replace=TRUE))[genes]
    sim2_case[,i] <- table(sample(genes, size=ls, prob=v2_case, replace=TRUE))[genes]
  }  
  
  # follow-up processing of simulated count tables
  # combine case and control samples for sim1 and sim2
  sim1 <- cbind(sim1_control,sim1_case)
  sim2 <- cbind(sim2_control,sim2_case)
  # replace NA values using 0
  sim1[is.na(sim1)] <- 0
  sim2[is.na(sim2)] <- 0
  # remove all zero genes
  sim1_final <- sim1[rowSums(sim1)>0,]
  sim2_final <- sim2[rowSums(sim2)>0,]
  
  # return
  return(list(sim1_final,sim2_final))
}


#======================================================================================================================#
### arguments
command_args=commandArgs(trailingOnly=T)
#command_args=c(50, 1.02, 100, 1000000, 10, "data/GuptaA_2019_ctrl_count.rds", "data/FengQ_2015_ctrl_count.rds",
#               50,"TSS",10,"ranger")

### parameters
sim_cluster=as.numeric(command_args[1])  # number of clusters to do the simulation
ep=c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.4, 0.6, 0.8, 1)           # population effect
ed=as.numeric(command_args[2])           # disease effect
ns=as.numeric(command_args[3])           # sample size
ls=as.numeric(command_args[4])           # library size
nd=as.numeric(command_args[5])           # number of disease related genes
population1=readRDS(command_args[6])     # simulation template 1
population2=readRDS(command_args[7])     # simulation template 2
norm_cluster=as.numeric(command_args[8]) # number of clusters to do the normalization
norm_method=command_args[9]              # normalization method
pred_cluster=as.numeric(command_args[10])# number of clusters to do the prediction
pred_method=command_args[11]             # prediction method


#======================================================================================================================#
### simulation ###
# template data processing, make the names for genes in two populations
rownames(population1) <- make.names(rownames(population1))
rownames(population2) <- make.names(rownames(population2))

# choose disease related gene
set.seed(1234)
inter_gene <- intersect(rownames(population1),rownames(population2))
disease_gene <- sample(inter_gene,nd,replace=F)
enriched_gene <- disease_gene[1:(nd/2)]
depleted_gene <- disease_gene[(nd/2+1):nd]

# simulation
if(!dir.exists("./sim_data")) dir.create("./sim_data")
for(i in ep){
  cl<- makeCluster(sim_cluster,setup_strategy="sequential")      
  registerDoParallel(cl)
  sim_tabs <- foreach(m=1:100) %dopar% {
    sim.func(population1=population1,population2=population2,alpha=i,ed=ed,ns=ns,ls=ls,
             enriched_gene=enriched_gene,depleted_gene=depleted_gene)
  }
  stopCluster(cl)
  saved_file <- paste0("./sim_data/sim_nd",nd,"_ed",ed,"_alpha",i,".rds")
  saveRDS(sim_tabs, saved_file)
  print(paste("nd=",nd,", ed=",ed,", alpha=",i))
}


#======================================================================================================================#
### normalization ###
for(i in ep){
  # simulated data
  sim_data_file <- paste0("./sim_data/sim_nd",nd,"_ed",ed,"_alpha",i,".rds")
  sim_data <- readRDS(sim_data_file)
    
  # normalization
  cl<- makeCluster(norm_cluster,setup_strategy="sequential")      
  registerDoParallel(cl)
  norm_tabs <- foreach(x=1:100, .packages=c("metagenomeSeq","edgeR","DESeq2","GUniFrac","compositions","sva")) %dopar% {
    norm.func(p1=sim_data[[x]][[1]],p2=sim_data[[x]][[2]],norm_method=norm_method)
  }
  stopCluster(cl)
    
  # save the results
  if(!dir.exists("sim_data_norm")) dir.create("sim_data_norm")
  if(!dir.exists(paste0("./sim_data_norm/",norm_method))) dir.create(paste0("./sim_data_norm/",norm_method))
  saved_file <- paste0("./sim_data_norm/",norm_method,"/norm_nd",nd,"_ed",ed,"_alpha",i,"_",norm_method,".rds")
  saveRDS(norm_tabs, saved_file)
  print(paste("nd=",nd,", ed=",ed,", alpha=",i,", norm_method=",norm_method))
}


#======================================================================================================================#
### prediction ###
# data frame for saving the results
pred_df <- data.frame(alpha=NA,ed=NA,nd=NA,norm_method=NA,pred_method=NA,auc_value=NA)
pred_df <- pred_df[-1,]

# prediction
for(i in ep){
  
  # normalized count tabs
  norm_data_file <- paste0("./sim_data_norm/",norm_method,"/norm_nd",nd,"_ed",ed,"_alpha",i,"_",norm_method,".rds")
  norm_data <- readRDS(norm_data_file)
  
  # prediction
  cl<- makeCluster(pred_cluster,setup_strategy="sequential")      
  registerDoParallel(cl)
  auc_values <- foreach(x=1:100,.packages=c("caret","pROC"),.combine=c) %dopar% {
    sim.pred.func(trn=norm_data[[x]][[1]],tst=norm_data[[x]][[2]],pred_method=pred_method)
  }
  stopCluster(cl)
  
  # summarize the results
  df <- data.frame(alpha=i,ed=ed,nd=nd,norm_method=norm_method,pred_method=pred_method,auc_value=auc_values)
  pred_df <- rbind(pred_df,df)
  print(paste("nd=",nd,", ed=",ed,", alpha=",i,", norm_method=",norm_method,", pred_method=",pred_method))
}

# save the results
if(!dir.exists("./sim_data_pred")) dir.create("./sim_data_pred")
saved_file <- paste0("./sim_data_pred/auc_nd",nd,"_ed",ed,"_",norm_method,"_",pred_method,".rds")
saveRDS(pred_df, saved_file)








