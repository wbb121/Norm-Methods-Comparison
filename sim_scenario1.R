########################################################################################################################
### simulation scenario 1: different background distributions of taxa in populations
### parameters: population effect ep=0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.4, 0.6, 0.8, 1
#               disease effect ed=1.02, 1.04,1.06
#               library size ls=1000000
#               sample size n=100
#               disease related gene nd=10, with half enriched, half depleted
### 2023/09/28
########################################################################################################################

setwd("/home/wangbb/normalization_comparison")

### packages
all(sapply(c("DirichletReg","caret","pROC"), require, character.only=TRUE))

### load the helper function
source("helper.R")

#======================================================================================================================#
### arguments ###
#======================================================================================================================#
command_args=commandArgs(trailingOnly=T)
#command_args=c(100,1000000,10,"data/GuptaA_2019_ctrl_count.rds","data/FengQ_2015_ctrl_count.rds","TSS","rfr")

### parameters
population_effects=c(0,0.05,0.1,0.15,0.2,0.25,0.4,0.6,0.8,1)  # population effect
disease_effects=c(1.02,1.04,1.06)                             # disease effect
sample_size=as.numeric(command_args[1])     # sample size
library_size=as.numeric(command_args[2])    # library size
num_genes=as.numeric(command_args[3])       # number of phenotype related gene 
count1=readRDS(command_args[4])             # template data1
count2=readRDS(command_args[5])             # template data2
norm_method=command_args[6]                 # normalization method
pred_method=command_args[7]                 # prediction method

# make the names for genes in two populations
rownames(count1) <- make.names(rownames(count1))
rownames(count2) <- make.names(rownames(count2))
# choose phenotype related gene
inter_genes <- intersect(rownames(count1),rownames(count2))
set.seed(1234); selected_genes <- sample(inter_genes,num_genes,replace=F)
enriched_genes <- selected_genes[1:(num_genes/2)]
depleted_genes <- selected_genes[(num_genes/2+1):num_genes]


#======================================================================================================================#
### functions ###
#======================================================================================================================#
### function for simulating the count table
#alpha=0.2;disease_effect=1.02;seed=1234
sim.count.func <- function(count1,count2,sample_size,library_size,alpha,disease_effect,enriched_genes,depleted_genes,seed){
  
  # the probability vectors from real data
  tab1 <- merge.func(count1,count2)[,colnames(count1)]
  tab2 <- merge.func(count1,count2)[,colnames(count2)]
  prob1 <- rowSums(tab1)/sum(rowSums(tab1))
  prob2 <- rowSums(tab2)/sum(rowSums(tab2))
  
  # pseudo probability vectors with population effect alpha
  prob1 <- alpha*prob1+(1-alpha)*prob2
  prob1 <- prob1/sum(prob1)
  prob1 <- prob1[prob1!=0]
  prob2 <- prob2[prob2!=0]
  
  # generate the case and control vectors with disease effect for different populations
  prob1_ctrl <- prob1
  prob1_case <- prob1
  prob1_case[enriched_genes] <- prob1_case[enriched_genes]*disease_effect
  prob1_case[depleted_genes] <- prob1_case[depleted_genes]/disease_effect
  prob1_case <- prob1_case/sum(prob1_case)
  prob2_ctrl <- prob2
  prob2_case <- prob2
  prob2_case[enriched_genes] <- prob2_case[enriched_genes]*disease_effect
  prob2_case[depleted_genes] <- prob2_case[depleted_genes]/disease_effect
  prob2_case <- prob2_case/sum(prob2_case)
  
  # generate the dirichlet distributed probability vectors for each sample based on the base probability vectors
  set.seed(seed); prob1s_ctrl <- rdirichlet(sample_size/2,1e6*prob1_ctrl)
  set.seed(seed); prob1s_case <- rdirichlet(sample_size/2,1e6*prob1_case)
  set.seed(seed); prob2s_ctrl <- rdirichlet(sample_size/2,1e6*prob2_ctrl)
  set.seed(seed); prob2s_case <- rdirichlet(sample_size/2,1e6*prob2_case)
  
  # simulate the counts table using multinomial distribution MN(library_size,prob)
  sim1_ctrl <- as.data.frame(matrix(NA,nrow=length(prob1),ncol=sample_size/2,dimnames=list(names(prob1),paste0("p1_ctrl_",1:(sample_size/2)))))
  sim1_case <- as.data.frame(matrix(NA,nrow=length(prob1),ncol=sample_size/2,dimnames=list(names(prob1),paste0("p1_case_",1:(sample_size/2)))))
  sim2_ctrl <- as.data.frame(matrix(NA,nrow=length(prob2),ncol=sample_size/2,dimnames=list(names(prob2),paste0("p2_ctrl_",1:(sample_size/2)))))
  sim2_case <- as.data.frame(matrix(NA,nrow=length(prob2),ncol=sample_size/2,dimnames=list(names(prob2),paste0("p2_case_",1:(sample_size/2)))))
  for(i in 1:(sample_size/2)){
    set.seed(i); sim1_ctrl[,i] <- rmultinom(1, size=library_size, prob=prob1s_ctrl[i,])
    set.seed(i); sim1_case[,i] <- rmultinom(1, size=library_size, prob=prob1s_case[i,])
    set.seed(i); sim2_ctrl[,i] <- rmultinom(1, size=library_size, prob=prob2s_ctrl[i,])
    set.seed(i); sim2_case[,i] <- rmultinom(1, size=library_size, prob=prob2s_case[i,])
  }
  sim1 <- cbind(sim1_ctrl,sim1_case); sim1 <- sim1[rowSums(sim1)>0,]
  sim2 <- cbind(sim2_ctrl,sim2_case); sim2 <- sim2[rowSums(sim2)>0,]
  
  # return
  return(list(sim1,sim2))
}


#======================================================================================================================#
### simulation ###
#======================================================================================================================#
if(!dir.exists("scenario1")) dir.create("scenario1")
if(!dir.exists("scenario1/sim_data")) dir.create("scenario1/sim_data")
for(alpha in population_effects){
  for(disease_effect in disease_effects){
    sim_tabs <- list()
    for(i in 1:100){
      sim_tabs[[i]] <- sim.count.func(count1=count1,count2=count2,sample_size=sample_size,library_size=library_size,alpha=alpha,
                                      disease_effect=disease_effect,enriched_genes=enriched_genes,depleted_genes=depleted_genes,seed=i)
    }
    saveRDS(sim_tabs,paste0("scenario1/sim_data/sim_alpha",alpha,"_ed",disease_effect,".rds"))
    print(paste0("alpha=",alpha,",disease_effect=",disease_effect))
  }
}
# After generating simulation data, comment out the simulation part so that data will not be generated repeatedly.


#======================================================================================================================#
### normalization ###
#======================================================================================================================#
if(!dir.exists("scenario1/norm_data")) dir.create("scenario1/norm_data")
if(!dir.exists(paste0("scenario1/norm_data/",norm_method))) dir.create(paste0("scenario1/norm_data/",norm_method))
for(alpha in population_effects){
  for(disease_effect in disease_effects){
    sim_tabs <- readRDS(paste0("scenario1/sim_data/sim_alpha",alpha,"_ed",disease_effect,".rds"))
    norm_tabs <- list()
    for(i in 1:100){
      norm_tabs[[i]] <- norm.func(p1=sim_tabs[[i]][[1]],p2=sim_tabs[[i]][[2]],norm_method=norm_method)
    }
    saveRDS(norm_tabs,paste0("scenario1/norm_data/",norm_method,"/norm_alpha",alpha,"_ed",disease_effect,"_",norm_method,".rds"))
    print(paste0("alpha=",alpha,",disease_effect=",disease_effect,",norm_method=",norm_method))
  }
}


#======================================================================================================================#
### prediction ###
#======================================================================================================================#
if(!dir.exists("scenario1/pred_results")) dir.create("scenario1/pred_results")
if(!dir.exists(paste0("scenario1/pred_results/",norm_method))) dir.create(paste0("scenario1/pred_results/",norm_method))
for(alpha in population_effects){
  for(disease_effect in disease_effects){
    norm_tabs <- readRDS(paste0("scenario1/norm_data/",norm_method,"/norm_alpha",alpha,"_ed",disease_effect,"_",norm_method,".rds"))
    pred_df <- data.frame(alpha=numeric(0),ed=numeric(0),num_genes=numeric(0),
                          norm_method=character(0),pred_method=character(0),auc_value=numeric(0))
    for(i in 1:100){
      pred_res <- sim.pred.func(trn=norm_tabs[[i]][[1]],tst=norm_tabs[[i]][[2]],pred_method=pred_method)
      pred_df[i,] <- c(alpha=alpha,ed=disease_effect,num_genes=num_genes,
                       norm_method=norm_method,pred_method=pred_method,auc_value=pred_res)
    }
    saveRDS(pred_df, paste0("scenario1/pred_results/",norm_method,"/norm_alpha",alpha,"_ed",disease_effect,"_",norm_method,"_",pred_method,".rds"))
    print(paste0("alpha=",alpha,",disease_effect=",disease_effect,",norm_method=",norm_method,",pred_method=",pred_method))
  }
}























