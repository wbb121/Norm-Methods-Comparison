########################################################################################################################
### simulation in scenario 3, binary phenotype prediction
### run the code in SYS-7089P-TR4T sever
### 2023/09/09
########################################################################################################################

### set working directory
setwd("/home/wangbb/normalization_comparison")

### packages
all(sapply(c("DirichletReg","caret","pROC"), require, character.only=TRUE))

### load the helper function
source("helper.R")


#======================================================================================================================#
### arguments ###
#======================================================================================================================#
command_args=commandArgs(trailingOnly=T)
#command_args=c(100,1000000,10,"data/FengQ_2015_ctrl_count.rds","TSS","rfr")

# parameters
num_genes_overlap=c(2,4,6,8,10)             # number of overlapped genes
disease_effects=c(1.02,1.04,1.06)           # disease effect
sample_size=as.numeric(command_args[1])     # sample size
library_size=as.numeric(command_args[2])    # library size
num_genes=as.numeric(command_args[3])       # number of phenotype related gene 
count=readRDS(command_args[4])              # template data
norm_method=command_args[5]                 # normalization method
pred_method=command_args[6]                 # prediction method

# make the names for genes in two populations
rownames(count) <- make.names(rownames(count))

# choose phenotype related gene
set.seed(1234); selected_genes <- sample(rownames(count),num_genes,replace=F)
set.seed(1234); extra_genes <- sample(rownames(count)[!rownames(count)%in%selected_genes],num_genes,replace=F)
selected_genes_df <- as.data.frame(matrix(NA,nrow=length(num_genes_overlap),ncol=num_genes,
                                          dimnames=list(paste0("overlap",num_genes_overlap),paste0("genes",1:num_genes))))
for(i in 1:length(num_genes_overlap)){
  if(num_genes_overlap[i]<num_genes) 
    selected_genes_df[i,] <- c(selected_genes[1:(num_genes_overlap[i]/2)],
                               extra_genes[1:((num_genes-num_genes_overlap[i])/2)],
                               selected_genes[(num_genes/2+1):(num_genes/2+num_genes_overlap[i]/2)],
                               extra_genes[((num_genes-2)/2+1):((num_genes-2)/2+(num_genes-num_genes_overlap[i])/2)])
  else selected_genes_df[i,] <- selected_genes
}


#======================================================================================================================#
### functions ###
#======================================================================================================================#
### function for simulating count table
sim.count.func <- function(count,sample_size,library_size,selected_genes1,selected_genes2,disease_effect,seed){
  
  #selected_genes1=selected_genes;selected_genes2=as.character(selected_genes_df[5,])
  
  # the probability vectors from real data  
  prob <- rowSums(count)/sum(rowSums(count))
  
  # generate the case and control vectors with disease effect for different populations
  prob1_ctrl <- prob
  prob1_case <- prob
  prob1_case[selected_genes1[1:(num_genes/2)]] <- prob1_case[selected_genes1[1:(num_genes/2)]]*disease_effect
  prob1_case[selected_genes1[(num_genes/2+1):num_genes]] <- prob1_case[selected_genes1[(num_genes/2+1):num_genes]]/disease_effect
  prob1_case <- prob1_case/sum(prob1_case)
  prob2_ctrl <- prob
  prob2_case <- prob
  prob2_case[selected_genes2[1:(num_genes/2)]] <- prob2_case[selected_genes2[1:(num_genes/2)]]*disease_effect
  prob2_case[selected_genes2[(num_genes/2+1):num_genes]] <- prob2_case[selected_genes2[(num_genes/2+1):num_genes]]/disease_effect
  prob2_case <- prob2_case/sum(prob2_case)
  
  # generate the dirichlet distributed probability vectors for each sample based on the base probability vectors
  set.seed(seed); prob1s_ctrl <- rdirichlet(sample_size/2,1e6*prob1_ctrl)
  set.seed(seed); prob1s_case <- rdirichlet(sample_size/2,1e6*prob1_case)
  set.seed(seed+100); prob2s_ctrl <- rdirichlet(sample_size/2,1e6*prob2_ctrl)
  set.seed(seed+100); prob2s_case <- rdirichlet(sample_size/2,1e6*prob2_case)
  
  # simulate the counts table using multinomial distribution MN(library_size,prob)
  sim1_ctrl <- as.data.frame(matrix(NA,nrow=length(prob),ncol=sample_size/2,dimnames=list(names(prob),paste0("p1_ctrl_",1:(sample_size/2)))))
  sim1_case <- as.data.frame(matrix(NA,nrow=length(prob),ncol=sample_size/2,dimnames=list(names(prob),paste0("p1_case_",1:(sample_size/2)))))
  sim2_ctrl <- as.data.frame(matrix(NA,nrow=length(prob),ncol=sample_size/2,dimnames=list(names(prob),paste0("p2_ctrl_",1:(sample_size/2)))))
  sim2_case <- as.data.frame(matrix(NA,nrow=length(prob),ncol=sample_size/2,dimnames=list(names(prob),paste0("p2_case_",1:(sample_size/2)))))
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
if(!dir.exists("scenario3")) dir.create("scenario3")
if(!dir.exists("scenario3/sim_data")) dir.create("scenario3/sim_data")
# simulate the count table
for(disease_effect in disease_effects){
  for(j in 1:length(num_genes_overlap)){
    sim_tabs <- list()
    for(i in 1:100){
      sim_tabs[[i]] <- sim.count.func(count=count,sample_size=sample_size,library_size=library_size,
                                      selected_genes1=selected_genes,
                                      selected_genes2=as.character(selected_genes_df[j,]),
                                      disease_effect=disease_effect,seed=i)
    }
    saveRDS(sim_tabs,paste0("scenario3/sim_data/sim_ep",disease_effect,"_overlap",num_genes_overlap[j],".rds"))
    print(paste0("disease_effect=",disease_effect,",overlap_genes=",num_genes_overlap[j]))
  }
}
# After generating simulation data, comment out the simulation part so that data will not be generated repeatedly.


#======================================================================================================================#
### normalization ###
#======================================================================================================================#
if(!dir.exists("scenario3/norm_data")) dir.create("scenario3/norm_data")
if(!dir.exists(paste0("scenario3/norm_data/",norm_method))) dir.create(paste0("scenario3/norm_data/",norm_method))
for(disease_effect in disease_effects){
  for(j in 1:length(num_genes_overlap)){
    sim_tabs <- readRDS(paste0("scenario3/sim_data/sim_ep",disease_effect,"_overlap",num_genes_overlap[j],".rds"))
    norm_tabs <- list()
    for(i in 1:100){
      norm_tabs[[i]] <-  norm.func(p1=sim_tabs[[i]][[1]],p2=sim_tabs[[i]][[2]],norm_method=norm_method)
    }
    saveRDS(norm_tabs,paste0("scenario3/norm_data/",norm_method,"/norm_ep",disease_effect,"_overlap",num_genes_overlap[j],"_",norm_method,".rds"))
    print(paste0("disease_effect=",disease_effect,",overlap_genes=",num_genes_overlap[j],",norm_method=",norm_method))
  }
}


#======================================================================================================================#
### prediction ###
#======================================================================================================================#
if(!dir.exists("scenario3/pred_results")) dir.create("scenario3/pred_results")
if(!dir.exists(paste0("scenario3/pred_results/",norm_method))) dir.create(paste0("scenario3/pred_results/",norm_method))
for(disease_effect in disease_effects){
  for(j in 1:length(num_genes_overlap)){
    norm_tabs <- readRDS(paste0("scenario3/norm_data/",norm_method,"/norm_ep",disease_effect,"_overlap",num_genes_overlap[j],"_",norm_method,".rds"))
    pred_df <- data.frame(disease_effect=numeric(0),overlap_genes=numeric(0),norm_method=character(0),
                          pred_method=character(0),auc_value=numeric(0))
    for(i in 1:100){
      pred_res <- sim.pred.func(trn=norm_tabs[[i]][[1]],tst=norm_tabs[[i]][[2]],pred_method=pred_method)
      pred_df[i,] <- c(disease_effect=disease_effect,overlap_genes=num_genes_overlap[j],
                       norm_method=norm_method,pred_method=pred_method,auc_value=pred_res)
    }
    saveRDS(pred_df, paste0("scenario3/pred_results/",norm_method,"/pred_ed",disease_effect,"_overlap",num_genes_overlap[j],"_",norm_method,"_",pred_method,".rds"))
    print(paste0("disease_effect=",disease_effect,",overlap_genes=",num_genes_overlap[j],",norm_method=",norm_method,",pred_method=",pred_method))
  }
}




