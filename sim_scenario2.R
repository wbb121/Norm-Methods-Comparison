########################################################################################################################
### simulation scenario 2: different batch effects in studies with the same background distribution of taxa in populations
### parameters: batch mean batch_means=c(0,500,1000)
#               batch variance batch_vars=c(1,2,4)
#               disease effect ed=1.02, 1.04,1.06
#               library size ls=1000000
#               sample size n=100
#               disease related gene nd=10, with half enriched, half depleted
### 2023/09/28
########################################################################################################################

### set working directory
setwd("/home/wangbb/normalization/simulations/")

### packages
all(sapply(c("DirichletReg","MCMCpack","caret","pROC"), require, character.only=TRUE))

### load the helper function
source("helper.R")

#======================================================================================================================#
### arguments ###
#======================================================================================================================#
command_args=commandArgs(trailingOnly=T)
#command_args=c(100,1000000,10,"data/FengQ_2015_ctrl_count.rds","TSS","rfr")

# parameters
batch_means=c(0,500,1000)                   # batch effect on the mean
batch_vars=c(1,2,4)                         # batch effect on the variance
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
enriched_genes <- selected_genes[1:(num_genes/2)]
depleted_genes <- selected_genes[(num_genes/2+1):num_genes]


#======================================================================================================================#
### functions ###
#======================================================================================================================#
### functions related to inverse gamma distribution
# transform the mean and variance to scale and shape parameters
mv2ab <- function(m, v){
  a <- 2 + m^2/v
  b <- m * (a-1)
  return(list(alpha=a, beta=b))
}
# transform the scale and shape parameters to mean and variance
ab2mv <- function(a, b){
  m <- b / (a-1)
  v <- b^2 / ((a-1)^2*(a-2))
  return(list(mean=m, var=v))
}


### function for simulating count table
sim.count.func <- function(count,sample_size,library_size,enriched_genes,depleted_genes,disease_effect,seed){
  
  # the probability vector from real data
  prob <- rowSums(count)/sum(rowSums(count))
  
  # generate the case and control vectors with disease effect for different populations
  prob_ctrl <- prob
  prob_case <- prob
  prob_case[enriched_genes] <- prob_case[enriched_genes]*disease_effect
  prob_case[depleted_genes] <- prob_case[depleted_genes]/disease_effect
  prob_case <- prob_case/sum(prob_case)
  
  # generate the probability vectors for each sample based on the base probability vectors
  set.seed(seed); prob1s_ctrl <- rdirichlet(sample_size/2,1e6*prob_ctrl)
  set.seed(seed); prob1s_case <- rdirichlet(sample_size/2,1e6*prob_case)
  set.seed(seed+100); prob2s_ctrl <- rdirichlet(sample_size/2,1e6*prob_ctrl)
  set.seed(seed+100); prob2s_case <- rdirichlet(sample_size/2,1e6*prob_case)
  
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


### function for simulating batch effect
sim.batch.func <- function(count_table,hyper_pars,seed){
  
  #count_table=sim.count.func(count=count,sample_size=sample_size,library_size=library_size,enriched_genes=enriched_genes,depleted_genes=depleted_genes,disease_effect=disease_effect,seed=seed)[[1]];
  #hyper_pars=list(hyper_mu=batch_mean, hyper_sd=sqrt(0.01),hyper_alpha=mv2ab(m=batch_var,v=0.01)$alpha,hyper_beta=mv2ab(m=batch_var,v=0.01)$beta)
  
  # Simulate batch parameters from hyper-pars
  set.seed(seed); gamma <- rnorm(nrow(count_table), mean=hyper_pars$hyper_mu, sd=hyper_pars$hyper_sd)
  set.seed(seed); delta2 <- rinvgamma(nrow(count_table), shape=hyper_pars$hyper_alpha, scale=hyper_pars$hyper_beta)
  
  # status
  status=as.factor(gsub("p[0-9]_([a-z]+)_[0-9]+","\\1",colnames(count_table)))
  
  # Simulate batch effect
  # fit linear model to data with no batch parameters, calculate residual variance
  X <- model.matrix(~Condition, data=data.frame(Condition=status))
  beta <- solve(t(X) %*% X) %*% t(X) %*% t(count_table)
  resid <- count_table - t(X %*% beta)
  # spike-in batch variance: multiply by condition adjusted data with delta
  resid_varbatch <- resid*sqrt(delta2)
  # construct mean batch parameter design matrix for adding gamma
  X_batch <- model.matrix(~-1+Batch, data=data.frame(Batch=rep(1,ncol(count_table))))
  # new data with added batch effect
  new_count_table <- t(cbind(X, X_batch) %*% rbind(beta, gamma)) + resid_varbatch
  # round the new data, making it integers
  new_count_table <- round(new_count_table,0)
  # if the count<0, set it 0
  new_count_table[new_count_table<0] <- 0
  
  # return
  return(new_count_table)
  
}


#======================================================================================================================#
### simulation ###
#======================================================================================================================#
if(!dir.exists("scenario2/sim_data")) dir.create("scenario2/sim_data")
# simulate the count table
for(disease_effect in disease_effects){
  sim_tabs <- list()
  for(i in 1:100){
    sim_tabs[[i]] <- sim.count.func(count=count,sample_size=sample_size,library_size=library_size,
                                    enriched_genes=enriched_genes,depleted_genes=depleted_genes,
                                    disease_effect=disease_effect,seed=i)
  }
  saveRDS(sim_tabs,paste0("scenario2/sim_data/sim_ed",disease_effect,".rds"))
  print(paste0("disease_effect=",disease_effect))
}


# simulate the batch effect
for(disease_effect in disease_effects){
  sim_tabs <- readRDS(paste0("scenario2/sim_data/sim_ed",disease_effect,".rds"))
  for(batch_mean in batch_means){
    for(batch_var in batch_vars){
      # hyper-parameters for simulating the batch effect
      # the mean of gene g caused by batch: gamma_g~N(mu,sd2)
      # the variance of gene g caused by batch: delta_g~InvGamma(alpha,beta)
      hyper_pars <- list(hyper_mu=batch_mean, hyper_sd=sqrt(0.01),
                         hyper_alpha=mv2ab(m=batch_var,v=0.01)$alpha,
                         hyper_beta=mv2ab(m=batch_var,v=0.01)$beta)
      # simulate the batch effect
      sim_tabs_adjust <- list()
      for(i in 1:100){
        sim_tabs_adjust_trn <- sim.batch.func(count_table=sim_tabs[[i]][[1]],hyper_pars=hyper_pars,seed=i)
        sim_tabs_adjust_tst <- sim_tabs[[i]][[2]]
        sim_tabs_adjust[[i]] <- list(sim_tabs_adjust_trn,sim_tabs_adjust_tst)
      }
      saveRDS(sim_tabs_adjust,paste0("scenario2/sim_data/sim_ed",disease_effect,"_mean",batch_mean,"_var",batch_var,".rds"))
      print(paste0("disease_effect=",disease_effect,",batch_mean=",batch_mean,",batch_var=",batch_var))
    }
  }
}

#======================================================================================================================#
### normalization ###
#======================================================================================================================#
if(!dir.exists("scenario2/norm_data")) dir.create("scenario2/norm_data")
if(!dir.exists(paste0("scenario2/norm_data/",norm_method))) dir.create(paste0("scenario2/norm_data/",norm_method))
for(disease_effect in disease_effects){
  for(batch_mean in batch_means){
    for(batch_var in batch_vars){
      sim_tabs <- readRDS(paste0("scenario2/sim_data/sim_ed",disease_effect,"_mean",batch_mean,"_var",batch_var,".rds"))
      norm_tabs <- list()
      for(i in 1:100){
        norm_tabs[[i]] <-  norm.func(p1=sim_tabs[[i]][[1]],p2=sim_tabs[[i]][[2]],norm_method=norm_method)
      }
      saveRDS(norm_tabs, paste0("scenario2/norm_data/",norm_method,"/norm_ed",disease_effect,"_mean",batch_mean,"_var",batch_var,"_",norm_method,".rds"))
      print(paste0("disease_effect=",disease_effect,",batch_mean=",batch_mean,",batch_var=",batch_var,",norm_method=",norm_method))
    }
  }
}


#======================================================================================================================#
### prediction ###
#======================================================================================================================#
if(!dir.exists("scenario2/pred_results")) dir.create("scenario2/pred_results")
if(!dir.exists(paste0("scenario2/pred_results/",norm_method))) dir.create(paste0("scenario2/pred_results/",norm_method))
for(disease_effect in  disease_effects){
  for(batch_mean in batch_means){
    for(batch_var in batch_vars){
      norm_tabs <- readRDS(paste0("scenario2/norm_data/",norm_method,"/norm_ed",disease_effect,"_mean",batch_mean,"_var",batch_var,"_",norm_method,".rds"))
      pred_df <- data.frame(ed=numeric(0),batch_mean=numeric(0),batch_var=numeric(0),num_genes=numeric(0),
                            norm_method=character(0),pred_method=character(0),auc_value=numeric(0))
      for(i in 1:100){
        pred_res <- sim.pred.func(trn=norm_tabs[[i]][[1]],tst=norm_tabs[[i]][[2]],pred_method=pred_method)
        pred_df[i,] <- c(ed=disease_effect,batch_mean=batch_mean,batch_var=batch_var,num_genes=num_genes,
                         norm_method=norm_method,pred_method=pred_method,auc_value=pred_res)
      }
      saveRDS(pred_df,paste0("scenario2/pred_results/",norm_method,"/pred_ed",disease_effect,"_mean",batch_mean,"_var",batch_var,"_",norm_method,"_",pred_method,".rds"))
      print(paste0("batch_mean=",batch_mean,",batch_var=",batch_var,",norm_method=",norm_method,",pred_method=",pred_method))
    }
  }
}











