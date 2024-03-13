########################################################################################################################
### helper for the comparsions of normalization methods for metagenomic cross-study phenotype prediction
### 2023/09/28
########################################################################################################################


#======================================================================================================================#
### function for merging two count table
merge.func <- function(table1,table2){
  table <- merge(table1,table2,by="row.names",all=T)
  rownames(table) <- table$Row.names
  table <- table[,-grep("Row.names",colnames(table))]
  table[is.na(table)] <- 0
  return(table)
}


#======================================================================================================================#
### function for obtaining the count table from curatedMetagenomicData
obtain.count.table <- function(study_name,metadata){
  study_list <- curatedMetagenomicData(paste0(study_name,".relative_abundance"), dryrun=F, counts=T)
  study_obj <- study_list[[1]]
  study_count_table <- as.data.frame(assays(study_obj)[[1]])
  # filter the samples according to the metadata
  study_count_table <- study_count_table[,colnames(study_count_table)%in%metadata$sample_id]
  # remove the all zero features
  study_count_table <- study_count_table[rowSums(study_count_table)>0,]
  # return
  return(study_count_table)
}


#======================================================================================================================#
### function for permanova analysis
permanova.func <- function(distance,meta_df,variable){
  meta_df <- meta_df[!is.na(meta_df[,variable]),]
  distance <- as.dist(as.matrix(distance)[meta_df$sample_id,meta_df$sample_id])
  formulas <- as.formula(paste0("distance~",variable))
  permanova <- adonis2(formulas,data=meta_df,permutations=1000,parallel=50)
  c(variable=variable,r2=round(permanova[variable,"R2"],3),pvalue=round(permanova[variable,"Pr(>F)"],3))
}


#======================================================================================================================#
### function for normalize the data
norm.func <- function(p1,p2,norm_method){
  
  # remove the all zero genes
  p1 <- p1[rowSums(p1)>0,]
  p2 <- p2[rowSums(p2)>0,]
  
  # remove samples with only 1 non-zero features
  p1 <- p1[,colSums(p1!=0)>1]
  p2 <- p2[,colSums(p2!=0)>1]
  
  # make names for genes
  rownames(p1) <- make.names(rownames(p1))
  rownames(p2) <- make.names(rownames(p2))
  
  # TSS, for samples
  if(norm_method=="TSS"){
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    norm_p2 <- as.data.frame(apply(p2,2,function(x) x/sum(x)))
    # let p2 have the same genes as p1 
    merged <- merge.func(norm_p1,norm_p2)
    final_p1 <- merged[rownames(norm_p1),colnames(norm_p1)]
    final_p2 <- merged[rownames(norm_p1),colnames(norm_p2)]
  }
  
  # UQ, for samples
  if(norm_method=="UQ"){
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/quantile(x[x>0])["75%"]))
    norm_p2 <- as.data.frame(apply(p2,2,function(x) x/quantile(x[x>0])["75%"]))
    # let p2 have the same genes as p1 
    merged <- merge.func(norm_p1,norm_p2)
    final_p1 <- merged[rownames(norm_p1),colnames(norm_p1)]
    final_p2 <- merged[rownames(norm_p1),colnames(norm_p2)]
  }
  
  # MED, for samples
  if(norm_method=="MED"){
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/median(x[x>0])))
    norm_p2 <- as.data.frame(apply(p2,2,function(x) x/median(x[x>0])))
    # let p2 have the same genes as p1 
    merged <- merge.func(norm_p1,norm_p2)
    final_p1 <- merged[rownames(norm_p1),colnames(norm_p1)]
    final_p2 <- merged[rownames(norm_p1),colnames(norm_p2)]
  }
  
  # CSS, for samples
  if(norm_method=="CSS"){
    require(metagenomeSeq)
    # function for CSS normalization
    css.func <- function(tab){
      tab_mrexperiment <- newMRexperiment(tab)
      tab_css <- cumNorm(tab_mrexperiment)
      tab_norm <- MRcounts(tab_css, norm=T)
      as.data.frame(tab_norm)
    }
    # css normalization
    norm_p1 <- css.func(p1)
    norm_p2 <- css.func(p2)
    # let p2 have the same genes as p1 
    merged <- merge.func(norm_p1,norm_p2)
    final_p1 <- merged[rownames(norm_p1),colnames(norm_p1)]
    final_p2 <- merged[rownames(norm_p1),colnames(norm_p2)]
  }
  
  # TMM, for samples, need to choose a reference
  # performed normalization on the training data, and then performed addon normalization of the test data onto the training data, 
  # to ensure that the normalization of the training data does not in any sense depend on the testing data.
  if(norm_method=="TMM"){
    require(edgeR)
    # function for TMM normalization
    tmm.func <- function(tab){
      tab <- as.matrix(tab)
      tab_dge <- DGEList(counts=tab)
      tab_tmm_dge <- calcNormFactors(tab_dge, method="TMM")
      tab_norm <- cpm(tab_tmm_dge)
      as.data.frame(tab_norm)
    }
    # TMM normalization
    norm_p1 <- tmm.func(p1)
    norm_p2 <- tmm.func(merge.func(p1,p2))[rownames(p1),colnames(p2)]
    # let p2 have the same genes as p1 
    final_p1 <- norm_p1
    final_p2 <- norm_p2
  }
  
  # TMM+, for samples, need to choose a reference
  # add a pseudo count 1 to zero values
  # performed normalization on the training data, and then performed addon normalization of the test data onto the training data, 
  # to ensure that the normalization of the training data does not in any sense depend on the testing data.
  if(norm_method=="TMM+"){
    require(edgeR)
    # replaced all the 0 abundances with 1 (need integers as input)
    p1[p1==0] <- 1
    p2[p2==0] <- 1
    # function for TMM normalization
    tmm.func <- function(tab){
      tab <- as.matrix(tab)
      tab_dge <- DGEList(counts=tab)
      tab_tmm_dge <- calcNormFactors(tab_dge, method="TMM")
      tab_norm <- cpm(tab_tmm_dge)
      as.data.frame(tab_norm)
    }
    # TMM normalization
    norm_p1 <- tmm.func(p1)
    norm_p2 <- tmm.func(merge.func(p1,p2))[,colnames(p2)]
    # let p2 have the same genes as p1 
    final_p1 <- norm_p1
    final_p2 <- norm_p2[rownames(p1),]
  }
  
  # TMMwsp, for samples, need to choose a reference first
  # TMM with singleton pairing, an alternative for highly sparse data
  # performed normalization on the training data, and then performed addon normalization of the test data onto the training data, 
  # to ensure that the normalization of the training data does not in any sense depend on the testing data.
  if(norm_method=="TMMwsp"){
    require(edgeR)
    # function for TMMwsp normalization
    tmmwsp.func <- function(tab){
      tab <- as.matrix(tab)
      tab_dge <- DGEList(counts=tab)
      tab_tmm_dge <- calcNormFactors(tab_dge, method="TMMwsp")
      tab_norm <- cpm(tab_tmm_dge)
      as.data.frame(tab_norm)
    }
    # TMM normalization
    norm_p1 <- tmmwsp.func(p1)
    norm_p2 <- tmmwsp.func(merge.func(p1,p2))[rownames(p1),colnames(p2)]
    # let p2 have the same genes as p1 
    final_p1 <- norm_p1
    final_p2 <- norm_p2
  }
  
  
  # RLE+, for samples
  # add a pseudo count 1 to zero values
  # performed normalization on the training data, and then performed addon normalization of the test data onto the training data, 
  # to ensure that the normalization of the training data does not in any sense depend on the testing data.
  if(norm_method=="RLE+"){
    require(DESeq2)
    # replaced all the 0 abundances with 1 (need integers as input)
    p1[p1==0] <- 1
    p2[p2==0] <- 1
    # function for RLE normalization
    rle.func <- function(tab){
      metadata <- data.frame(class=factor(1:ncol(tab)))
      tab_dds <- DESeqDataSetFromMatrix(countData=tab,colData=metadata,design=~class) 
      tab_rle_dds <- estimateSizeFactors(tab_dds)
      tab_norm <- counts(tab_rle_dds, normalized=TRUE)
      as.data.frame(tab_norm)
    }
    # RLE normalization
    norm_p1 <- rle.func(p1)
    norm_p2 <- rle.func(merge.func(p1,p2))[rownames(p1),colnames(p2)]
    # let p2 have the same genes as p1 
    final_p1 <- norm_p1
    final_p2 <- norm_p2
  }
  
  # RLE_poscounts, for samples
  # use the non-zero counts to calculate the geometric mean (type="poscounts")
  # performed normalization on the training data, and then performed addon normalization of the test data onto the training data, 
  # to ensure that the normalization of the training data does not in any sense depend on the testing data.
  if(norm_method=="RLE_poscounts"){
    require(DESeq2)
    # function for RLE with poscounts estimator normalization
    rle.poscounts.func <- function(tab){
      metadata <- data.frame(class=factor(1:ncol(tab)))
      tab_dds <- DESeqDataSetFromMatrix(countData=tab,colData=metadata,design=~class) 
      tab_rle_dds <- estimateSizeFactors(tab_dds,type="poscounts")
      tab_norm <- counts(tab_rle_dds, normalized=TRUE)
      as.data.frame(tab_norm)
    }
    # RLE normalization
    norm_p1 <- rle.poscounts.func(p1)
    norm_p2 <- rle.poscounts.func(merge.func(p1,p2))[rownames(p1),colnames(p2)]
    # let p2 have the same genes as p1 
    final_p1 <- norm_p1
    final_p2 <- norm_p2
  }
  
  # GMPR, for samples
  # switching the two steps in RLE(DESeq2) normalization
  # performed normalization on the training data, and then performed addon normalization of the test data onto the training data, 
  # to ensure that the normalization of the training data does not in any sense depend on the testing data.
  if(norm_method=="GMPR"){
    require(GUniFrac)
    # function for GMPR normalization
    gmpr.func <- function(tab){
      gmpr_size_factor <- GMPR(tab)
      tab_norm <- as.data.frame(t(t(tab)/gmpr_size_factor))
      as.data.frame(tab_norm)
    }
    # GMPR normalization
    norm_p1 <- gmpr.func(p1)
    norm_p2 <- gmpr.func(p2)
    #norm_p2 <- gmpr.func(merge.func(p1,p2))[rownames(p1),colnames(p2)]
    # let p2 have the same genes as p1 
    #final_p1 <- norm_p1
    #final_p2 <- norm_p2
    merged <- merge.func(norm_p1,norm_p2)
    final_p1 <- merged[rownames(norm_p1),colnames(norm_p1)]
    final_p2 <- merged[rownames(norm_p1),colnames(norm_p2)]
  }
  
  # CLR+, for samples
  # based on TSS normalized data, add a pseudo count 0.65*minimum to zero values
  if(norm_method=="CLR+"){
    require(compositions)
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    norm_p2 <- as.data.frame(apply(p2,2,function(x) x/sum(x)))
    # replaced all the 0 abundances with 0.65 times minimum non-zero abundance
    norm_p1[norm_p1==0] <- min(norm_p1[norm_p1!=0])*0.65
    norm_p2[norm_p2==0] <- min(norm_p2[norm_p2!=0])*0.65
    # clr transformation
    trans_p1 <- as.data.frame(apply(norm_p1,2,function(x) clr(x)))
    trans_p2 <- as.data.frame(apply(norm_p2,2,function(x) clr(x)))
    # let p2 have the same genes as p1 
    merged <- merge.func(trans_p1,trans_p2)
    final_p1 <- merged[rownames(trans_p1),colnames(trans_p1)]
    final_p2 <- merged[rownames(trans_p1),colnames(trans_p2)]
  }
  
  # CLR_poscounts, for samples
  # based on TSS normalized data, use the positive counts only
  if(norm_method=="CLR_poscounts"){
    require(compositions)
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    norm_p2 <- as.data.frame(apply(p2,2,function(x) x/sum(x)))
    # clr transformation
    trans_p1 <- as.data.frame(apply(norm_p1,2,function(x) clr(x)))
    trans_p2 <- as.data.frame(apply(norm_p2,2,function(x) clr(x)))
    # let p2 have the same genes as p1 
    merged <- merge.func(trans_p1,trans_p2)
    final_p1 <- merged[rownames(trans_p1),colnames(trans_p1)]
    final_p2 <- merged[rownames(trans_p1),colnames(trans_p2)]
  }
  
  # logcpm, for samples
  if(norm_method=="logcpm"){
    require(edgeR)
    # replaced all the 0 abundances with 1 (need integers as input)
    p1[p1==0] <- 1
    p2[p2==0] <- 1
    # logcpm normalization
    norm_p1 <- as.data.frame(cpm(p1,log=TRUE))
    norm_p2 <- as.data.frame(cpm(p2,log=TRUE))
    # let p2 have the same genes as p1 
    merged <- merge.func(norm_p1,norm_p2)
    final_p1 <- merged[rownames(norm_p1),colnames(norm_p1)]
    final_p2 <- merged[rownames(norm_p1),colnames(norm_p2)]
  }
  
  # LOG, for genes
  # based on TSS normalized data, add a pseudo count 0.65*minimum to zero values
  if(norm_method=="LOG"){
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    norm_p2 <- as.data.frame(apply(p2,2,function(x) x/sum(x)))
    # replaced all the 0 abundances with 0.65 times minimum non-zero abundance
    norm_p1[norm_p1==0] <- min(norm_p1[norm_p1!=0])*0.65
    norm_p2[norm_p2==0] <- min(norm_p2[norm_p2!=0])*0.65
    # log transformation
    trans_p1 <- log(norm_p1)
    trans_p2 <- log(norm_p2)
    # let p2 have the same genes as p1 
    merged <- merge.func(trans_p1,trans_p2)
    final_p1 <- merged[rownames(trans_p1),colnames(trans_p1)]
    final_p2 <- merged[rownames(trans_p1),colnames(trans_p2)]
  }
  
  # AST, for genes
  # based on TSS normalized data
  if(norm_method=="AST"){
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    norm_p2 <- as.data.frame(apply(p2,2,function(x) x/sum(x)))
    # AST transformation
    trans_p1 <- asin(sqrt(norm_p1))
    trans_p2 <- asin(sqrt(norm_p2))
    # let p2 have the same genes as p1 
    merged <- merge.func(trans_p1,trans_p2)
    final_p1 <- merged[rownames(trans_p1),colnames(trans_p1)]
    final_p2 <- merged[rownames(trans_p1),colnames(trans_p2)]
  }
  
  # STD, for genes
  # based on TSS normalized data, substract the mean and devided by the standard deviation
  # performed transformation on the training data, and then performed addon transformation of the test data onto the training data, 
  # to ensure that the transformation of the training data does not in any sense depend on the testing data.
  if(norm_method=="STD"){
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    norm_p2 <- as.data.frame(apply(p2,2,function(x) x/sum(x)))
    # function for STD
    std.func <- function(tab){
      tab_trans <- as.data.frame(t(apply(tab,1,function(x) (x-mean(x))/sd(x))))
      tab_trans
    }
    # STD transformation
    trans_p1 <- std.func(norm_p1)
    trans_p2 <- std.func(merge.func(norm_p1,norm_p2))[,colnames(norm_p2)]
    # let p2 have the same genes as p1 
    final_p1 <- trans_p1
    final_p2 <- trans_p2[rownames(trans_p1),]
  }
  
  # rank, for genes
  # based on TSS normalization data
  # performed transformation on the training data, and then performed addon transformation of the test data onto the training data, 
  # to ensure that the transformation of the training data does not in any sense depend on the testing data.
  if(norm_method=="rank"){
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    norm_p2 <- as.data.frame(apply(p2,2,function(x) x/sum(x)))
    # function for rank transformation
    rank.func <- function(tab){
      tab <- as.matrix(tab)
      # a small noise term is added before data transformation to handle the ties
      noise <- matrix(rnorm(nrow(tab)*ncol(tab),mean=0,sd=10^-10),nrow=nrow(tab),ncol=ncol(tab))
      tab_noised <- tab+noise
      tab_trans <- as.data.frame(t(apply(tab_noised,1,rank)))
      tab_trans
    }
    # rank transformation
    trans_p1 <- rank.func(norm_p1)
    trans_p2 <- rank.func(merge.func(norm_p1,norm_p2))[,colnames(norm_p2)]
    # let p2 have the same genes as p1 
    final_p1 <- trans_p1
    final_p2 <- trans_p2[rownames(trans_p1),]
  }
  
  # blom, for genes
  # based on TSS normalization data
  # performed transformation on the training data, and then performed addon transformation of the test data onto the training data, 
  # to ensure that the transformation of the training data does not in any sense depend on the testing data.
  if(norm_method=="blom"){
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    norm_p2 <- as.data.frame(apply(p2,2,function(x) x/sum(x)))
    # function for blom transformation
    blom.func <- function(tab){
      tab <- as.matrix(tab)
      # a small noise term is added before data transformation to handle the ties
      noise <- matrix(rnorm(nrow(tab)*ncol(tab),mean=0,sd=10^-10),nrow=nrow(tab),ncol=ncol(tab))
      tab_noised <- tab+noise
      c <- 3/8
      tab_trans <- as.data.frame(t(apply(tab_noised,1,function(x) qnorm((rank(x)-c)/(ncol(tab_noised)-2*c+1)))))
      as.data.frame(tab_trans)
    }
    # blom transformation
    trans_p1 <- blom.func(norm_p1)
    trans_p2 <- blom.func(merge.func(norm_p1,norm_p2))[,colnames(norm_p2)]
    # let p2 have the same genes as p1 
    final_p1 <- trans_p1
    final_p2 <- trans_p2[rownames(trans_p1),]
  }
  
  # VST, for genes
  if(norm_method=="VST"){
    require(DESeq2)
    # replaced all the 0 abundances with 1 (need integers as input)
    p1[p1==0] <- 1
    p2[p2==0] <- 1
    # VST transformation
    trans_p1 <- as.data.frame(varianceStabilizingTransformation(as.matrix(p1),blind=FALSE,fitType="mean"))
    trans_p2 <- as.data.frame(varianceStabilizingTransformation(as.matrix(merge.func(p1,p2)),blind=FALSE,fitType="mean"))[,colnames(p2)]
    # let p2 have the same genes as p1 
    final_p1 <- trans_p1
    final_p2 <- trans_p2[rownames(trans_p1),]
  }
  
  # NPN, for genes
  # based on TSS normalized data
  # performed transformation on the training data, and then performed addon transformation of the test data onto the training data, 
  # to ensure that the transformation of the training data does not in any sense depend on the testing data.
  if(norm_method=="NPN"){
    require(huge)
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    norm_p2 <- as.data.frame(apply(p2,2,function(x) x/sum(x)))
    # NPN transformation
    trans_p1 <- as.data.frame(t(huge.npn(t(norm_p1),npn.func="truncation")))
    trans_p2 <- as.data.frame(t(huge.npn(t(merge.func(norm_p1,norm_p2)),npn.func="truncation")))[,colnames(norm_p2)]  
    # let p2 have the same genes as p1 
    final_p1 <- trans_p1
    final_p2 <- trans_p2[rownames(trans_p1),]
  }
  
  # QN, batch correction
  # based on LOG transformed TSS normalized data, using traning as reference 
  if(norm_method=="QN"){
    require(preprocessCore)
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    norm_p2 <- as.data.frame(apply(p2,2,function(x) x/sum(x)))
    # replaced all the 0 abundances with 0.65 times minimum non-zero abundance
    norm_p1[norm_p1==0] <- min(norm_p1[norm_p1!=0])*0.65
    norm_p2[norm_p2==0] <- min(norm_p2[norm_p2!=0])*0.65
    # log transformation
    trans_p1 <- log(norm_p1)
    trans_p2 <- log(norm_p2)
    # reference quantiles, testing set
    ref_quantiles <- normalize.quantiles.determine.target(x=as.matrix(trans_p1))
    # quantile normalize norm_p1 and norm_p2, using the ref_quantiles as reference
    correct_p1 <- normalize.quantiles(as.matrix(trans_p1))
    dimnames(correct_p1) <- dimnames(p1)
    merged <- merge.func(trans_p1,trans_p2)
    norm_merged <- normalize.quantiles.use.target(as.matrix(merged),target=ref_quantiles)
    dimnames(norm_merged) <- dimnames(merged)
    correct_p2 <- as.data.frame(norm_merged[,colnames(p2)])
    # let p2 have the same genes as p1 
    final_p1 <- correct_p1
    final_p2 <- correct_p2[rownames(correct_p1),]
  }
  
  # FSQN, batch correction
  # based on LOG transformed TSS normalized data
  if(norm_method=="FSQN"){
    require(FSQN)
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    norm_p2 <- as.data.frame(apply(p2,2,function(x) x/sum(x)))
    # replaced all the 0 abundances with 0.65 times minimum non-zero abundance
    norm_p1[norm_p1==0] <- min(norm_p1[norm_p1!=0])*0.65
    norm_p2[norm_p2==0] <- min(norm_p2[norm_p2!=0])*0.65
    # log transformation
    trans_p1 <- log(norm_p1)
    trans_p2 <- log(norm_p2)
    # FSQN
    correct_p1 <- trans_p1
    correct_p2 <- merge.func(trans_p1,trans_p2)[rownames(trans_p1),colnames(trans_p2)]
    correct_p2 <- t(quantileNormalizeByFeature(matrix_to_normalize=as.matrix(t(correct_p2)),
                                               target_distribution_matrix=as.matrix(t(correct_p1))))
    # let p2 have the same genes as p1 
    final_p1 <- correct_p1
    final_p2 <- correct_p2
  }
  
  # bmc, batch correction
  # based on LOG transformed TSS normalized data
  if(norm_method=="BMC"){
    require(pamr)
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    norm_p2 <- as.data.frame(apply(p2,2,function(x) x/sum(x)))
    # replaced all the 0 abundances with 0.65 times minimum non-zero abundance
    norm_p1[norm_p1==0] <- min(norm_p1[norm_p1!=0])*0.65
    norm_p2[norm_p2==0] <- min(norm_p2[norm_p2!=0])*0.65
    # log transformation
    trans_p1 <- log(norm_p1)
    trans_p2 <- log(norm_p2)
    # bmc correction
    merged <- merge.func(trans_p1,trans_p2)
    batch_factor <- factor(c(rep(1,ncol(trans_p1)),rep(2,ncol(trans_p2))))
    correct_merged <- as.data.frame(pamr.batchadjust(list(x=as.matrix(merged),batchlabels=batch_factor))$x)
    # let p2 have the same genes as p1 
    final_p1 <- correct_merged[rownames(trans_p1),colnames(trans_p1)]
    final_p2 <- correct_merged[rownames(trans_p1),colnames(trans_p2)]
  }
  
  # limma, batch correction
  # based on LOG transformed TSS normalized data
  if(norm_method=="limma"){
    require(limma)
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    norm_p2 <- as.data.frame(apply(p2,2,function(x) x/sum(x)))
    # replaced all the 0 abundances with 0.65 times minimum non-zero abundance
    norm_p1[norm_p1==0] <- min(norm_p1[norm_p1!=0])*0.65
    norm_p2[norm_p2==0] <- min(norm_p2[norm_p2!=0])*0.65
    # log transformation
    trans_p1 <- log(norm_p1)
    trans_p2 <- log(norm_p2)
    #limma
    merged <- merge.func(trans_p1,trans_p2)
    batch_factor <- factor(c(rep(1,ncol(trans_p1)),rep(2,ncol(trans_p2))))
    correct_merged <- as.data.frame(removeBatchEffect(merged, batch=batch_factor))
    # let p2 have the same genes as p1 
    final_p1 <- correct_merged[rownames(trans_p1),colnames(trans_p1)]
    final_p2 <- correct_merged[rownames(trans_p1),colnames(trans_p2)]
  }
  
  # combat, batch correction
  # based on LOG transformed TSS normalized data, using traning as reference 
  if(norm_method=="combat"){
    require(sva)
    # TSS normalization
    norm_p1 <- as.data.frame(apply(p1,2,function(x) x/sum(x)))
    norm_p2 <- as.data.frame(apply(p2,2,function(x) x/sum(x)))
    # replaced all the 0 abundances with 0.65 times minimum non-zero abundance
    norm_p1[norm_p1==0] <- min(norm_p1[norm_p1!=0])*0.65
    norm_p2[norm_p2==0] <- min(norm_p2[norm_p2!=0])*0.65
    # log transformation
    trans_p1 <- log(norm_p1)
    trans_p2 <- log(norm_p2)
    # combat
    merged <- merge.func(trans_p1,trans_p2)
    batch_factor <- factor(c(rep(1,ncol(trans_p1)),rep(2,ncol(trans_p2))))
    correct_merged <- as.data.frame(ComBat(merged, batch=batch_factor,ref.batch=1))
    # let p2 have the same genes as p1 
    final_p1 <- correct_merged[rownames(trans_p1),colnames(trans_p1)]
    final_p2 <- correct_merged[rownames(trans_p1),colnames(trans_p2)]
  }
  
  # conqur_trn, batch correction
  # directly worked on the count table, using training as reference
  if(norm_method=="conqur"){
    require(ConQuR)
    require(foreach)
    merged <- merge.func(p1,p2)
    batch_factor <- factor(c(rep(1,ncol(p1)),rep(2,ncol(p2))))
    covariates <- data.frame(covariate=rep(1,ncol(merged)))   # no additional information could be used for each dataset
    correct_merged <- ConQuR(tax_tab=as.data.frame(t(merged)),batchid=batch_factor,batch_ref="1",covariates=covariates,simple_match=T)
    correct_merged <- as.data.frame(t(correct_merged))
    ## normalize the corrected count data
    #norm_merged <- as.data.frame(apply(correct_merged,2,function(x) x/sum(x)))
    # let p2 have the same genes as p1 
    final_p1 <- correct_merged[rownames(p1),colnames(p1)]
    final_p2 <- correct_merged[rownames(p1),colnames(p2)]
  }
  
  # return
  list(final_p1,final_p2)
}


#======================================================================================================================#
### function for predictions of simulated data
sim.pred.func  <- function(trn,tst,pred_method){
  
  ## rare taxa filtering, remove genes with >90% zeros
  #filter_trn <- trn[!(rowSums(trn==0) > ncol(trn)*0.9),]
  #filter_tst <- tst[!(rowSums(tst==0) > ncol(tst)*0.9),]
  
  # transform the trn and tst, making samples in rows and genes in columns, add the disease status for train
  pred_trn <- as.data.frame(t(trn))
  pred_tst <- as.data.frame(t(tst))
  pred_trn$status <- as.factor(gsub("p1_([a-z]+)_[0-9]+","\\1",rownames(pred_trn)))
  
  ### random forest
  if(pred_method=="rf"){
    ctrl <- trainControl(method="repeatedcv", number=10, search="grid", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions=TRUE)
    rf <- train(status~., data=pred_trn, method="rf", num.trees=1000, metric="ROC",trControl=ctrl)
    rf_pred <- predict(rf, pred_tst, type="prob")
    rf_pred$real <- gsub("p2_([a-z]+)_[0-9]+","\\1",rownames(pred_tst))
    auc_value <- roc(rf_pred$real,rf_pred$control,levels=c("control", "case"), direction=">")$auc
  }
  
  ### random forest
  if(pred_method=="ranger"){
    ctrl <- trainControl(method="repeatedcv", number=10, search="grid", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions=TRUE)
    rf <- train(status~., data=pred_trn, method="ranger", num.trees=1000, metric="ROC",trControl=ctrl)
    rf_pred <- predict(rf, pred_tst, type="prob")
    rf_pred$real <- gsub("p2_([a-z]+)_[0-9]+","\\1",rownames(pred_tst))
    # AUC
    auc_value <- roc(rf_pred$real,rf_pred$control,levels=c("control", "case"), direction=">")$auc
    # confusion matrix
    TN <- sum(rf_pred$real=="ctrl"&rf_pred$pred=="ctrl")  # True Negatives
    FP <- sum(rf_pred$real=="ctrl"&rf_pred$pred=="case")  # False Positives
    FN <- sum(rf_pred$real=="case"&rf_pred$pred=="ctrl")  # False Negatives
    TP <- sum(rf_pred$real=="case"&rf_pred$pred=="case")  # True Positives
    # accuracy, sensitivity, and specificity
    accuracy <- (TN+TP)/sum(TN+FP+FN+TP)
    sensitivity <- TP/(TP+FN)
    specificity <- TN/(TN+FP)
  }
  
  #### lasso (logistic regression)
  #if(pred_method=="lr"){
  #  parameters <- seq(-4,4,0.1)
  #  ctrl <- trainControl(method="repeatedcv", number=10, search="grid", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions=TRUE)
  #  lr <- train(status~., data=pred_trn, method="glmnet", family="binomial", metric="ROC",
  #              tuneGrid=expand.grid(alpha=1, lambda=parameters),trControl=ctrl)
  #  lr_pred <- predict(lr, pred_tst, type="prob")
  #  lr_pred$real <- gsub("p2_([a-z]+)_[0-9]+","\\1",rownames(pred_tst))
  #  auc_value <- roc(lr_pred$real,lr_pred$control,levels=c("control", "case"), direction=">")$auc
  #}
  
  #### svm (linear kernal)
  #if(pred_method=="svm"){
  #  ctrl <- trainControl(method="repeatedcv", number=10, search="grid", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions=TRUE)
  #  svm <- train(status~., data=pred_trn, method="svmLinear2", metric="ROC",trControl=ctrl)
  #  svm_pred <- predict(svm, pred_tst, type="prob")
  #  svm_pred$real <- gsub("p2_([a-z]+)_[0-9]+","\\1",rownames(pred_tst))
  #  auc_value <- roc(svm_pred$real,svm_pred$control,levels=c("control", "case"), direction=">")$auc
  #}
  
  # return
  return(c(auc_value,accuracy,sensitivity,specificity))
}


#======================================================================================================================#
### function for predictions of real data
pred.func  <- function(trn,tst,meta,pred_method){

  ## rare taxa filtering, remove genes with >90% zeros
  #filter_trn <- trn[!(rowSums(trn==0) > ncol(trn)*0.9),]
  #filter_tst <- tst[!(rowSums(tst==0) > ncol(tst)*0.9),]

  # transform the train and test data, add the phenotype for train data
  pred_trn <- as.data.frame(t(trn))
  pred_tst <- as.data.frame(t(tst))
  pred_trn$status <- as.factor(meta[rownames(pred_trn),"phenotype"])
  
  # GMPR had the problem of NA values, remove samples with all NAs
  pred_trn <- pred_trn[!is.na(pred_trn[,1]),]
  pred_tst <- pred_tst[!is.na(pred_tst[,1]),]
  
  ### random forest
  if(pred_method=="rf"){
    ctrl <- trainControl(method="repeatedcv", number=10, search="grid", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions=TRUE)
    rf <- train(status~., data=pred_trn, method="rf", num.trees=1000, metric="ROC",trControl=ctrl)
    rf_pred <- predict(rf, pred_tst, type="prob")
    rf_pred$real <- meta[rownames(pred_tst),"phenotype"]
    auc_value <- roc(rf_pred$real,rf_pred$control,levels=c("control", "case"), direction=">")$auc
  }
  
  ### random forest, realize in ranger
  if(pred_method=="ranger"){
    ctrl <- trainControl(method="repeatedcv", number=10, search="grid", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions=TRUE)
    rf <- train(status~., data=pred_trn, method="ranger", num.trees=1000, metric="ROC",trControl=ctrl)
    rf_pred <- predict(rf, pred_tst, type="prob")
    rf_pred$real <- meta[rownames(pred_tst),"phenotype"]
    # AUC
    auc_value <- roc(rf_pred$real,rf_pred$control,levels=c("control", "case"), direction=">")$auc
    # confusion matrix
    TN <- sum(rf_pred$real=="ctrl"&rf_pred$pred=="ctrl")  # True Negatives
    FP <- sum(rf_pred$real=="ctrl"&rf_pred$pred=="case")  # False Positives
    FN <- sum(rf_pred$real=="case"&rf_pred$pred=="ctrl")  # False Negatives
    TP <- sum(rf_pred$real=="case"&rf_pred$pred=="case")  # True Positives
    # accuracy, sensitivity, and specificity
    accuracy <- (TN+TP)/sum(TN+FP+FN+TP)
    sensitivity <- TP/(TP+FN)
    specificity <- TN/(TN+FP)
  }
  
  #### lasso (logistic regression)
  #if(pred_method=="lr"){
  #  parameters <- seq(-4,4,0.1)
  #  ctrl <- trainControl(method="repeatedcv", number=10, search="grid", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions=TRUE)
  #  lr <- train(status~., data=pred_trn, method="glmnet", family="binomial", metric="ROC",
  #              tuneGrid=expand.grid(alpha=1, lambda=parameters),trControl=ctrl)
  #  lr_pred <- predict(lr, pred_tst, type="prob")
  #  lr_pred$real <- meta[rownames(pred_tst),"phenotype"]
  #  auc_value <- roc(lr_pred$real,lr_pred$control,levels=c("control", "case"), direction=">")$auc
  #}
  
  #### svm (linear kernal)
  #if(pred_method=="svm"){
  #  ctrl <- trainControl(method="repeatedcv", number=10, search="grid", summaryFunction=twoClassSummary, classProbs=TRUE, savePredictions=TRUE)
  #  svm <- train(status~., data=pred_trn, method="svmLinear2", metric="ROC",trControl=ctrl)
  #  svm_pred <- predict(svm, pred_tst, type="prob")
  #  svm_pred$real <- meta[rownames(pred_tst),"phenotype"]
  #  auc_value <- roc(svm_pred$real,svm_pred$control,levels=c("control", "case"), direction=">")$auc
  #}
  
  # return
  return(c(auc_value,accuracy,sensitivity,specificity))
  
}




