########################################################################################################################
### obtain CRC and IBD data from curatedMetagenomicData
### 2023/06/05
########################################################################################################################

setwd("/home/wangbb/normalization_comparison")

### packages
library(curatedMetagenomicData)

### load the helper function
source("helper.R")

#================================================================================================================#
### datasets related to CRC ###
crc_studies <- unique(sampleMetadata[sampleMetadata$study_condition=="CRC","study_name"])
crc_studies <- crc_studies[!is.na(crc_studies)]


### obtain the metadata of control and CRC samples from above mentioned datasets
crc_meta <- sampleMetadata[sampleMetadata$study_name%in%crc_studies & sampleMetadata$study_condition%in%c("control","CRC"),]
rownames(crc_meta) <- crc_meta$sample_id


### obtain the count table of control and CRC samples from above mentioned datasets
crc_count <- data.frame()
for(i in crc_studies){
  count_table <- obtain.count.table(study_name=i,metadata=crc_meta)
  crc_count <- merge.func(crc_count,count_table)
}


### crc data processing
table(crc_meta$study_name,crc_meta$study_condition)
# remove "HanniganGD_2017" due to its small sample size
crc_meta <- crc_meta[crc_meta$study_name!="HanniganGD_2017",]
crc_count <- crc_count[,colnames(crc_count)%in%crc_meta$sample_id]
crc_count <- crc_count[rowSums(crc_count)>0,]
# remove "ThomasAM_2019_c" due to the overlap with "YachidaS_2019"
crc_meta <- crc_meta[crc_meta$study_name!="ThomasAM_2019_c",]
crc_count <- crc_count[,colnames(crc_count)%in%crc_meta$sample_id]
crc_count <- crc_count[rowSums(crc_count)>0,]
# combine "ThomasAM_2018a" and "ThomasAM_2018b" to "ThomasAM_2018"
crc_meta$study_name <- sub("ThomasAM_2018a","ThomasAM_2018",crc_meta$study_name)
crc_meta$study_name <- sub("ThomasAM_2018b","ThomasAM_2018",crc_meta$study_name)
# simplify the study_name
crc_meta$study_name_simplified <- gsub("([A-Z][a-z]+)[A-Z]+_[0-9]+", "\\1", crc_meta$study_name)


### characteristics of CRC datasets
crc_studies <- unique(crc_meta$study_name) # 8 CRC datasets:"FengQ_2015","GuptaA_2019","ThomasAM_2018","VogtmannE_2016","WirbelJ_2018","YachidaS_2019","YuJ_2015","ZellerG_2014"
crc_character <- data.frame(study_name=NA,country=NA,control_num=NA,case_num=NA,dna_extraction_kit=NA,
                            sequencing_platform=NA,reference=NA)
for(i in 1:length(crc_studies)){
  df <- crc_meta[crc_meta$study_name==crc_studies[i],]
  crc_character[i,"study_name"] <- crc_studies[i]
  crc_character[i,"country"] <- paste(unique(df$country),collapse="/")
  crc_character[i,"control_num"] <- sum(df$study_condition=="control")
  crc_character[i,"case_num"] <- sum(df$study_condition=="CRC")
  crc_character[i,"dna_extraction_kit"] <- paste(unique(df$DNA_extraction_kit),collapse="/")
  crc_character[i,"sequencing_platform"] <- paste(unique(df$sequencing_platform),collapse="/")
  crc_character[i,"reference"] <- paste(unique(df$PMID),collapse="/")
}
write.csv(crc_character,"tables/crc_data_characteristics.csv")


### save the data
saveRDS(crc_meta,"data/crc_meta.rds")
saveRDS(crc_count,"data/crc_count.rds")


#================================================================================================================#
### datasets related to IBD ###
ibd_studies <- unique(sampleMetadata[sampleMetadata$study_condition=="IBD","study_name"])
ibd_studies <- ibd_studies[!is.na(ibd_studies)]


### obtain the metadata of control and IBD samples from above mentioned datasets
ibd_meta <- sampleMetadata[sampleMetadata$study_name%in%c(ibd_studies,"LifeLinesDeep_2016"),]
ibd_meta <- ibd_meta[ibd_meta$study_condition=="control"|ibd_meta$study_condition=="IBD",]
rownames(ibd_meta) <- ibd_meta$sample_id


### obtain the count table of control and IBD samples from above mentioned datasets
ibd_count <- data.frame()
for(i in c(ibd_studies,"LifeLinesDeep_2016")){
  count_table <- obtain.count.table(study_name=i,metadata=ibd_meta)
  ibd_count <- merge.func(ibd_count,count_table)
}


### ibd data processing
table(ibd_meta$study_name,ibd_meta$study_condition)
# remove "LiJ_2014" due to it small sample size
ibd_meta <- ibd_meta[ibd_meta$study_name!="LiJ_2014",]
ibd_count <- ibd_count[,colnames(ibd_count)%in%ibd_meta$sample_id]
ibd_count <- ibd_count[rowSums(ibd_count)>0,]
# combine "LifeLinesDeep_2016" and "VilaAV_2018" to "VilaAV_2018"
ibd_meta$study_name <- gsub("LifeLinesDeep_2016","VilaAV_2018",ibd_meta$study_name)
# one sample in "HMP_2019_ibdmdb" had no abundance profiles, remove it from ibd_meta
ibd_meta <- ibd_meta[ibd_meta$sample_id%in%colnames(ibd_count),]
# simplify the study_name
ibd_meta$study_name_simplified <- gsub("([A-Z][a-z]+)[A-Z]+_[0-9]+", "\\1", ibd_meta$study_name)
ibd_meta$study_name_simplified <- sub("HMP_2019_ibdmdb","HMP",ibd_meta$study_name_simplified)


# Characteristics of IBD datasets
ibd_studies <- unique(ibd_meta$study_name) # 5 ibd datasets: "HallAB_2017","HMP_2019_ibdmdb","IjazUZ_2017","VilaAV_2018","NielsenHB_2014"
ibd_character <- data.frame(study_name=NA,country=NA,control_num=NA,case_num=NA,subject_num=NA,
                            dna_extraction_kit=NA,sequencing_platform=NA,reference=NA)
for(i in 1:length(ibd_studies)){
  df <- ibd_meta[ibd_meta$study_name==ibd_studies[i],]
  ibd_character[i,"study_name"] <- ibd_studies[i]
  ibd_character[i,"country"] <- paste(unique(df$country),collapse="/")
  ibd_character[i,"control_num"] <- sum(df$study_condition=="control")
  ibd_character[i,"case_num"] <- sum(df$study_condition=="IBD")
  ibd_character[i,"subject_num"] <- length(unique(df$subject_id))
  ibd_character[i,"dna_extraction_kit"] <- paste(unique(df$DNA_extraction_kit),collapse="/")
  ibd_character[i,"sequencing_platform"] <- paste(unique(df$sequencing_platform),collapse="/")
  ibd_character[i,"reference"] <- paste(unique(df$PMID),collapse="/")
}
write.csv(ibd_character,"tables/ibd_data_characteristics.csv")


### save the data
saveRDS(ibd_meta,"data/ibd_meta.rds")
saveRDS(ibd_count,"data/ibd_count.rds")

















