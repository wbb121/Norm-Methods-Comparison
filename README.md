# Comparison of the effectiveness of different normalization methods for metagenomic cross-study phenotype prediction under heterogeneity

This repository contains R scripts for running the analysis described in our manuscript. 

##### [helper.R](https://github.com/wbb121/Norm-Methods-Comparison/blob/main/helper.R):

Including functions to merge count tables, obtain count tables from curatedMetagenomicData, perform PERMANOVA analysis, and normalize the data.

##### [real_data_obtain.R](https://github.com/wbb121/Norm-Methods-Comparison/blob/main/real_data_obtain.R):

Obtain metadata and count tables of datasets related to CRC (colorectal cancer) and IBD (inflammatory bowel disease) from curatedMetagenomicData.

##### [sim_data_analysis.R](https://github.com/wbb121/Norm-Methods-Comparison/blob/main/sim_data_analysis.R)：

Validate the performance of different normalization methods on binary phenotype prediction in simulated datasets.

Use the following command to run the script:

```shell
Rscript sim_data_analysis.R sim_cluster ed ns ls nd population1 population2 norm_cluster norm_method pred_cluster pred_method
```

+ **sim_cluster**: number of clusters to do the simulation
+ **ed**: disease effect
+ **ns**: sample size of one simulated dataset
+ **ls**: library size of each sample
+ **nd**: number of disease related genes
+ **population1**: simulation template 1
+ **population2**: simulation template 2
+ **norm_cluster**: number of clusters to do the normalization
+ **norm_method**: normalization method
+ **pred_cluster**: number of clusters to do the prediction
+ **pred_method**: prediction method

##### [real_data_analysis.R](https://github.com/wbb121/Norm-Methods-Comparison/blob/main/real_data_analysis.R):

Validate the performance of different normalization methods on binary phenotype prediction in CRC/IBD datasets.

Use the following command to run the script:

```shell
Rscript real_data_analysis.R meta count norm_method pred_cluster pred_method
```

+ **meta**: metadata of all datasets related to a certain disease (CRC/IBD)
+ **count**: count table of all datasets related to a certain disease (CRC/IBD)
+ **norm_method**: normalization method
+ **pred_cluster**: number of clusters to do the prediction
+ **pred_method**: prediction method

##### [pred_res_summ.R](https://github.com/wbb121/Norm-Methods-Comparison/blob/main/pred_res_summ.R):

Summarize the prediction results from simulation and real data.

##### [figures.R](https://github.com/wbb121/Norm-Methods-Comparison/blob/main/figures.R): 

Draw the figures in the manuscript.



