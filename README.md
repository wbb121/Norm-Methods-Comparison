# Comparison of the effectiveness of different normalization methods for metagenomic cross-study phenotype prediction under heterogeneity

This repository contains R scripts for running the analysis described in our manuscript. 

### [helper.R](https://github.com/wbb121/Norm-Methods-Comparison/blob/main/helper.R):

Including functions to merge count tables, obtain count tables from curatedMetagenomicData, perform PERMANOVA analysis, and normalize the data.

### [real_data_obtain.R](https://github.com/wbb121/Norm-Methods-Comparison/blob/main/real_data_obtain.R):

Obtain metadata and count tables of datasets related to CRC (colorectal cancer) and IBD (inflammatory bowel disease) from curatedMetagenomicData.

### [sim_scenario1.R](https://github.com/wbb121/Norm-Methods-Comparison/blob/main/sim_scenario1.R)：

Validate the performance of different normalization methods on binary phenotype prediction in simulation scenario 1.

Use the following command to run the script:

```shell
Rscript sim_scenario1.R sample_size library_size num_genes count1 count2 norm_method pred_method
```

+ **sample_size**: sample size of one simulated dataset
+ **library_size**: library size of each sample
+ **num_genes**: number of disease related genes
+ **count1**: simulation template 1
+ **count2**: simulation template 2
+ **norm_method**: normalization method
+ **pred_method**: prediction method

### [sim_scenario2.R](https://github.com/wbb121/Norm-Methods-Comparison/blob/main/sim_scenario2.R)：

Validate the performance of different normalization methods on binary phenotype prediction in simulation scenario 2.

Use the following command to run the script:

```shell
Rscript sim_scenario2.R sample_size library_size num_genes count norm_method pred_method
```

+ **sample_size**: sample size of one simulated dataset
+ **library_size**: library size of each sample
+ **num_genes**: number of disease related genes
+ **count**: simulation template 
+ **norm_method**: normalization method
+ **pred_method**: prediction method

### [sim_scenario3.R](https://github.com/wbb121/Norm-Methods-Comparison/blob/main/sim_scenario3.R)：

Validate the performance of different normalization methods on binary phenotype prediction in simulation scenario 3.

Use the following command to run the script:

```shell
Rscript sim_scenario3.R sample_size library_size num_genes count norm_method pred_method
```

+ **sample_size**: sample size of one simulated dataset
+ **library_size**: library size of each sample
+ **num_genes**: number of disease related genes
+ **count**: simulation template 
+ **norm_method**: normalization method
+ **pred_method**: prediction method

### [real_data_analysis.R](https://github.com/wbb121/Norm-Methods-Comparison/blob/main/real_data_analysis.R):

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

### [pred_res_summ.R](https://github.com/wbb121/Norm-Methods-Comparison/blob/main/pred_res_summ.R):

Summarize the prediction results from simulation and real data.

### [figures.R](https://github.com/wbb121/Norm-Methods-Comparison/blob/main/figures.R): 

Draw the figures in the manuscript.



