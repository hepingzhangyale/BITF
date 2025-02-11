## Ball Impurity: Measuring Heterogeneity in General Metric Spaces

Menglu Che, Ting Li, Wenliang Pan, Xueqin Wang, Heping Zhang

Data in various domains, such as neuroimaging and network data analysis, often come in complex forms without possessing a Hilbert structure. The complexity necessitates innovative approaches for effective analysis. We propose a novel measure of heterogeneity, ball impurity, which is designed to work with complex non-Euclidean objects. Our approach extends the notion of impurity to general metric spaces, providing a versatile tool for feature selection and tree models. The ball impurity measure exhibits desirable properties, such as the triangular inequality, and is computationally tractable, enhancing its practicality and usefulness. Extensive experiments on synthetic data and real data from the UK Biobank validate the efficacy of our approach in capturing data heterogeneity. Remarkably, our results compare favorably with stateof-the-art methods in metric spaces, highlighting the potential of ball impurity as a valuable tool for addressing complex data analysis tasks.

This repository provides all the R codes needed to reproduce the tables and figures in the paper. 

Helper funcitions
-------
The helper functions required to be sourced are located under ./utility_functions/*, including 
* `BallDistanceVector.cpp`
* `BallImpurity.cpp`
* `BImp_bin.R` : computationally improved Ball impurity computation based on binary search.
* `BImp_subsample.R` : computation of Ball impurity using the subsample strategy, which greatly reduces the computation time.
* `splitDataset.R` : splitting the data set according to feature values.
* `biRegressionTree_Euclidean.R` and `biRegressionTree_NE.R` : Ball tree construction of Euclidean and Non-euclidean response data.
* `SNP_IR_full.R` and `SNP_IR_subsample.R` : computation of the ball impurity reduction of a single SNP using the full sample, or subsampling strategy.
* `simu_rank_matrix.R` : used in simulation for ranking the candidate features by their associations with the response.

Code
----
The R codes needed to reproduce the numerical demonstrations and simulation part in the paper are organized by folders named after each Table and Figure number. The random seeds used in the simulation are provided in the R scripts. All scripts can be used directly. Each code is  Figure 5 shares the same code as Table 3 thus not repeated. 

Given data privacy and capacity constraints, we did not put the original genotype data online. We generated a demo dataset. We generated a demo data set `geno_simu.RData` for demonstration and reproduces our results. The other non-sensitive files required to complete the simulation is provided as `emoji_smile_adverse.csv`.

Real Data Application
----
The `\Real_data_application` folder contains the codes required to perform the real data application in Chapter 5. 

A cvs file containing all the baseline characteristics, specifically including fields "eid", "31-0.0", "34-0.0", "52-0.0", "53-0.0", "53-2.0", "22001-0.0", "22006-0.0", "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", is needed (in our case, named as ukb50595.csv), together with all the plink data files (.bed, .bim, and .fam) of the genomic data of Chromosome 1-22, and the phenotype data of Field 25753. Each subject should have a `.txt` file containing a upper triangular part of a 55x55 matrix. These data can be obtained from UK Biobank via their standard data access procedure. After obtaining the data,  follow the workflow:

1. Run "Initialized_geno_blkwise.R" and "initialize_pheno.R" to save the genotype and phenotye data into files required later. The genotype data are organized into blocks of size 1MB, and the discovery and verification sets will be saved separately. The full distance matrices will be calculated for both the discovery and verification set.

2. For each, compute the ball impurity reduction as demonstrated in "Genome-wide-BI.R". The R script will generate ball impurity reduction results for a block of size 1MB, on both the discovery and verification set.

3. Run "summarize_final.R" to summarize and plot the manhattan plots. 


For the Ball Tree of Chromosome 19, first filter the .bed, .bim, .fam files through plink command: 
plink --bfile your_input_file --indep-pairwise 500 50 0.2 --out pruned_data

Then compute the ball impurity reduction of the 3740 SNPS individually as  demonstrated in "Genome-wide-BI.R". Iteratively repeat the scanning of each SNP to obtain the tree.  
