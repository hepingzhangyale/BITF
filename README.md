## Ball Impurity: Measuring Heterogeneity in General Metric Spaces

Menglu Che, Ting Li, Wenliang Pan, Xueqin Wang, Heping Zhang

Data in various domains, such as neuroimaging and network data analysis, often come in complex forms without possessing a Hilbert structure. The complexity necessitates innovative approaches for effective analysis. We propose a novel measure of heterogeneity, ball impurity, which is designed to work with complex non-Euclidean objects. Our approach extends the notion of impurity to general metric spaces, providing a versatile tool for feature selection and tree models. The ball impurity measure exhibits desirable properties, such as the triangular inequality, and is computationally tractable, enhancing its practicality and usefulness. Extensive experiments on synthetic data and real data from the UK Biobank validate the efficacy of our approach in capturing data heterogeneity. Remarkably, our results compare favorably with stateof-the-art methods in metric spaces, highlighting the potential of ball impurity as a valuable tool for addressing complex data analysis tasks.

This repository provides the R codes to employ the ball impurity to select effective splits and construct tree models, as well as reproduce the tables and figures in the paper. 


## Code
---------
Required codes are located under ./code/*, including 
* `BallDistanceVector.cpp` : C++ code to fastly compute the distance matrix of complex data objects
* `BallImpurity.cpp`: C++ code to fastly compute the ball impurity of a set of complex data objects.
  
Note： with the library `Rcpp` the above cpp codes can be directly compiled in R.

* `BImp_bin.R` : algorithmically improved Ball impurity computation based on binary search.
* `biRegressionTree.R` : Ball tree construction
* `SNP_IR_full.R` and `SNP_IR_subsample.R` : computation of the ball impurity reduction of a single SNP using the full sample, or subsampling strategy.


## Instructions to Reproduce the Numerical Demonstrations and Simulations (Tables 1-3, Figure 3,5,6)
--------
### Simulations based on randomly generated data
 `./simulation/1-reproduction_tables1-3.R` reproduces the simulation results in the paper (Tables 1-3, Figures 1-5). 

### Simulation based on real genotype data from UKB
`./simulation/2-real_geno_data_based_simu.R` is a demo of the simulation based on real UKB data in the paper (Figure 6). 

Given data privacy and capacity constraints, we did not put the original genotype data online. The other non-sensitive files required to complete the simulation is provided as `emoji_smile_adverse.csv`.

## Real Data Application
----
Given data privacy and capacity constraints, we did not put the original genotype data from UKB online. A general, demonstrative set of codes is provided in the `\Real_data_application` folder, contains the codes required to perform the real data application in Chapter 5. 

Required input files:
* A cvs file containing all the baseline characteristics, (in our case, named as ukb50595.csv),
* The plink data files (.bed, .bim, and .fam) of the genomic data.
* The phenotype, or complicated data objects files. For the resting state fMRI partial correlaiton matrices, a `txt` file is provided for each subject. But in general, it can be tailored to any data objects, provided as a list in R. 

### Step 1. 

* Before R programming, prune the plink genetic data as needed. For example, the pruning used to obtain the 3740 SNPs in Chromosome 19 for the ball tree in Section 5.4 of the paper used the following filter:

plink --bfile your_input_file --indep-pairwise 500 50 0.2 --out pruned_data

* Run "Initialization.R" to transform the genotype data into the format of a snpMatrix, and compute the distance matrix of the phenotye data. The genotype and phenotype data will be saved as `geno_transed.RData` and `pheno.RData` respectively. The full distance matrice of the phenotype data objects will be saved as `full_dist_mat.RData`.

### Step 2. 
Compute the ball impurity reduction of all the SNPs in the snpMatrix as demonstrated in "Genome-wide-BI.R". The default method is to use the subsampling strategy. 

### Step 3. 
Run "summarize_final.R" to summarize and generate the manhattan plots. 

### Step 4 (Tree construction)
For smaller data sets, `biRegressionTree.R` in the `\code` folder can be directly utilized to obtain the ball tree. For larger data sets, e.g. the ball tree in Section 5.4 of the paper, we recommend iteratively performing Step 2 and 3 to mannually obtain the tree. 
