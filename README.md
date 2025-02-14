## Ball Impurity: Measuring Heterogeneity in General Metric Spaces

Menglu Che, Ting Li, Wenliang Pan, Xueqin Wang, Heping Zhang

Data in various domains, such as neuroimaging and network data analysis, often come in complex forms without possessing a Hilbert structure. The complexity necessitates innovative approaches for effective analysis. We propose a novel measure of heterogeneity, ball impurity, which is designed to work with complex non-Euclidean objects. Our approach extends the notion of impurity to general metric spaces, providing a versatile tool for feature selection and tree models. The ball impurity measure exhibits desirable properties, such as the triangular inequality, and is computationally tractable, enhancing its practicality and usefulness. Extensive experiments on synthetic data and real data from the UK Biobank validate the efficacy of our approach in capturing data heterogeneity. Remarkably, our results compare favorably with stateof-the-art methods in metric spaces, highlighting the potential of ball impurity as a valuable tool for addressing complex data analysis tasks.

This repository provides the R codes needed to reproduce the tables and figures in the paper. 


## Code
---------
Required codes are located under ./code/*, including 
* `BallDistanceVector.cpp` : C++ code to fastly compute the distance matrix of complex data objects
* `BallImpurity.cpp`: C++ code to fastly compute the ball impurity of a set of complex data objects.
Note the above codes requires the library `Rcpp` to compile. 
* `BImp_bin.R` : algorithmically improved Ball impurity computation based on binary search.
* `biRegressionTree.R` : Ball tree construction
* `SNP_IR_full.R` and `SNP_IR_subsample.R` : computation of the ball impurity reduction of a single SNP using the full sample, or subsampling strategy.

## Data
--------
Data used for simulation and real data preprocessing are located under  ./data/*.

### Simulation Based on Real Data
1. `geno_simu.RData`: We generated a demo data set `geno_simu.RData` for demonstration and reproduces our results. 
2. The 64x65 signal matrix (Figure S1 in supplementary material) required to complete the simulation is provided as `emoji_smile_adverse.csv`.

   ### Real UKB Data Application
1. `kinship_3degree.Rdata`: Subjects within 3 degree kinship are provided in the file for the preprocessing of removing relatives. 
2. `num_blk.csv`: In preprocessing of the UKB data, we divided the genotype data into blocks of size 1MB. The number of blocks of each chromosome is provided for convenient dividing. 




## Instructions to Reproduce the Numerical Demonstrations (Figures 1-2, 4)
--------------
Required codes are located under ./numerical_demonstrations/*., which generates all subfigures of Figures 1-2 and 4 of the manuscript.




## Instructions to Perform Simulations
--------
### Step 1
 `./simulation/1-reproduction_tables1-3.R` reproduces the simulation results in the paper (Tables 1-3, Figures 3, 5, or respectively Models 2.1a to 3.3). The model number need to be provided as an input of the wrapper function `reproduction`. For Models 3.1 and 3.2, an outlier ratio number (in percentage) is also needed as a second input parameter. 

### Step 2
`./simulation/2-real_geno_data_based_simu.R` reproduces the results of the simulation based on real UKB data in the paper (Figure 6). 

Given data privacy and capacity constraints, we did not put the original genotype data online. We generated a demo dataset. We generated a demo data set `geno_simu.RData` for demonstration and reproduces our results. The other non-sensitive files required to complete the simulation is provided as `emoji_smile_adverse.csv`.

## Real Data Application
----
The `\Real_data_application` folder contains the codes required to perform the real data application in Chapter 5. 

A cvs file containing all the baseline characteristics, specifically including fields "eid", "31-0.0", "34-0.0", "52-0.0", "53-0.0", "53-2.0", "22001-0.0", "22006-0.0", "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", is needed (in our case, named as ukb50595.csv), together with all the plink data files (.bed, .bim, and .fam) of the genomic data of Chromosome 1-22, and the phenotype data of Field 25753. Each subject should have a `.txt` file containing a upper triangular part of a 55x55 matrix. These data can be obtained from UK Biobank via their standard data access procedure. After obtaining the data,  follow the workflow:

### Step 1. 
Run "Initialization.R" to save the genotype and phenotye data into files required later. The genotype data are organized into blocks of size 1MB, and the discovery and verification sets will be saved separately. The full distance matrices will be calculated for both the discovery and verification set.

### Step 2. 
For each, compute the ball impurity reduction as demonstrated in "Genome-wide-BI.R". The R script will generate ball impurity reduction results for a block of size 1MB, on both the discovery and verification set.

### Step 3. 
Run "summarize_final.R" to summarize and plot the manhattan plots. 


For the Ball Tree of Chromosome 19, first filter the .bed, .bim, .fam files through plink command: 
plink --bfile your_input_file --indep-pairwise 500 50 0.2 --out pruned_data

Then compute the ball impurity reduction of the 3740 SNPS individually as  demonstrated in "Genome-wide-BI.R". Iteratively repeat the scanning of each SNP to obtain the tree. 
