rm(list = ls()); gc(reset = TRUE)
source("code/SNP_IR_full.R")
source("code/SNP_IR_subsample.R")
library(snpStats)

geno.data<-"geno_transed.RData"
load(geno.data)
snpsum <- col.summary(geno)

geno <- t(as(geno, "numeric"))
coding.test <- snpsum[, 8] - snpsum[, 6]
geno[which(coding.test>0), ] <- 2 - geno[which(coding.test>0), ]
geno <- as.data.frame(t(geno))

# impute NAs with MAF
for(j in 1:ncol(geno)){
  geno_imp = sample(c(0,1,2),size = sum(is.na(geno.dis[, j])), prob = unlist(c(snpsum[j,c(8,7,6)])), replace =T)
  geno[which(is.na(geno.dis[,j])),j] = geno_imp
}

rm(list=c("basechar","call.cut","call.keep","coding.test","geno_imp","hwe.keep","HWE.cut.off",,"map.call.keep","map.use","minor","snpsum","snpsum.call.keep"))
n.snp=nrow(map)
###########initialize the BI Reduction results matrix
PG.mat=matrix(rep(0,n.snp),nrow=n.snp)
row.names(PG.mat)<-row.names(map)
#colnames(PG.mat)<-c("discovery","verification")
###########compute the BI reduction on the discovery set
load("../pheno.RData")
load("../full_dist_mat.RData")
for(k in 1:n.snp){
PG.mat[k,1]<-SNP_IR_subsample(pheno,geno[,k],10,10,full_dist_mat)
}
rm(list=c("geno","full_dist_mat","pheno"))

save.image(file=paste0("GWBI_results.RData"))
