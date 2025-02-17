rm(list = ls()); gc(reset = TRUE)
require(snpStats) 
library(igraph)
library(data.table)
wd="" ########please set working directory as needed to store the data to be used in later analysis.
setwd(wd)
maintable=fread("ukb50595.csv")

bulk=read.table("ukb50595.bulk")

load("data/kinship_3degree.Rdata")
kinship_3degree[,1]=as.character(kinship_3degree[,1])
kinship_3degree[,2]=as.character(kinship_3degree[,2])
Kinship_network=simplify(graph_from_edgelist(as.matrix(kinship_3degree),directed=F))

delete_relative <- function(data)
{
  cases_kinship_network=induced_subgraph(Kinship_network, v=V(Kinship_network)$name %in% data$eid)
  S_G=delete.vertices(cases_kinship_network,degree(cases_kinship_network)==0)
  rm(cases_kinship_network)
  cases_kikedout_kinship_id=as.numeric(V(S_G)$name)
  cases_id=setdiff(data$eid,cases_kikedout_kinship_id)
  new_data=data[data$eid %in% cases_id,]
  return(new_data)
}
maintable=maintable[which(maintable$eid %in% bulk$V1),]

################optional, to stratify the sample by age and sex ##############
# colnames(maintable)<-c("eid","sex","birth_year","birth_month","assessment_date","image_date","genetic_sex","ethnicity","pc1","pc2","pc3",
#                        "pc4","pc5")
# main<-maintable[-which(is.na(maintable$genetic_sex)|maintable$sex!=maintable$genetic_sex),]
# main=delete_relative(main)
# main$age=year(main$image_date)-main$birth_year+(month(main$image_date)-main$birth_month)/12
# rm(maintable,kinship_3degree,Kinship_network)
# blk_info=read.csv("/../num_blk.csv",header=TRUE)
# basechar=main[,c(1,2,8,14,9:13)]
# set.seed(202204)
# G1=basechar[which(basechar$age<=65&basechar$sex==0),]
# G2=basechar[which(basechar$age>65&basechar$sex==0),]
# G3=basechar[which(basechar$age<=65&basechar$sex==1),]
# G4=basechar[which(basechar$age>65&basechar$sex==1),]
# d1=sample(G1$eid,size=ceiling(nrow(G1)/2),replace = FALSE)
# d2=sample(G2$eid,size=ceiling(nrow(G2)/2),replace = FALSE)
# d3=sample(G3$eid,size=ceiling(nrow(G3)/2),replace = FALSE)
# d4=sample(G4$eid,size=ceiling(nrow(G4)/2),replace = FALSE)
# v1=setdiff(G1$eid,d1)
# v2=setdiff(G2$eid,d2)
# v3=setdiff(G3$eid,d3)
# v4=setdiff(G4$eid,d4)
# discover=sort(c(d1,d2,d3,d4))
# verify=sort(c(v1,v2,v3,v4))
# base.dis=basechar[which(basechar$eid%in%discover),] ###baseline characteristics of the discovery set
# base.ver=basechar[which(basechar$eid%in%verify),] ####baseline charac. of the verification set
################ the genotype #################

##For illustration purppose we used chr1.bed, chr1.bim, chr1.fam as the file names
      bedfile="chr1.bed"
      geno.use<<-list()
      map.use <<- list()
      bimfile="chr1.bim"
      famfile="chr1.fam"
      data.chr<-read.plink(bedfile,bimfile,famfile,na.strings="0")
      geno0=data.chr$genotypes
      geno=geno0[which(rownames(geno0)%in%basechar$eid),] ## extract those overlapping with basechar
      eid_keep=rownames(geno) ###left 40123 subjects
      
      ###baseline characteristics matrix, with columns: eid, sex, ethnicity, age, pc1-5.
      ########Quality Control below 
      
      snpsum=col.summary(geno)
      ##remove snps with call rate less than 90%
      call.cut<-0.9
      call.keep <- with(snpsum, (!is.na(Call.rate) & Call.rate>= call.cut))
      call.keep[is.na(call.keep)] <- FALSE
      geno <- geno[, call.keep]
      map.call.keep <- data.chr$map[call.keep, ]
      snpsum.call.keep <- snpsum[call.keep, ]
      # remove SNPs with HWE test p-value less than 10-7
      HWE.cut.off <- 10^-7
      snpsum.control <- col.summary(geno)
      hwe.keep <- with(snpsum.control, (!is.na(z.HWE) & (abs(z.HWE)<abs(qnorm(HWE.cut.off/2)))))
      hwe.keep[is.na(hwe.keep)] <- FALSE
      geno.use <- geno[, hwe.keep]
      map.use <<- map.call.keep[hwe.keep, ]
      geno=geno.use[which(rownames(geno.use)%in%basechar$eid),]
      save(basechar, geno, geno.use,map.use,file = paste0(wd, "geno_transed.RData"))

###########initialize the responses#############

n.pheno=nrow(basechr)

pheno=array(0,dim=c(55,55,n.pheno))
#field="25753_2_0"
for(m in 1:n.dis){
  id=base.dis$eid[m]
  field=bulk$V2[which(bulk$V1==id)]
  field=field[1]
  ##########
  filename=paste0("/../basket-50595/field25753/",id,"_",field,".txt") 
  ######## to be modified to the specific folder of pheno data
  partial_cor=fread(filename)
  B <- diag(0,55)
  B[upper.tri(B,diag = F)]<- as.numeric(partial_cor)
  B <- B+t(B)
  pheno[,,m]<-B
}
save(pheno,file = paste0(wd, "pheno.RData"))

full_dist_mat=BallDistanceVector(pheno)
save(full_dist_mat,file = paste0(wd, "full_dist_mat.RData"))