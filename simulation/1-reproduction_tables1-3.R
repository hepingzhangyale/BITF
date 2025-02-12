# Simulation of Ball Impurity 
# This file provides the data generation process and the random seed to reproduce
# Tables 1-3 corresponding to Models 2.1a, 2.1b, 2.2, 2.3, 2.4, 3.1, 3.2, and 3.3
# To reproduce the corresponding simulation results, please specify the model number
# in the form of a character array, e.g. "2.1" as the input of the wrapper. 
# For Models 3.1-3.3, an outlier ratio (percentage) is also needed as "ops",
# the default is set to 1.

rm(list=ls())
parent.path<- ### specify your parent path
source(paste0(parent.path,"simulation/simu_rank_matrix.R"))
reproduction<-function(Model="2.1a",ops=1){
  if(Model=="2.1a"){
    library(Rcpp)
    library(energy) #for Distance Correlation
    library(Ball)
    library(VGAM)
    sourceCpp(paste0(parent.path,"code/BallDistanceVector.cpp"))
    sourceCpp(paste0(parent.path,"code/BallImpurity.cpp"))
    # calculate the frequency of the true signals being selected 
    Sel.rate<-function(n,d=10,true.v,rank.mtx) {
      # Input
      # n        :  the sample size
      # d        :  coeficient of cutoffs
      # true.v   :  the true variables index
      # rank.mtx :  the ranked index matrix by screening method for the 1000 replications
      #             each column corresponds the ranked index in one replication.
      # Output
      # rate     :  the proportions that every single active predictor is selected 
      #             for a given model size, which is defauted c[n/log(n)], in the 1,000 replications.
      rank.mtx.sel<-rank.mtx[1:d,]
      r<-min(dim(rank.mtx)[2],length(rank.mtx)) #repeat times
      p0<-nrow(true.v)
      R<-matrix(0,p0,r)
      rate<-rep(NA,p0)
      for (i in 1:p0) {
        for (j in 1:r) {R[i,j]<-(min(abs(rank.mtx.sel[,j]-true.v[i,j]))==0) }
        rate[i]<-mean(R[i,])
      }
      return(rate)
    }
    # calculate the frequency all true signals being selected.
    Sel.rate.all<-function(n,d=10,true.v,rank.mtx) {
      # Input
      # n        :  the sample size
      # d        :  coeficient of cutoffs
      # true.v   :  the true variables index
      # rank.mtx :  the ranked index matrix by screening method for the 1000 replications
      #             each column corresponds the ranked index in one replication.
      # Output
      # rate     :  the proportions that every single active predictor is selected 
      #             for a given model size, which is defauted c[n/log(n)], in the 1,000 replications.
      rank.mtx.sel<-rank.mtx[1:d,]
      r<-min(dim(rank.mtx)[2],length(rank.mtx))
      p0<-nrow(true.v)
      R<-matrix(0,p0,r)
      rate<-rep(NA,p0)
      for (i in 1:p0) {
        for (j in 1:r) {R[i,j]<-(min(abs(rank.mtx.sel[,j]-true.v[i,j]))==0) }
      }
      rate<-mean(apply(R,2,min))
      return(rate)
    }
    n=200
    p=100
    e_X=0.05 #noise level on X
    e_Y=0.05 #noise level on Y
    b_X=4  ##coefficient generating Z
    b_Y=3  ## coefficient generating Y
    repeat_times=100 #repeat times
    
    IG_result_weighted_sum=matrix(NA,nrow=p,ncol=repeat_times)
    BCor_result=matrix(NA,nrow=p,ncol=repeat_times)
    Dcor_result=matrix(NA,nrow=p,ncol=repeat_times)
    true_set_rec=matrix(NA,nrow=4,ncol=repeat_times)
    
    for (r in 1:repeat_times)
    {
      #Generate X
      set.seed(20250124+r)
      X=matrix(rbinom(n*p,1,0.5),nrow=n, ncol=p)
      true_set=sample(1:100,size=4)
      #add noise on X
      Z=apply(X,c(1,2),function(x) b_X*x+rnorm(1,mean=0,sd=e_X))
      Z=scale(Z)
      #Generate Y
      Y=list()
      for(i in 1:n)
      {
        Y1=b_Y*Z[i,true_set[1]]*Z[i,true_set[2]]*Z[i,true_set[3]]+b_Y*Z[i,true_set[4]]+rnorm(1,mean=0,sd=e_Y)
        
        Y[[i]]=c(Y1)
      }
      distance_matrix=BallDistanceVector(Y)
      
      Y_matrix=matrix(unlist(Y),nrow=n, ncol=length(Y[[1]]), byrow=T)
      Y_matrix=scale(Y_matrix)
      for(i in 1:n)
      {
        Y[[i]]=Y_matrix[i,]
      }
      
      d=quantile(distance_matrix,0.5,na.rm=T)
      BI0=BallImpurity(n,distance_matrix,1:n,d)
      #  IG_max_abs=rep(NA,p) #Information Gain
      IG_weighted_sum=rep(NA,p) #Information Gain
      BCor=rep(NA,p) #Ball Correlation    
      Dcor=rep(NA,p) #Distance Correlation
      for(j in 1:p)
      {
        tmp1=which(X[,j]==1)
        tmp0=which(X[,j]==0)
        #Information Gain on Ball Impurity: BI_0(we take 1 here for convenience)-1/n(BI_x1*n1+BI_x2)
        IG_weighted_sum[j]=BI0-(BallImpurity(n,distance_matrix,tmp1,d)*length(tmp1)/n+BallImpurity(n,distance_matrix,tmp0,d)*length(tmp0)/n)
        
        BCor[j]= bcor(X[,j],Y_matrix)
        Dcor[j]= dcor(X[,j],Y_matrix)
      }
      IG_result_weighted_sum[,r]=order(IG_weighted_sum,decreasing=T)
      BCor_result[,r]=order(BCor,decreasing=T)
      Dcor_result[,r]=order(Dcor,decreasing=T)
      true_set_rec[,r]=true_set
    }
    result=list()
    select_number=10 #2*floor(n/log(n))
    result<-c(result,
              list("Indicator P_m of BI"=c(Sel.rate(n,select_number,true_set_rec,IG_result_weighted_sum),"Indicator P_a of BI"=Sel.rate.all(n,select_number,true_set_rec,IG_result_weighted_sum))),
              list("Indicator P_m of BCor"=c(Sel.rate(n,select_number,true_set_rec,BCor_result),"Indicator P_a of BCor"=Sel.rate.all(n,select_number,true_set_rec,BCor_result))),
              list("Indicator P_m of Dcor"=c(Sel.rate(n,select_number,true_set_rec,Dcor_result),"Indicator P_a of Dcor"=Sel.rate.all(n,select_number,true_set_rec,Dcor_result))))
    save.image(file="Model2.1a.Rdata")
    #IG describes the difference of 3 distributions!!!
  }else 
    if(Model=="2.1b"){
    library(Rcpp)
    library(energy) #for Distance Correlation
    library(Ball)
    library(VGAM)
    sourceCpp(paste0(parent.path,"code/BallDistanceVector.cpp"))
    sourceCpp(paste0(parent.path,"code/BallImpurity.cpp"))
    # calculate the frequency of the true signals being selected 
    Sel.rate<-function(n,d=10,true.v,rank.mtx) {
    # Input
    # n        :  the sample size
    # d        :  coeficient of cutoffs
    # true.v   :  the true variables index
    # rank.mtx :  the ranked index matrix by screening method for the 1000 replications
    #             each column corresponds the ranked index in one replication.
    # Output
    # rate     :  the proportions that every single active predictor is selected 
    #             for a given model size, which is defauted c[n/log(n)], in the 1,000 replications.
    rank.mtx.sel<-rank.mtx[1:d,]
    r<-min(dim(rank.mtx)[2],length(rank.mtx)) #repeat times
    p0<-nrow(true.v)
    R<-matrix(0,p0,r)
    rate<-rep(NA,p0)
    for (i in 1:p0) {
      for (j in 1:r) {R[i,j]<-(min(abs(rank.mtx.sel[,j]-true.v[i,j]))==0) }
      rate[i]<-mean(R[i,])
    }
    return(rate)
  }
  # calculate the frequency all true signals being selected.
    Sel.rate.all<-function(n,d=10,true.v,rank.mtx) {
    # Input
    # n        :  the sample size
    # d        :  coeficient of cutoffs
    # true.v   :  the true variables index
    # rank.mtx :  the ranked index matrix by screening method for the 1000 replications
    #             each column corresponds the ranked index in one replication.
    # Output
    # rate     :  the proportions that every single active predictor is selected 
    #             for a given model size, which is defauted c[n/log(n)], in the 1,000 replications.
    rank.mtx.sel<-rank.mtx[1:d,]
    r<-min(dim(rank.mtx)[2],length(rank.mtx))
    p0<-nrow(true.v)
    R<-matrix(0,p0,r)
    rate<-rep(NA,p0)
    for (i in 1:p0) {
      for (j in 1:r) {R[i,j]<-(min(abs(rank.mtx.sel[,j]-true.v[i,j]))==0) }
    }
    rate<-mean(apply(R,2,min))
    return(rate)
  }
    n=200
    p=1000
    e_X=0.05 #noise level on X
    e_Y=0.05 #noise level on Y
    b_X=4  ##coefficient generating Z
    b_Y=3  ## coefficient generating Y
    repeat_times=100 #repeat times
    
    IG_result_weighted_sum=matrix(NA,nrow=p,ncol=repeat_times)
    BCor_result=matrix(NA,nrow=p,ncol=repeat_times)
    Dcor_result=matrix(NA,nrow=p,ncol=repeat_times)
    true_set_rec=matrix(NA,nrow=4,ncol=repeat_times)
    
    for (r in 1:repeat_times)
    {
      #Generate X
      set.seed(20250124+r)
      X=matrix(rbinom(n*p,1,0.5),nrow=n, ncol=p)
      true_set=sample(1:100,size=4)
      #add noise on X
      Z=apply(X,c(1,2),function(x) b_X*x+rnorm(1,mean=0,sd=e_X))
      Z=scale(Z)
      #Generate Y
      Y=list()
      for(i in 1:n)
      {
        Y1=b_Y*Z[i,true_set[1]]*Z[i,true_set[2]]*Z[i,true_set[3]]+b_Y*Z[i,true_set[4]]+rpareto(1,scale=0.05,shape=2)
        
        Y[[i]]=c(Y1)
      }
      distance_matrix=BallDistanceVector(Y)
      
      Y_matrix=matrix(unlist(Y),nrow=n, ncol=length(Y[[1]]), byrow=T)
      Y_matrix=scale(Y_matrix)
      for(i in 1:n)
      {
        Y[[i]]=Y_matrix[i,]
      }
      
      d=quantile(distance_matrix,0.5,na.rm=T)
      BI0=BallImpurity(n,distance_matrix,1:n,d)
      IG_weighted_sum=rep(NA,p) #Information Gain
      BCor=rep(NA,p) #Ball Correlation    
      Dcor=rep(NA,p) #Distance Correlation
      for(j in 1:p)
      {
        tmp1=which(X[,j]==1)
        tmp0=which(X[,j]==0)
        #Information Gain on Ball Impurity: BI_0(we take 1 here for convenience)-1/n(BI_x1*n1+BI_x2)
        IG_weighted_sum[j]=BI0-(BallImpurity(n,distance_matrix,tmp1,d)*length(tmp1)/n+BallImpurity(n,distance_matrix,tmp0,d)*length(tmp0)/n)
        
        BCor[j]= bcor(X[,j],Y_matrix)
        Dcor[j]= dcor(X[,j],Y_matrix)
      }
      IG_result_weighted_sum[,r]=order(IG_weighted_sum,decreasing=T)
      BCor_result[,r]=order(BCor,decreasing=T)
      Dcor_result[,r]=order(Dcor,decreasing=T)
      true_set_rec[,r]=true_set
    }
    
    result=list()
    
    select_number=10 #2*floor(n/log(n))
    result<-c(result,
              list("Indicator P_m of BI"=c(Sel.rate(n,select_number,true_set_rec,IG_result_weighted_sum),"Indicator P_a of BI"=Sel.rate.all(n,select_number,true_set_rec,IG_result_weighted_sum))),
              list("Indicator P_m of BCor"=c(Sel.rate(n,select_number,true_set_rec,BCor_result),"Indicator P_a of BCor"=Sel.rate.all(n,select_number,true_set_rec,BCor_result))),
              list("Indicator P_m of Dcor"=c(Sel.rate(n,select_number,true_set_rec,Dcor_result),"Indicator P_a of Dcor"=Sel.rate.all(n,select_number,true_set_rec,Dcor_result))))
    result
    save.image(file="Model2.1b.Rdata")
  }else 
    if(Model=="2.2"){
    library(Rcpp)
    library(energy) #for Distance Correlation
    library(Ball)
    library(VGAM)
    sourceCpp(paste0(parent.path,"code/BallDistanceVector.cpp"))
    sourceCpp(paste0(parent.path,"code/BallImpurity.cpp"))
    # calculate the frequency of the true signals being selected 
    Sel.rate<-function(n,d=10,true.v,rank.mtx) {
      # Input
      # n        :  the sample size
      # d        :  coeficient of cutoffs
      # true.v   :  the true variables index
      # rank.mtx :  the ranked index matrix by screening method for the 1000 replications
      #             each column corresponds the ranked index in one replication.
      # Output
      # rate     :  the proportions that every single active predictor is selected 
      #             for a given model size, which is defauted c[n/log(n)], in the 1,000 replications.
      rank.mtx.sel<-rank.mtx[1:d,]
      r<-min(dim(rank.mtx)[2],length(rank.mtx)) #repeat times
      p0<-nrow(true.v)
      R<-matrix(0,p0,r)
      rate<-rep(NA,p0)
      for (i in 1:p0) {
        for (j in 1:r) {R[i,j]<-(min(abs(rank.mtx.sel[,j]-true.v[i,j]))==0) }
        rate[i]<-mean(R[i,])
      }
      return(rate)
    }
    # calculate the frequency all true signals being selected.
    Sel.rate.all<-function(n,d=10,true.v,rank.mtx) {
      # Input
      # n        :  the sample size
      # d        :  coeficient of cutoffs
      # true.v   :  the true variables index
      # rank.mtx :  the ranked index matrix by screening method for the 1000 replications
      #             each column corresponds the ranked index in one replication.
      # Output
      # rate     :  the proportions that every single active predictor is selected 
      #             for a given model size, which is defauted c[n/log(n)], in the 1,000 replications.
      rank.mtx.sel<-rank.mtx[1:d,]
      r<-min(dim(rank.mtx)[2],length(rank.mtx))
      p0<-nrow(true.v)
      R<-matrix(0,p0,r)
      rate<-rep(NA,p0)
      for (i in 1:p0) {
        for (j in 1:r) {R[i,j]<-(min(abs(rank.mtx.sel[,j]-true.v[i,j]))==0) }
      }
      rate<-mean(apply(R,2,min))
      return(rate)
    }
    
    n=200
    p=1000
    e_X=0.1 #noise level on X
    e_Y=0.1 #noise level on Y
    b_X=4  ##coefficient generating Z
    noise_level2=0.1
    #b_Y=3  ## coefficient generating Y
    repeat_times=100 #repeat times
    
    IG_result_weighted_sum=matrix(NA,nrow=p,ncol=repeat_times)
    BCor_result=matrix(NA,nrow=p,ncol=repeat_times)
    Dcor_result=matrix(NA,nrow=p,ncol=repeat_times)
    true_set_rec=matrix(NA,nrow=4,ncol=repeat_times)
    
    for (r in 1:repeat_times)
    {
      #Generate X
      set.seed(20250124+r)
      X=matrix(rbinom(n*p,1,0.5),nrow=n, ncol=p)
      true_set=sample(1:1000,size=4)
      #add noise on X
      Z=apply(X,c(1,2),function(x) b_X*x+rnorm(1,mean=0,sd=e_X))
      Z=scale(Z)
      #Generate Y
      Y=list()
      for(i in 1:n)
      {
        Y1=(4-Z[i,true_set[1]]^2)*Z[i,true_set[3]]*Z[i,true_set[4]]+Z[i,true_set[2]]*Z[i,true_set[4]]^2+rnorm(1,mean=0,sd=noise_level2)
        Y2=(4-Z[i,true_set[2]]^2)*Z[i,true_set[3]]*Z[i,true_set[4]]+Z[i,true_set[1]]*Z[i,true_set[4]]^2+rnorm(1,mean=0,sd=noise_level2)
        Y[[i]]=c(Y1,Y2)
      }
      distance_matrix=BallDistanceVector(Y)
      
      Y_matrix=matrix(unlist(Y),nrow=n, ncol=length(Y[[1]]), byrow=T)
      Y_matrix=scale(Y_matrix)
      for(i in 1:n)
      {
        Y[[i]]=Y_matrix[i,]
      }
      
      d=quantile(distance_matrix,0.5,na.rm=T)/6
      BI0=BallImpurity(n,distance_matrix,1:n,d)
      #  IG_max_abs=rep(NA,p) #Information Gain
      IG_weighted_sum=rep(NA,p) #Information Gain
      BCor=rep(NA,p) #Ball Correlation    
      Dcor=rep(NA,p) #Distance Correlation
      for(j in 1:p)
      {
        tmp1=which(X[,j]==1)
        tmp0=which(X[,j]==0)
        #Information Gain on Ball Impurity: BI_0(we take 1 here for convenience)-1/n(BI_x1*n1+BI_x2)
        #IG_weighted_sum[j]=BI0-BImp(n,distance_matrix,1:n,tmp1,d)*length(tmp1)/n+BImp(n,distance_matrix,1:n,tmp0,d)*length(tmp0)/n
        IG_weighted_sum[j]=BI0-(BallImpurity(n,distance_matrix,tmp1,d)*length(tmp1)/n+BallImpurity(n,distance_matrix,tmp0,d)*length(tmp0)/n)
        
        BCor[j]= bcor(X[,j],Y_matrix)
        Dcor[j]= dcor(X[,j],Y_matrix)
      }
      IG_result_weighted_sum[,r]=order(IG_weighted_sum,decreasing=T)
      BCor_result[,r]=order(BCor,decreasing=T)
      Dcor_result[,r]=order(Dcor,decreasing=T)
      true_set_rec[,r]=true_set
    }
    
    result=list()
    
    select_number=10 #2*floor(n/log(n))
    result<-c(result,
              list("Indicator P_m of BI"=c(Sel.rate(n,select_number,true_set_rec,IG_result_weighted_sum),"Indicator P_a of BI"=Sel.rate.all(n,select_number,true_set_rec,IG_result_weighted_sum))),
              list("Indicator P_m of BCor"=c(Sel.rate(n,select_number,true_set_rec,BCor_result),"Indicator P_a of BCor"=Sel.rate.all(n,select_number,true_set_rec,BCor_result))),
              list("Indicator P_m of Dcor"=c(Sel.rate(n,select_number,true_set_rec,Dcor_result),"Indicator P_a of Dcor"=Sel.rate.all(n,select_number,true_set_rec,Dcor_result))))
    result
    save.image(file="Model2.2.Rdata")
  }else 
    if(Model=="2.3"){
    library(Ball)
    library(Rcpp)
    sourceCpp(paste0(parent.path,"code/BallImpurity.cpp"))
    sourceCpp(paste0(parent.path,"code/BallDistanceVector.cpp"))
    source(paste0(parent.path,"code/biRegressionTree.R"))
    rep=100
    n=200          #### sample size
    signal_level1 <- 1
    signal_level2 <- 0.5
    noise_level1 <- 0.2
    noise_level2 <- 0.25
    beta0=c(0.5,1)*signal_level2
    results=NULL#
    rank_mat=NULL
    for(i in 1:rep){
      seed=202279+i
      set.seed(seed)
      X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
      #### the original predictors are binary with p=0.5
      sig = sample(1:p_X,size=n_sig,replace=F)
      #### randomly select n_sig of them to be active predictors
      Z=X[,sig]*signal_level1+matrix(rnorm(n*n_sig,mean=0,sd=noise_level1),nrow=n)
      ### transform X into Z by magnifying by signal_level1 and adding perturbation 
      ## with sd=noise_level1
      transformed_cov=cbind(Z[,1]*Z[,2],Z[,3])
      
      Y=transformed_cov%*%beta0
      rank_matrix_BCor=simu_rank_matrix_euclid(Y,X,sig,method="BCor")$rank
      bcorsig=sig[which(rank_matrix_BCor<=10)]
      biTree<- buildTree(list(X=X,Y=Y), leafFunc=regLeaf, method='BallImpurity', ops=c(0, minLeafSize=20), maxDepth=5, lossChoice="loss2")
      newTree=PruneTreeComplete(biTree,0.10,depth=8)
      bisig=unique(getnodes(newTree))
      if(length(bisig)>10){bisig=bisig[1:10]}
      rank_matrix_BT=c(as.integer(sig[1]%in%btsig),as.integer(sig[2]%in%btsig),as.integer(sig[3]%in%btsig))
      rank_mat=rbind(rank_mat,c(rank_matrix_BT))
    }
    ratesi=apply((rank_mat<=10)&(rank_mat>0),2,sum)
    
    ratesallbt=sum(apply((rank_mat[,(1:n_sig)]>0)&(rank_mat[,(1:n_sig)]<=10),1,sum)>=3)
    results=rbind(results,c(ratesi,ratesallbt))
    save.image("Model2.3.RData")
  }else 
    if(Model=="2.4"){
    library(Ball)
    library(Rcpp)
    sourceCpp(paste0(parent.path,"code/BallImpurity.cpp"))
    sourceCpp(paste0(parent.path,"code/BallDistanceVector.cpp"))
    source(paste0(parent.path,"code/biRegressionTree.R"))
    n=200          #### sample size
    p_X <- 3    #### number of columns in X
    n_sig=3
    signal_level1 <- 1
    signal_level2 <- 3
    e_Y=2
    beta0=c(1,1)*signal_level2
    results=NULL#
    rank_mat=NULL
    rep=1
    for(i in 1:rep){
      seed=202279+i
      set.seed(seed)
      X=mvrnorm(n=200,mu=c(5,10,15),Sigma = matrix(c(3,1,-1,1,4,-0.5,1,-0.5,2),ncol=3))
      Z=cbind(X[,1]>=5,X[,2]>=11,X[,3]>=14)
      transformed_cov=cbind(Z[,1]*Z[,2],Z[,3])
      Y=transformed_cov%*%beta0+rnorm(1,mean=0,sd=e_Y)
      Y=scale(Y)
      D=BallDistanceVector(Y)
      biTree<- buildTree(list(X=X,Y=Y), leafFunc=regLeaf, method='BallImpurity',D, ops=c(0, minLeafSize=10), maxDepth=4, lossChoice="loss2")
    }
    biTree<-PruneTreeComplete(biTree,1.5,4)
    save.image("model2.4.RData")
  }else 
    if(Model=="3.1"){
    library(snow)
    library(snowfall)
    library(tictoc)
    library(movMF)
    library(MASS)
    sim_BI4=function(r){
      seed=startingseed+r
      set.seed(seed)
      X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
      sig = sample(1:p_X,size=n_sig,replace=F)
      Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=znoise_level),nrow=n)
      Y=matrix(NA,nrow=n,ncol=p_Y)
      out_id=sample(1:n,size=floor(n*outlier_ratio))
      disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
      if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
        disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
      }else{
        disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
      disturb=disturb_sample[disturb_filtered_id,]
      for(i in 1:dim(X)[1])
      {
        mu=Z[i,,drop=FALSE]%*%B_True
        if(i %in% out_id){
          mu=mu+disturb[match(i,out_id),]
        }
        Y[i,]=rmovMF(1,theta=mu)
      }
      #### Add noise
      Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(ynoise_level^2,5)))
      res_r=simu_rank_matrix(Y,X,sig,method="BI",dweight=4)$rank
      return(res_r)
    }
    sim_BCor=function(r){
      seed=startingseed+r
      set.seed(seed)
      X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
      sig = sample(1:p_X,size=n_sig,replace=F)
      Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=znoise_level),nrow=n)
      Y=matrix(NA,nrow=n,ncol=p_Y)
      out_id=sample(1:n,size=floor(n*outlier_ratio))
      disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
      if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
        disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
      }else{
        disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
      disturb=disturb_sample[disturb_filtered_id,]
      for(i in 1:dim(X)[1])
      {
        mu=Z[i,,drop=FALSE]%*%B_True
        if(i %in% out_id){
          mu=mu+disturb[match(i,out_id),]
        }
        Y[i,]=rmovMF(1,theta=mu)
      }
      #### Add noise
      Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(ynoise_level^2,5)))
      res_r=simu_rank_matrix(Y,X,sig,method="BCor")$rank
      return(res_r)
    }
    sim_DCor=function(r){
      seed=startingseed+r
      set.seed(seed)
      X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
      sig = sample(1:p_X,size=n_sig,replace=F)
      #print(sig)
      Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=znoise_level),nrow=n)
      Y=matrix(NA,nrow=n,ncol=p_Y)
      out_id=sample(1:n,size=floor(n*outlier_ratio))
      disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
      if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
        disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
      }else{
        disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
      disturb=disturb_sample[disturb_filtered_id,]
      for(i in 1:dim(X)[1])
      {
        mu=Z[i,,drop=FALSE]%*%B_True
        if(i %in% out_id){
          mu=mu+disturb[match(i,out_id),]
        }
        Y[i,]=rmovMF(1,theta=mu)
      }
      #### Add noise
      Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(ynoise_level^2,5)))
      res_r=simu_rank_matrix(Y,X,sig,method="DCor")$rank
      return(res_r)
    }
    
    sim_maxBI4=function(r){
      seed=startingseed+r
      set.seed(seed)
      X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
      sig = sample(1:p_X,size=n_sig,replace=F)
      Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=znoise_level),nrow=n)
      Y=matrix(NA,nrow=n,ncol=p_Y)
      out_id=sample(1:n,size=floor(n*outlier_ratio))
      disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
      if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
        disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
      }else{
        disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
      disturb=disturb_sample[disturb_filtered_id,]
      for(i in 1:dim(X)[1])
      {
        mu=Z[i,,drop=FALSE]%*%B_True
        if(i %in% out_id){
          mu=mu+disturb[match(i,out_id),]
        }
        Y[i,]=rmovMF(1,theta=mu)
      }
      #### Add noise
      Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(ynoise_level^2,5)))
      res_r=simu_rank_matrix_euclid(Y,X,sig,method="BI",dweight=4,dist_method = "maximum")$rank
      return(res_r)
    }
    sim_maxBCor=function(r){
      seed=startingseed+r
      set.seed(seed)
      X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
      sig = sample(1:p_X,size=n_sig,replace=F)
      Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=znoise_level),nrow=n)
      Y=matrix(NA,nrow=n,ncol=p_Y)
      out_id=sample(1:n,size=floor(n*outlier_ratio))
      disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
      if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
        disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
      }else{
        disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
      disturb=disturb_sample[disturb_filtered_id,]
      for(i in 1:dim(X)[1])
      {
        mu=Z[i,,drop=FALSE]%*%B_True
        if(i %in% out_id){
          mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
        }
        Y[i,]=rmovMF(1,theta=mu)
      }
      #### Add noise
      Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(ynoise_level^2,5)))
      res_r=simu_rank_matrix_euclid(Y,X,sig,method="BCor",dist_method = "maximum")$rank
      return(res_r)
    }
    sim_maxDCor=function(r){
      seed=startingseed+r
      set.seed(seed)
      X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
      sig = sample(1:p_X,size=n_sig,replace=F)
      Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=znoise_level),nrow=n)
      Y=matrix(NA,nrow=n,ncol=p_Y)
      out_id=sample(1:n,size=floor(n*outlier_ratio))
      disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
      if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
        disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
      }else{
        disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
      disturb=disturb_sample[disturb_filtered_id,]
      for(i in 1:dim(X)[1])
      {
        mu=Z[i,,drop=FALSE]%*%B_True
        if(i %in% out_id){
          mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
        }
        Y[i,]=rmovMF(1,theta=mu)
      }
      #### Add noise
      Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(ynoise_level^2,5)))
      res_r=simu_rank_matrix_euclid(Y,X,sig,method="DCor",dist_method = "maximum")$rank
      return(res_r)
    }
    rep=100
    n=200          #### sample size
    p_X <- 500    #### number of columns in X
    n_sig=3
    p_Y<-5
    
    startingseed=202411
    znoise_level=0.2
    ynoise_level=0.3
    B_True=matrix(c(0,1,1,1,0,0,2,2,0,0,0,3,0,0,0),nrow=3,byrow=T)
    mu0=c(1,0,0,0,0)#
    outlier_ratio=0.01 #to be varied
    rank_mat=list(nh_BI4=matrix(rep(0,rep*n_sig),nrow=rep),nh_BCor=matrix(rep(0,rep*n_sig),nrow=rep),nh_DCor=matrix(rep(0,rep*n_sig),nrow=rep),
                  max_BI4=matrix(rep(0,rep*n_sig),nrow=rep),max_BCor=matrix(rep(0,rep*n_sig),nrow=rep),max_DCor=matrix(rep(0,rep*n_sig),nrow=rep))
    
    tic()
    sfInit(cpus=3,parallel=TRUE) 
    sfExportAll()
    sfLibrary(MASS)
    sfLibrary(movMF)
    sfLibrary("Ball")
    rank_mat$nh_BI4=sfSapply(seq(1:100),sim_BI4)
    sfStop()
    toc()
    
    tic()
    sfInit(cpus=3,parallel=TRUE) 
    sfExportAll()
    sfLibrary(MASS)
    sfLibrary(movMF)
    sfLibrary(Ball)
    rank_mat$nh_BCor=sfSapply(seq(1:100),sim_BCor)
    sfStop()
    toc()
    
    tic()
    sfInit(cpus=3,parallel=TRUE) 
    sfExportAll()
    sfLibrary(MASS)
    sfLibrary(movMF)
    sfLibrary(Ball)
    rank_mat$nh_DCor=sfSapply(seq(1:100),sim_DCor)
    sfStop()
    toc()
    
    
    
    tic()
    sfInit(cpus=3,parallel=TRUE) 
    sfExportAll()
    sfLibrary(MASS)
    sfLibrary(movMF)
    sfLibrary(Ball)
    rank_mat$max_BI4=sfSapply(seq(1:100),sim_maxBI4)
    sfStop()
    toc()
    
    tic()
    sfInit(cpus=3,parallel=TRUE) 
    sfExportAll()
    sfLibrary(MASS)
    sfLibrary(movMF)
    rank_mat$max_BI8=sfSapply(seq(1:100),sim_maxBI8)
    sfStop()
    toc()
    
    tic()
    sfInit(cpus=3,parallel=TRUE) 
    sfExportAll()
    sfLibrary(MASS)
    sfLibrary(movMF)
    rank_mat$max_BCor=sfSapply(seq(1:100),sim_maxBCor)
    sfStop()
    toc()
    
    tic()
    sfInit(cpus=3,parallel=TRUE) 
    sfExportAll()
    sfLibrary(MASS)
    sfLibrary(movMF)
    rank_mat$max_DCor=sfSapply(seq(1:100),sim_maxDCor)
    sfStop()
    toc()
    
    thres=21
    rowSums(rank_mat$nh_BI4<thres)
    sum(colSums(rank_mat$nh_BI4<thres)>2)
    rowSums(rank_mat$nh_BCor<thres)
    sum(colSums(rank_mat$nh_BCor<thres)>2)
    rowSums(rank_mat$nh_DCor<thres)
    sum(colSums(rank_mat$nh_DCor<thres)>2)
    
    rowSums(rank_mat$max_BI4<thres)
    sum(colSums(rank_mat$max_BI4<thres)>2)
    rowSums(rank_mat$max_BCor<thres)
    sum(colSums(rank_mat$max_BCor<thres)>2)
    rowSums(rank_mat$max_DCor<thres)
    sum(colSums(rank_mat$max_DCor<thres)>2)
    
    save.image(paste0("md31_outlier_",ops,".RData"))
  }else 
    if(Model=="3.2"){
    sim_BI4=function(r){
      seed=startingseed+r
      set.seed(seed)
      X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
      sig = sample(1:p_X,size=n_sig,replace=F)
      #print(sig)
      Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
      mu0=c(1,0,0,0,0)#*signal_level
      Y=matrix(NA,nrow=n,ncol=p_Y)
      out_id=sample(1:n,size=floor(n*outlier_ratio))
      disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
      if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
        disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
      }else{
        disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
      disturb=disturb_sample[disturb_filtered_id,]
      for(i in 1:dim(X)[1])
      {
        mu=Z[i,,drop=FALSE]%*%B_True
        if(i %in% out_id){
          mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
        }
        Y[i,]=rmovMF(1,theta=mu)
      }
      #### Add noise
      Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
      D=nhdist(Y)
      res_r=simu_rank_matrix(Y,X,sig,D,method="BI",dweight=4)$rank
      return(res_r)
    }
    sim_BCor=function(r){
      seed=startingseed+r
      set.seed(seed)
      X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
      sig = sample(1:p_X,size=n_sig,replace=F)
      #print(sig)
      Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
      mu0=c(1,0,0,0,0)#*signal_level
      Y=matrix(NA,nrow=n,ncol=p_Y)
      out_id=sample(1:n,size=floor(n*outlier_ratio))
      disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
      if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
        disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
      }else{
        disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
      disturb=disturb_sample[disturb_filtered_id,]
      for(i in 1:dim(X)[1])
      {
        mu=Z[i,,drop=FALSE]%*%B_True
        if(i %in% out_id){
          mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
        }
        Y[i,]=rmovMF(1,theta=mu)
      }
      #### Add noise
      Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
      D=nhdist(Y)
      res_r=simu_rank_matrix(Y,X,sig,D,method="BCor")$rank
      return(res_r)
    }
    sim_DCor=function(r){
      seed=startingseed+r
      set.seed(seed)
      X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
      sig = sample(1:p_X,size=n_sig,replace=F)
      #print(sig)
      Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
      mu0=c(1,0,0,0,0)#*signal_level
      Y=matrix(NA,nrow=n,ncol=p_Y)
      out_id=sample(1:n,size=floor(n*outlier_ratio))
      disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
      if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
        disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
      }else{
        disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
      disturb=disturb_sample[disturb_filtered_id,]
      for(i in 1:dim(X)[1])
      {
        mu=Z[i,,drop=FALSE]%*%B_True
        if(i %in% out_id){
          mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
        }
        Y[i,]=rmovMF(1,theta=mu)
      }
      #### Add noise
      Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
      D=nhdist(Y)
      res_r=simu_rank_matrix(Y,X,sig,D,method="DCor")$rank
      return(res_r)
    }
    sim_euBI4=function(r){
      seed=startingseed+r
      set.seed(seed)
      X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
      sig = sample(1:p_X,size=n_sig,replace=F)
      #print(sig)
      Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
      mu0=c(1,0,0,0,0)#*signal_level
      Y=matrix(NA,nrow=n,ncol=p_Y)
      out_id=sample(1:n,size=floor(n*outlier_ratio))
      disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
      if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
        disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
      }else{
        disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
      disturb=disturb_sample[disturb_filtered_id,]
      for(i in 1:dim(X)[1])
      {
        mu=Z[i,,drop=FALSE]%*%B_True
        if(i %in% out_id){
          mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
        }
        Y[i,]=rmovMF(1,theta=mu)
      }
      #### Add noise
      Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
      D=as.matrix(dist(Y,method = "euclidean"))
      res_r=simu_rank_matrix(Y,X,sig,D,method="BI",dweight=4)$rank
      return(res_r)
    }
    sim_euBCor=function(r){
      seed=startingseed+r
      set.seed(seed)
      X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
      sig = sample(1:p_X,size=n_sig,replace=F)
      #print(sig)
      Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
      mu0=c(1,0,0,0,0)#*signal_level
      Y=matrix(NA,nrow=n,ncol=p_Y)
      out_id=sample(1:n,size=floor(n*outlier_ratio))
      disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
      if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
        disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
      }else{
        disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
      disturb=disturb_sample[disturb_filtered_id,]
      for(i in 1:dim(X)[1])
      {
        mu=Z[i,,drop=FALSE]%*%B_True
        if(i %in% out_id){
          mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
        }
        Y[i,]=rmovMF(1,theta=mu)
      }
      #### Add noise
      Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
      D=as.matrix(dist(Y,method = "euclidean"))
      res_r=simu_rank_matrix(Y,X,sig,D,method="BCor")$rank
      return(res_r)
    }
    sim_euDCor=function(r){
      seed=startingseed+r
      set.seed(seed)
      X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
      sig = sample(1:p_X,size=n_sig,replace=F)
      #print(sig)
      Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
      mu0=c(1,0,0,0,0)#*signal_level
      Y=matrix(NA,nrow=n,ncol=p_Y)
      out_id=sample(1:n,size=floor(n*outlier_ratio))
      disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
      if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
        disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
      }else{
        disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
      disturb=disturb_sample[disturb_filtered_id,]
      for(i in 1:dim(X)[1])
      {
        mu=Z[i,,drop=FALSE]%*%B_True
        if(i %in% out_id){
          mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
        }
        Y[i,]=rmovMF(1,theta=mu)
      }
      #### Add noise
      Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
      D=as.matrix(dist(Y,method = "euclidean"))
      res_r=simu_rank_matrix(Y,X,sig,D,method="DCor")$rank
      return(res_r)
    }
    
    sim_maxBI4=function(r){
      seed=startingseed+r
      set.seed(seed)
      X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
      sig = sample(1:p_X,size=n_sig,replace=F)
      #print(sig)
      Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
      mu0=c(1,0,0,0,0)#*signal_level
      Y=matrix(NA,nrow=n,ncol=p_Y)
      out_id=sample(1:n,size=floor(n*outlier_ratio))
      disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
      if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
        disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
      }else{
        disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
      disturb=disturb_sample[disturb_filtered_id,]
      for(i in 1:dim(X)[1])
      {
        mu=Z[i,,drop=FALSE]%*%B_True
        if(i %in% out_id){
          mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
        }
        Y[i,]=rmovMF(1,theta=mu)
      }
      #### Add noise
      Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
      D=as.matrix(dist(Y,method = "maximum"))
      res_r=simu_rank_matrix(Y,X,sig,D,method="BI",dweight=4)$rank
      return(res_r)
    }
    sim_maxBCor=function(r){
      seed=startingseed+r
      set.seed(seed)
      X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
      sig = sample(1:p_X,size=n_sig,replace=F)
      #print(sig)
      Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
      mu0=c(1,0,0,0,0)#*signal_level
      Y=matrix(NA,nrow=n,ncol=p_Y)
      out_id=sample(1:n,size=floor(n*outlier_ratio))
      disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
      if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
        disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
      }else{
        disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
      disturb=disturb_sample[disturb_filtered_id,]
      for(i in 1:dim(X)[1])
      {
        mu=Z[i,,drop=FALSE]%*%B_True
        if(i %in% out_id){
          mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
        }
        Y[i,]=rmovMF(1,theta=mu)
      }
      #### Add noise
      Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
      D=as.matrix(dist(Y,method = "maximum"))
      res_r=simu_rank_matrix(Y,X,sig,D,method="BCor")$rank
      return(res_r)
    }
    sim_maxDCor=function(r){
      seed=startingseed+r
      set.seed(seed)
      X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
      sig = sample(1:p_X,size=n_sig,replace=F)
      #print(sig)
      Z=X[,sig]*1+matrix(rnorm(n*n_sig,mean=0,sd=0.2),nrow=n)
      mu0=c(1,0,0,0,0)#*signal_level
      Y=matrix(NA,nrow=n,ncol=p_Y)
      out_id=sample(1:n,size=floor(n*outlier_ratio))
      disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
      if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
        disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))
      }else{
        disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
      disturb=disturb_sample[disturb_filtered_id,]
      for(i in 1:dim(X)[1])
      {
        mu=Z[i,,drop=FALSE]%*%B_True
        if(i %in% out_id){
          mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
        }
        Y[i,]=rmovMF(1,theta=mu)
      }
      #### Add noise
      Y=Y+mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(0.3^2,5)))
      D=as.matrix(dist(Y,method = "maximum"))
      res_r=simu_rank_matrix(Y,X,sig,D,method="DCor")$rank
      return(res_r)
    }
    rep=100
    n=200          #### sample size
    p_X <- 500    #### number of columns in X
    n_sig=2
    p_Y<-5
    
    outlier_ratio=0.1
    startingseed=202479
    B_True=matrix(c(0,2/3,-1/3,1,0,0,-1/3,-1/3,4/3,0)*1.5,nrow=2,byrow=T)
    
    rank_mat=list(nh_BI4=matrix(rep(0,rep*n_sig),nrow=rep),nh_BCor=matrix(rep(0,rep*n_sig),nrow=rep),nh_DCor=matrix(rep(0,rep*n_sig),nrow=rep),
                  max_BI4=matrix(rep(0,rep*n_sig),nrow=rep),max_BCor=matrix(rep(0,rep*n_sig),nrow=rep),max_DCor=matrix(rep(0,rep*n_sig),nrow=rep))
    
    tic()
    sfInit(cpus=4,parallel=TRUE) 
    sfExportAll()
    sfLibrary(MASS)
    sfLibrary(movMF)
    sfLibrary(Ball)
    rank_mat$nh_BI4[1:30,]=sfSapply(seq(1:30),sim_BI4)
    sfStop()
    toc()
    sum(rank_mat$nh_BI4[1:100,1]<11)
    sum(rank_mat$nh_BI4[1:100,2]<11)
    sum(rowSums(rank_mat$nh_BI4[1:100,]<11)>1)
    
    tic()
    sfInit(cpus=4,parallel=TRUE) 
    sfExportAll()
    sfLibrary(MASS)
    sfLibrary(movMF)
    sfLibrary(Ball)
    rank_mat$nh_BCor[1:100,]=sfSapply(seq(1:100),sim_BCor)
    sfStop()
    toc()
    
    
    tic()
    sfInit(cpus=3,parallel=TRUE) 
    sfExportAll()
    sfLibrary(MASS)
    sfLibrary(movMF)
    sfLibrary(Ball)
    rank_mat$nh_DCor=sfSapply(seq(1:100),sim_DCor)
    sfStop()
    toc()
    
    
    
    
    tic()
    sfInit(cpus=4,parallel=TRUE) 
    sfExportAll()
    sfLibrary(MASS)
    sfLibrary(movMF)
    rank_mat$eu_BI4[1:100,]=sfSapply(seq(1:100),sim_euBI4)
    sfStop()
    toc()
    sum(rank_mat$eu_BI4[1:100,1]<11)
    sum(rank_mat$eu_BI4[1:100,2]<11)
    sum(rowSums(rank_mat$eu_BI4[1:100,]<11)>1)
    
    tic()
    sfInit(cpus=3,parallel=TRUE) 
    sfExportAll()
    sfLibrary(MASS)
    sfLibrary(movMF)
    rank_mat$eu_DCor[1:100,]=sfSapply(seq(1:100),sim_euDCor)
    sfStop()
    toc()
    sum(rank_mat$eu_DCor[1:100,1]<11)
    sum(rank_mat$eu_DCor[1:100,2]<11)
    sum(rowSums(rank_mat$eu_DCor[1:100,]<11)>1)
    
    tic()
    sfInit(cpus=3,parallel=TRUE) 
    sfExportAll()
    sfLibrary(MASS)
    sfLibrary(movMF)
    rank_mat$eu_BCor[1:100,]=sfSapply(seq(1:100),sim_euBCor)
    sfStop()
    toc()
    sum(rank_mat$eu_BCor[,1]<11)
    sum(rank_mat$eu_BCor[,2]<11)
    sum(rowSums(rank_mat$eu_BCor[],]<11)>1)
    
    
    tic()
    sfInit(cpus=4,parallel=TRUE) 
    sfExportAll()
    sfLibrary(MASS)
    sfLibrary(movMF)
    #rank_mat$nh_BCor=sfSapply(seq(1:100),sim_BCor)
    rank_mat$max_BI4[1:100,]=sfSapply(seq(1:100),sim_maxBI4)
    sfStop()
    toc()
    
    sum(rank_mat$max_BI4[1:100,1]<11)
    sum(rank_mat$max_BI4[,2]<11)
    sum(rowSums(rank_mat$max_BI4[,]<11)>1)
    
    
    tic()
    sfInit(cpus=3,parallel=TRUE) 
    sfExportAll()
    sfLibrary(MASS)
    sfLibrary(movMF)
    rank_mat$max_BCor[1:100,]=sfSapply(seq(1:100),sim_maxBCor)
    sfStop()
    toc()
    
    tic()
    sfInit(cpus=3,parallel=TRUE) 
    sfExportAll()
    sfLibrary(MASS)
    sfLibrary(movMF)
    rank_mat$max_DCor[1:100,]=sfSapply(seq(1:100),sim_maxDCor)
    sfStop()
    toc()
    
    thres=21
    rowSums(rank_mat$nh_BI4<thres)
    sum(colSums(rank_mat$nh_BI4<thres)>2)
    rowSums(rank_mat$nh_BCor<thres)
    sum(colSums(rank_mat$nh_BCor<thres)>2)
    rowSums(rank_mat$nh_DCor<thres)
    sum(colSums(rank_mat$nh_DCor<thres)>2)
    
    rowSums(rank_mat$max_BI4<thres)
    sum(colSums(rank_mat$max_BI4<thres)>2)
    rowSums(rank_mat$max_BCor<thres)
    sum(colSums(rank_mat$max_BCor<thres)>2)
    rowSums(rank_mat$max_DCor<thres)
    sum(colSums(rank_mat$max_DCor<thres)>2)
    save.image(paste0("md32_outlier_",ops,".RData"))
  }else 
    if(Model=="3.3"){
    library(Ball)
    library(BallImpurityFunc)
    library(movMF)
    rep=100
    n=200          #### sample size
    p_X <- 200    #### number of columns in X
    n_sig=3
    p_Y<-5
    #noise_level <- 0.2
    signal_level1 <- 3
    signal_level2 <- 3
    noise_level1 <- 0.2
    noise_level2 <- 0.3
    
    outlier_ratio_vector=c(0.10,0.05,0.01)
    results=NULL#
    for(outlier_ratio in outlier_ratio_vector){
      rank_mat=NULL
      for(r in 1:rep){
        seed=202279+r
        set.seed(seed)
        X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
        sig = sample(1:p_X,size=n_sig,replace=F)
        print(sig)
        Z=X[,sig]*signal_level1+matrix(rnorm(n*n_sig,mean=0,sd=noise_level1),nrow=n)
        transformed_cov=cbind(Z[,1]*Z[,2],Z[,3])
        q_Z=ncol(transformed_cov)
        B_True=matrix(c(0,0,1,2,2,0,0,1,0,0),nrow=2,byrow = F)#c(0,rep(0,p_Y-q_Z-1),rep(1,q_Z))*signal_level2
        mu0=c(1,0,0,0,0)#*signal_level
        B_True_vec <- as.vector(B_True*signal_level2)
        mu0_vec<-as.vector(mu0)
        ##Pheno
        library(movMF)
        Y=matrix(NA,nrow=dim(transformed_cov)[1],ncol=p_Y)
        out_id=sample(1:n,size=floor(n*outlier_ratio))
        library(MASS)
        disturb_sample=mvrnorm(n=n*10,mu=rep(0,5),Sigma = diag(rep(5,5)))
        if(length(which(apply(disturb_sample,1,max)>=5))>=length(out_id)){
          disturb_filtered_id=sample(which(apply(disturb_sample,1,max)>=5),length(out_id))}
        else{disturb_filtered_id=order(apply(disturb_sample,1,max),decreasing=T)[1:length(out_id)]}
        disturb=disturb_sample[disturb_filtered_id,]
        for(i in 1:dim(transformed_cov)[1])
        {
          mu=transformed_cov[i,,drop=FALSE]%*%B_True
          if(i %in% out_id){
            mu=mu+disturb[match(i,out_id),]#mvrnorm(n,mu=rep(0,5),Sigma=diag(rep(5,5)))
          }
          Y[i,]=rmovMF(1,theta=mu)}
        dist_Y=nhdist(Y)
        biTree<- buildTree(list(X=X,Y=Y), leafFunc=regLeaf, method='BallImpurity',D=dist_Y, ops = c(0, 4), maxDepth=4, leafDepth=1, lastSplitFeat=NULL, lossChoice='loss2')
        newTree=PruneTreeComplete(biTree,0.15,depth=4)
        bisig=unique(getnodes(newTree))
        if(length(bisig)>10){bisig=bisig[1:10]}
        rank_matrix_BI=c(as.integer(sig[1]%in%bisig),as.integer(sig[2]%in%bisig),as.integer(sig[3]%in%bisig))
        rank_mat=rbind(rank_mat,c(rank_matrix_BI))
      }
      ratesi=apply((rank_mat<=10)&(rank_mat>0),2,sum)
      ratesallbi=sum(apply(rank_mat[,(1:n_sig)]>0,1,sum)>=3)
      results=rbind(results,c(ratesi,ratesallbi))
    }
    save.image("md3.3.RData")
  }
}
