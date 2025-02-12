ImpRedFunc=function(Y,X,D,d,subsample=FALSE,rep=10){
  ### This function computes the impurity reduction (or equivalently, purity gain) of Y 
  ### corresponding to a covariate/feature X.
  ### Y is a list or vector of responses
  ### X is a vector of covariates
  ### D is the distance matrix of all the Y values
  ### d is the parameter in BI
  ### optional: subsample method can be called, with default number of rep=10
  ### default is to use the full sample
  ### output is a vector with length n_distinct(X)-1, corresponding to all the possible splits
  source("/code/splitDataset.R")
  source("/code/BImp_bin.R")
  N=length(X)
  all_var=unique(X)
  if(length(all_var)==1){  ##only one variant present
    break
    return(0)
  }
  n.cut=length(all_var)-1
  dataset=list(X=as.matrix(X),Y=Y)
  n.cut=length(all_var)-1 ##possible cuts
  cutpoints=sort(all_var,decreasing=TRUE)
  print(cutpoints)
  if(subsample==FALSE){
    BI_full=BImp_bin(N,D,1:N,1:N,d)
    PG.vec=rep(0,n.cut)###prepare the PG vector, each entry corresponding to a cut
    for(cut in 1:n.cut){
      split=splitDataset(dataset,ops1=1,ops2=cutpoints[cut])
      m1=length(split$leftInd)
      m2=N-m1
      BI_left=BImp_bin(N,D,1:N,split$leftInd,d)
      BI_right=BImp_bin(N,D,1:N,split$rightInd,d)
      PG_k=BI_full-m1/n*BI_left-m2/n*BI_right
      PG.vec[cut]<-PG_k}
  }else{
    BI_full=BImp_bin(N,D,1:N,1:N,d)
    PG.sub=matrix(rep(0,n.cut*rep),nrow=rep)  ##prepare the PG.sub matrix, each row is one replicate
    for(replicate in 1:rep){###replicate the subsample computation for 10 different subsamples
      ### and take the average PurityGain
      sub.id=NULL         ##start to select the subsample
      #print(all_var)
      set.seed(replicate+1000)
      
      sub.id=sample(1:N,size=ceiling(N/subsize),replace=F)
      print(sub.id)
      if(length(unique(snp.info[sub.id]))<length(all_var)){###if not all variants are selected, re-sample until it includes all variants
        resamp=0
        while(length(unique(snp.info[sub.id]))<length(all_var)){
          set.seed(replicate+1000+resamp)
          sub.id=sample(1:N,size=ceiling(N/subsize),replace=F)
          resamp=resamp+1
        }
      }
      #print(sub.id)
      n.sub=length(sub.id) 
      D=distance_matrix[sub.id,sub.id] ######get the n by N distance matrix
      snp.sub=snp.info[sub.id]
      pheno.sub=pheno.info[sub.id]
      dataset=list(X=as.matrix(snp.sub),Y=pheno.sub)
      print("begin BImp calculation:")
      BI_fullsubsample=BImp_bin(n.sub,D,1:n.sub,1:n.sub,d)
      print(paste0("completed BImp calculation, BImp full subsample equals ",BI_fullsubsample,"."))
      print(Sys.time())
      n.cut=length(all_var)-1 ##possible cuts
      cutpoints=sort(all_var,decreasing=TRUE)
      for(cut in 1:n.cut){
        split=splitDataset(dataset,ops1=1,ops2=cutpoints[cut])
        m1=length(split$leftInd)
        m2=n.sub-m1
        print(c(m1,m2))
        print("starts BI_left")
        BI_left=BImp_subsample(n.sub,D,1:n.sub,split$leftInd,d)#BallImpurity(n,distance_matrix,split$leftInd,d)
        print("BI_left:")
        print(BI_left)
        print("starts BI_right")
        BI_right=BImp_subsample(n.sub,D,1:n.sub,split$rightInd,d)#BallImpurity(n,distance_matrix,split$rightInd,d)
        print("BI_right:")
        print(BI_right)
        PG_k=BI_fullsubsample-m1/n.sub*BI_left-m2/n.sub*BI_right
        print(PG_k)
        PG.sub[replicate,cut]<-PG_k
      }
    }
    PG.vec=apply(PG.sub,2,mean) ## take the average of all replicates as the final PG vector
  }
 return(PG.vec)
}