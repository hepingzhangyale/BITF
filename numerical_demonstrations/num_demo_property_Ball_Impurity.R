##Property of Ball Impurity
##Menglu Che 20231119


#################################
############ Figure 1 ###########
#################################

rm(list=ls())
library(shapes)
library(Rcpp)
library(Ball)    # for Ball Correlation
library(energy)  # for Distance Correlation

sourceCpp("BallDistanceVector.cpp")
sourceCpp("BallImpurity.cpp")
repeat.time=100
#Gaussion
n=500
p=2
d=0.5
e_seq=c(20:50)/50
BI_seq_mean=rep(0,length(e_seq))
Dcov_seq_mean=rep(0,length(e_seq))
Bcov_seq_mean=rep(0,length(e_seq))
Sum_var_seq_mean=rep(0,length(e_seq))

for(r in 1:repeat.time)
{
  BI_seq=rep(NA,length(e_seq))
  Dcov_seq=rep(NA,length(e_seq))
  Bcov_seq=rep(NA,length(e_seq))
  Sum_var_seq=rep(NA,length(e_seq))
  for(i in 1:length(e_seq))
  {
    e=e_seq[i] #noise level sou
    Y=matrix(rnorm(n*p,mean=0, sd=e),nrow=n)
    #plot(Y)
    distance_matrix=BallDistanceVector(Y)
    BI_seq[i]=BallImpurity(n,distance_matrix,1:n,d)
    Dcov_seq[i]=dcov(Y,Y)
    Bcov_seq[i]=bcov(Y,Y)
    Sum_var_seq[i]=var(Y[,1])+var(Y[,2])
  }
  BI_seq_mean=BI_seq_mean+BI_seq
  Dcov_seq_mean=Dcov_seq_mean+Dcov_seq
  Bcov_seq_mean=Bcov_seq_mean+Bcov_seq
  Sum_var_seq_mean=Sum_var_seq_mean+Sum_var_seq
}

result=cbind(BI_seq_mean,Bcov_seq_mean, Dcov_seq_mean,Sum_var_seq_mean)
result=as.data.frame(scale(result))
min.y=min(result)
max.y=max(result)
pdf("Simulation3_nov24.pdf",width=7,height=7)
plot(e_seq,result$BI_seq_mean,pch=1,type="b",ylim=c(min.y,max.y),xlab="noise",ylab="",yaxt="n",cex.lab=2.5,cex.axis=2.5)
lines(e_seq,result$Bcov_seq_mean,pch=2,type="b",col="blue")
lines(e_seq,result$Dcov_seq_mean,pch=3,type="b",col="red")
lines(e_seq,result$Sum_var_seq_mean,pch=4,type="b",col="green")
legend("topleft", legend=c("BI","BCov", "DCov","SumVar"),
       col=c("black","blue", "red","green"), lty=1,pch=1:4, cex=2,
       box.lty=0)
dev.off()

##Mixtured model 1 Gaussion
n=500
p=2
Y=matrix(NA,nrow=n,ncol=p)
e=0.1
d=0.2
Pi_seq=seq(0,1,0.05)
n_seq=Pi_seq*n
BI_seq_mean=rep(0,length(Pi_seq))
Bcov_seq_mean=rep(0,length(Pi_seq))
Dcov_seq_mean=rep(0,length(Pi_seq))
Sum_var_seq_mean=rep(0,length(Pi_seq))

for(r in 1:repeat.time)
{
  BI_seq=rep(NA,length(Pi_seq))
  Bcov_seq=rep(NA,length(Pi_seq))
  Dcov_seq=rep(NA,length(Pi_seq))
  Sum_var_seq=rep(NA,length(Pi_seq))
  #### sample 2n Gaussian distributed points, and truncate 
  Y1=matrix(rep(0,n*4),nrow=2*n)
  Y2=matrix(rep(0,n*4),nrow=2*n)
  mu1=c(2,2)
  mu2=c(2,-2)
  for(j in 1:2){
    Y1[,j]=rnorm(2*n,mean=mu1[j],sd=e)
    Y2[,j]=rnorm(2*n,mean=mu2[j],sd=e)
  }
  d1=apply(Y1-rep(mu1,each=nrow(Y1)),1,norm,type = "2")
  d2=apply(Y2-rep(mu2,each=nrow(Y2)),1,norm,type = "2")
  id1=which(d1<=0.5)
  id2=which(d2<=0.5)
  for(i in 1:length(Pi_seq))
  {
    n1=n*Pi_seq[i]
    n2=n-n1
    Y.t=rbind(Y1[sample(id1,size=n1),],Y2[sample(id2,size=n2),])
    Y=list()
    for(k in 1:nrow(Y.t)){Y[[k]]<-Y.t[k,]}
    distance_matrix=BallDistanceVector(Y)
    BI_seq[i]=BallImpurity(n,distance_matrix,1:n,d)
    Dcov_seq[i]=dcov(Y.t,Y.t)
    Bcov_seq[i]=bcov(Y.t,Y.t)
    Sum_var_seq[i]=var(Y.t[,1])+var(Y.t[,2])
    #plot(Y)
  }
  print(BI_seq)
  BI_seq_mean=BI_seq_mean+BI_seq
  Bcov_seq_mean=Bcov_seq_mean+Bcov_seq
  Dcov_seq_mean=Dcov_seq_mean+Dcov_seq
  Sum_var_seq_mean=Sum_var_seq_mean+Sum_var_seq
}
result=cbind(BI_seq_mean,Bcov_seq_mean,Dcov_seq_mean,Sum_var_seq_mean)
result=as.data.frame(scale(result))
min.y=min(result)
max.y=max(result)
pdf("Simulation4_nov25.pdf",width=7,height=7)
plot(Pi_seq,result$BI_seq_mean,pch=1,type="b",ylim=c(min.y,max.y),xlab=expression(pi),ylab="",yaxt="n",cex.lab=2.5, cex.axis=2.5)
lines(Pi_seq,result$Bcov_seq_mean,pch=2,type="b",col="blue")
lines(Pi_seq,result$Dcov_seq_mean,pch=3,type="b",col="red")
lines(Pi_seq,result$Sum_var_seq_mean,pch=4,type="b",col="green")
legend("bottom", legend=c("BI","Bcov", "Dcov","SumVar"),
       col=c("black","blue","red","green"), lty=1,pch=1:4, cex=1.8,
       box.lty=0)
dev.off()

##Mixtured model 3 
 n=500
 p=2
 Y=matrix(NA,nrow=n,ncol=p)
 e=0.1
 d=0.2
 Pi_seq=seq(0,pi/2,length.out = 20)
 
 BI_seq_mean=rep(0,length(Pi_seq))
 
 Bcov_seq_mean=rep(0,length(Pi_seq))
 Dcov_seq_mean=rep(0,length(Pi_seq))
 Sum_var_seq_mean=rep(0,length(Pi_seq))
 
 for(r in 1:repeat.time)
 {
   BI_seq=rep(NA,length(Pi_seq))
   Bcov_seq=rep(NA,length(Pi_seq))
   Dcov_seq=rep(NA,length(Pi_seq))
   Sum_var_seq=rep(NA,length(Pi_seq))
   for(i in 1:length(Pi_seq))
   {
     Pi=runif(n,min=0,max=Pi_seq[i])+rbinom(n,1,0.5)*pi
     Y[,1]=cos(Pi)+rnorm(n,mean=0,sd=e)
     Y[,2]=sin(Pi)+rnorm(n,mean=0,sd=e)
     plot(Y,xlim=c(-1,1),ylim=c(-1,1))
     distance_matrix=BallDistanceVector(Y)
     BI_seq[i]=BallImpurity(n,distance_matrix,1:n,d)
     Dcov_seq[i]=dcov(Y,Y)
     Bcov_seq[i]=bcov(Y,Y)
     Sum_var_seq[i]=var(Y[,1])+var(Y[,2])
   }
   BI_seq_mean=BI_seq_mean+BI_seq
   Bcov_seq_mean=Bcov_seq_mean+Bcov_seq
   Dcov_seq_mean=Dcov_seq_mean+Dcov_seq
   Sum_var_seq_mean=Sum_var_seq_mean+Sum_var_seq
 }
 
 
 result=cbind(BI_seq_mean,Bcov_seq_mean,Dcov_seq_mean,Sum_var_seq_mean)
 result=as.data.frame(scale(result))
 min.y=min(result)
 max.y=max(result)
 pdf("Simulation5_nov24.pdf",width=7,height=7)
 plot(Pi_seq,result$BI_seq_mean,pch=1,type="b",ylim=c(min.y,max.y+0.5),xlab=expression(psi),ylab="",yaxt="n",cex.lab=2.5,cex.axis=2.5)
 lines(Pi_seq,result$Bcov_seq_mean,pch=2,type="b",col="blue")  
 lines(Pi_seq,result$Dcov_seq_mean,pch=3,type="b",col="red")
 lines(Pi_seq,result$Sum_var_seq_mean,pch=4,type="b",col="green")
 legend("topleft", legend=c("BI","BCov", "DCov","SumVar"),
        col=c("black","blue", "red","green"), lty=1,pch=1:4, cex=1.8,
        box.lty=0)
 dev.off()

 
##########################
####### Figure 2 #########
##########################
 rm(list=ls())
 library(Rcpp)
 sourceCpp("BallDistanceVector.cpp")
 sourceCpp("BallImpurity.cpp")
 pvec=c(0.1, 0.74, 1.47)
 
 n=500
 rep=100
 sd=0.1
 d=0.4
 set.seed(202303)
 results=NULL
 for(p in pvec){
   pheno_list = list()
   Z=runif(n,0,1)*p+pi*rbinom(n,size=1,prob=0.5)
   Yx=rep(0,n)
   Yy=rep(0,n)
   for(i in 1:n){
     eps=rnorm(2,0,sd)
     #theta=seq(0,2*pi,by=2*pi/(y_length-1))
     x=cos(Z[i])+eps[1]
     y=sin(Z[i])+eps[2]
     Yx[i]=x
     Yy[i]=y
     pheno_list[[i]] = matrix(cbind(x,y),ncol = 2)
   }
   pdf(paste0("cut_points",p,".pdf"),width =5,height=5)
   plot(Yx,Yy,pch=20,xlab="",ylab="",xlim=c(-1.3,1.3),ylim=c(-1.3,1.3),cex.lab=2.2, cex.axis=2.2)
   dev.off()
   
   alphavec=seq(0,pi,by=0.01)
   PG_vec=rep(0,length(alphavec))
   D=BallDistanceVector(pheno_list)#dist(pheno_list)
   for(k in 1:length(alphavec)){
     alpha=alphavec[k]
     set1=which(Yy-tan(alpha)*Yx>=0)
     set2=setdiff(1:n,set1)
     m1=length(set1)
     m2=length(set2)
     if(m1==0 | m2==0){
       PG_vec[k]=0
     }else{
       BI_full=BallImpurity(n,D,1:n,d)
       BI_left=BallImpurity(n,D,set1,d)
       BI_right=BallImpurity(n,D,set2,d)
       PG_vec[k]=BI_full-m1/n*BI_left-m2/n*BI_right
     }}
   results=rbind(results,PG_vec)
 }
 for(j in 1:3){
   pdf(paste0("new_PG_demo_md3_p=",pvec[j],".pdf"),width=5,height=5)
   plot(alphavec,results[j,],xlab=expression(alpha),ylab="",yaxt="n",pch=20,cex.lab=2.2,cex.axis=2)
   lines(alphavec,results[j,], col="red")
   dev.off()
 }
##############################
########## Figure 4 ##########
##############################

 library(movMF)
 library(circular)
 library(rgl)
 
 kappa_vector=c(50,20,10,5,2)
 j=1 ### plot the case for kappa=50. Change the index to get plots for other kappa values
 n_sig=2
 n=200
 p_X=100
 p_Y=3
 signal_level1=0.5;noise_level1=0.05
 set.seed(123)
 X=matrix(rbinom(n*p_X,size=1,prob=0.5),nrow=n)
 #### the original predictors are binary with p=0.5
 sig = sample(1:p_X,size=n_sig,replace=F)
 #### randomly select n_sig of them to be active predictors
 #print(sig)
 Z=X[,sig]+matrix(rnorm(n*n_sig,mean=0,sd=noise_level1),nrow=n)
 ### transform X into Z by magnifying by signal_level1 and adding perturbation 
 ## with sd=noise_level1
 transformed_cov=Z%*%matrix(c(2/3,1/3,-1/3,-1/3),nrow=2,byrow = T)
 ### transform Z into a more compact, 2-dim covariate vector transformed_cov
 ### 
 q_Z=ncol(transformed_cov)
 B_True=matrix(c(0,1,1,0,-1,1),nrow=2,byrow=T)  
 #B_True=c(0,rep(0,p_Y-q_Z-1),rep(1,q_Z))
 mu0=c(1,0,0)#signal_level
 mu0_vec<-as.vector(mu0)
 ##Pheno
 library(movMF)
 Y=matrix(NA,nrow=dim(transformed_cov)[1],ncol=length(B_True))
 for(i in 1:dim(transformed_cov)[1]){
   mu=c(0,transformed_cov[i,,drop = FALSE]) #%*% B_True
   mu=mu+mu0_vec#mu=colSums(t(transformed_cov[i,,drop = FALSE]) %*% t(as.matrix(B_True_vec)))+mu0_vec
   mu1=mu*kappa_vector[j]
   mu_final=mu/norm(mu,type="2")*kappa_vector[j]
   
   Y[i,]=rmovMF(1,mu_final)
 }
 X_vec=rbind(c(0,0),c(0,1),c(1,0),c(1,1))
 open3d()
 spheres3d(0,0,0,radius=1.0,lit=TRUE,color="lavender",front="lines")
 observer3d(0,0,7)
 for(k in 1:4){
   x=Y[which(X[,sig[1]]==X_vec[k,1]&X[,sig[2]]==X_vec[k,2]),]
   print(colMeans(x))
   points3d(x[,1],x[,2],x[,3],col=rainbow(4)[k],mode="point",size=6.5)
   #plot(circular(x[,1]),cex=0.8,bin=50,stack=TRUE,shrink=1.3,main=paste0("kappa=",kappa[i]))
 }
 par3d(windowRect = c(0, 0, 400, 400))
 
 mat <- par3d("userMatrix")
 # Rotate around the 3 axis to better show the points
 mat <- rotate3d(mat, angle = pi/2, x = 1, y = 0, z = 0)
 par3d(userMatrix = mat)
 mat <- rotate3d(mat, angle = pi/2, x = 0, y = -1, z = 0)
 par3d(userMatrix = mat)
 mat <- rotate3d(mat, angle = pi/16, x = 0, y = 0, z = 1)
 par3d(userMatrix = mat)
 
 legend3d("topright", legend = paste('Signal=', c('(0,0)', '(0,1)','(1,0)','(1,1)')), pch = 16, col = rainbow(4), cex=1.4, inset=c(0.01))
 highlevel()
 rgl.snapshot(filename="nov11_kappa50_4groups.png",fmt='png')
