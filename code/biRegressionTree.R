library(Rcpp)
sourceCpp("/code/BallImpurity.cpp")
sourceCpp("/code/BallDistanceVector.cpp")

splitDataset <- function(dataset, ops1, ops2) {
  # split dataset into left and right child nodes
  # input: 
  #   dataset contains covariates X and outcome Y
  #   When X is Euclidean，ops1 corresponds to the index of the feature to be used for split (X[, ops1]).
  #   ops2 is the feature value
  #  
  # output: left and right child nodes info
  leftSet <- list(); rightSet <- list()
  if (sum(class(dataset$X) == 'matrix') >= 1) {
    feat <- ops1; val <- ops2
    leftInd <- which(dataset$X[,feat] < val); rightInd <- which(dataset$X[,feat] >= val)
    leftSet$X <- dataset$X[leftInd,]; rightSet$X <- dataset$X[rightInd,]
  } else if (class(dataset$X) == 'array') {
    candidateVar <- ops1; D <- ops2
    n <- nrow(D); c1 <- candidateVar[1]; c2 <- candidateVar[2]; cInds <- c(1:n)[-c(c1,c2)]
    leftInd <- c(c1); rightInd <- c(c2)
    for (i in cInds) { if (D[c1, i] < D[c2, i]) {leftInd <- cbind(leftInd, i)} else {rightInd <- cbind(rightInd, i)} }
    leftSet$X <- dataset$X[,,leftInd,]; rightSet$X <- dataset$X[,,rightInd,]
    
    if (length(dim(leftInd)) == 1) {leftSet$X <- array(leftSet$X, c(dim(dataset$X)[1], dim(dataset$X)[2], 1, dim(dataset$X)[4]))}
    if (length(dim(rightInd)) == 1) {rightSet$X <- array(rightSet$X, c(dim(dataset$X)[1], dim(dataset$X)[2], 1, dim(dataset$X)[4]))}
    
    if (length(dim(leftSet$X)) < 4) {leftSet$X <- array(leftSet$X, c(dim(dataset$X)[1], dim(dataset$X)[2], dim(leftSet$X)[3], dim(dataset$X)[4]))}
    if (length(dim(rightSet$X)) < 4) {rightSet$X <- array(rightSet$X, c(dim(dataset$X)[1], dim(dataset$X)[2], dim(rightSet$X)[3], dim(dataset$X)[4]))}
  }
  
  if (sum(class(dataset$Y) == 'matrix') >= 1) {
    leftSet$Y <-dataset$Y[leftInd,]; rightSet$Y <-dataset$Y[rightInd,]
  } else if (class(dataset$Y) == 'numeric') {
    leftSet$Y <-dataset$Y[leftInd]; rightSet$Y <-dataset$Y[rightInd]
  } else if (class(dataset$Y) == 'array') {
    leftSet$Y <-dataset$Y[,,leftInd]; rightSet$Y <-dataset$Y[,,rightInd]
    if (length(dim(leftInd)) == 1) {leftSet$Y <- array(leftSet$Y, c(dim(dataset$Y)[1], dim(dataset$Y)[2], 1))}
    if (length(dim(rightInd)) == 1) {rightSet$Y <- array(rightSet$Y, c(dim(dataset$Y)[1], dim(dataset$Y)[2], 1))}
  } else if (class(dataset$Y) == 'list') {
    leftSet$Y=dataset[[2]][leftInd]
    rightSet$Y=dataset[[2]][rightInd]
    
  } 
  return(list(leftSet=leftSet, rightSet=rightSet, leftInd=leftInd, rightInd=rightInd))
  
}
PruningStat1<-function(D,ind1){
  ## the pruning statistic, in the spirit of ANOVA
  ## D is the distance matrix of parent node set
  ## ind1 is the indices of either child node set in the parent node set, the rest forms another child node in a binary tree
  ## ind2 is the indices of the other child node set
  SA=sum(D^2)-sum(diag(D^2))
  D1=D[ind1,ind1]
  S1=sum(D1^2)-sum(diag(D1^2))
  ind2=setdiff(seq(1:nrow(D)),ind1)
  D2=D[ind2,ind2]
  S2=sum(D2^2)-sum(diag(D2^2))
  df1=choose(nrow(D),2)
  numerator=SA/df1
  df2=(choose(length(ind1),2)+choose(nrow(D)-length(ind1),2))
  denom=(S1+S2)/df2
  fval=numerator/denom
  stat=1-pf(fval,df1,df2)
  return(stat)
}

regLeaf <- function(target, method = 'BallImpurity') {
  # 
  # input: Y values
  # output: mean values of Y in the leaf
  if (method == 'BallImpurity') {
    if (sum(class(target) == 'matrix') >= 1) { return(apply(target, 2, median)) }
    if (class(target) == 'numeric') {
      if (length(unique(target)) <= 2) {return(as.numeric(names(table(target))[table(target) == max(table(target))])[1])}
      else { return(median(target)) }
    }
    
    if (class(target) == 'array') { 
      n <- dim(target)[3]
      if (n == 1) {return(target)}
      temp <- lapply(1:n, function(i, Y_=target){Y_[,,i]})
      temp <- list(data = temp, size = dim(target[,,1]), name = 'landmark')
      class(temp) <- 'riemdata'
      temp <- km_median(temp) # 'intrinsic', 'extrinsic'
      return(temp$median)
    }
    if (class(target) == 'list') {
      if (length(target) == 1) {return(target[[1]])}
      temp <- list(data = target, size = dim(target[[1]]), name = 'landmark')
      class(temp) <- 'riemdata'
      temp <- km_median(temp) # 'intrinsic', 'extrinsic'
      return(temp$median)
    }
    
  } else {
    if (sum(class(target) == 'matrix') >= 1) {return(colMeans(target))}
    if (class(target) == 'numeric') {
      if (length(unique(target)) <= 2) {return(as.numeric(names(table(target))[table(target) == max(table(target))])[1])}
      else { return(mean(target)) }
    }
    if (class(target) == 'array') { return(apply(simplify2array(target), 1:2, mean)) }
    if (class(target) == 'list') {
      n <- length(target); k <- nrow(target[[1]]); m <- ncol(target[[1]])
      tmp <- array(as.numeric(unlist(target)), dim=c(k, m, n))
      return(apply(simplify2array(tmp), 1:2, mean))
    }
  }
}
bestSplit <- function(dataset, D, leafFunc=regLeaf, method='BallImpurity', ops=c(0, 4), lastSplitFeat=NULL) {
  # finding the best cut variable and cut value by CART 
  # input: dataset, leafFunc- the method to compute the leaf; ops - parameters specifying the tolerances on loss and N,  
  # lastSplitFeat - the node of last split
  # output: the best cut variable and cut value; if no split is performed, return the leaf value
  tolLoss <- ops[1]; tolN <- ops[2]
  X <- dataset$X; Y <- dataset$Y
  if (sum(class(Y) == 'matrix') >= 1) {tmpLen <- dim(unique(Y, MARGIN = 1))[1]
  }else if ((class(Y) == 'numeric') || (class(Y) == 'list')) {tmpLen <- length(unique(Y))
  }else if (class(Y) == 'array') {tmpLen <- dim(unique(Y, MARGIN = 3))[3]}
  if (tmpLen == 1) {return(leafFunc(Y, method = method))}

  n <- dim(D)[1]; d <- median(D)
  
  if (sum(class(X) == 'matrix') >= 1) {
    n <- dim(X)[1]; p <- dim(X)[2] #Euclidean
  } else if (class(X) == 'array') {
    p <- dim(X)[length(dim(X))]; n <- dim(X)[3]; D_X <- list() # non-Euclidean
    for (i in 1:p) { D_X[[i]] <- shapeDist(X[,,,i]) }
  }
  
  YHat <- leafFunc(Y, method = method); loss <-  BallImpurity(n,D,1:n,d)
  chosen_feature <- c(1:p)
 
  ## Euclidean 
  if (sum(class(X) == 'matrix') >= 1) {
    bestFeat <- chosen_feature[1]; bestVal <- X[1,bestFeat]; bestLoss <- 1e5
    for (i in 1:length(chosen_feature)) {
      feat <- chosen_feature[i]; x <- X[, feat]
      print(paste0("now computing feature ",feat))
      xSort <- sort(x); xSort <- unique(xSort)
      if (length(xSort) > 1) {
        cutVal <- (xSort[1:(length(xSort)-1)] + xSort[2:length(xSort)] )/2
        for (j in 1:length(cutVal)) {
          tmpSubset <- splitDataset(dataset, feat, cutVal[j])
          leftSet <- tmpSubset[[1]]; rightSet <- tmpSubset[[2]]
          leftInd <- tmpSubset[[3]]; rightInd <- tmpSubset[[4]]
          leftLen <- length(leftInd); rightLen <- length(rightInd)
          leftY <- leftSet$Y; rightY <- rightSet$Y
          if ((leftLen >= tolN) && (rightLen >= tolN)) {
            
            D_left <- D[leftInd, leftInd]; D_right <- D[rightInd, rightInd]
            nowLoss <- leftLen/n*BallImpurity(n,D,leftInd,d) + rightLen/n*BallImpurity(n,D,rightInd,d)
            
            canSplit <- TRUE
          } else {canSplit <- FALSE}
          
          if ((canSplit) && (nowLoss < bestLoss)) {
            bestLoss <- nowLoss; bestFeat <- feat; bestVal <- cutVal[j]
          }
        }
      }
      print(paste0("best feature is ",bestFeat))
    }
    tmpSubset <- splitDataset(dataset, bestFeat, bestVal)
    leftSet <- tmpSubset[[1]]; rightSet <- tmpSubset[[2]]
    leftInd <- tmpSubset[[3]]; rightInd <- tmpSubset[[4]]
    leftLen <- length(leftInd); rightLen <- length(rightInd)
    leftY <- leftSet$Y; rightY <- rightSet$Y
    if ((loss - bestLoss) < tolLoss) {print("loss not improving");return(YHat)}

    if ((leftLen < tolN) || (rightLen < tolN)) {print("min leaf size reached");return(YHat)}
    splitstat1  <- PruningStat1(D,leftInd)
    return(c(bestFeat, bestVal, bestLoss, splitstat1))
  }
  
  ##############
  ## Non_Euclidean X
  if (length(class(X))==1&&class(X)=="array") {
    bestFeat <- chosen_feature[1]; bestLoss <- 1e5; bestD_X <- D_X[[bestFeat]]
    #print(bestFeat)
    candidateVars <- which(bestD_X == max(bestD_X), arr.ind = TRUE)
    candidateVars <- getUniqueInd(candidateVars); bestCand <- candidateVars[1,]
    bestVal <- list(ind = bestCand, c1 = X[,,,bestFeat][,,bestCand[1]], c2 = X[,,,bestFeat][,,bestCand[2]]) 
    nCand <- nrow(candidateVars)
    
    for (i in 1:length(chosen_feature)) {
      feat <- chosen_feature[i]; tmpD_X <- D_X[[feat]]
      tempDx <- unique(as.numeric(tmpD_X)); tempDx <- tempDx[-which(tempDx == 0)]
      if (length(unique(X[,,,feat])) > 1){
        for (td in tempDx) {
          candidateVars <- which(tmpD_X == td, arr.ind = TRUE)
          candidateVars <- getUniqueInd(candidateVars)
          nCand <- nrow(candidateVars)
          for (j in 1:nCand) {
            candVar <- candidateVars[j,]
            tmpSubset <- splitDataset(dataset, candVar, tmpD_X)
            leftSet <- tmpSubset[[1]]; rightSet <- tmpSubset[[2]]
            leftInd <- tmpSubset[[3]]; rightInd <- tmpSubset[[4]]
            leftLen <- length(leftInd); rightLen <- length(rightInd)
            leftY <- leftSet$Y; rightY <- rightSet$Y
            if ((leftLen >= tolN) && (rightLen >= tolN)) {
              D_left <- D[leftInd, leftInd]; D_right <- D[rightInd, rightInd]
              nowLoss <- leftLen/n*BallImpurity(n,D,leftInd,d) + rightLen/n*BallImpurity(n,d,rightInd,d)
              canSplit <- TRUE
            } else {canSplit <- FALSE}
            
            if ((canSplit) && (nowLoss < bestLoss)) {
              bestLoss <- nowLoss; bestFeat <- feat; bestD_X <- tmpD_X
              bestVal <- list(ind = candVar, c1 = X[,,,bestFeat][,,candVar[1]], c2 = X[,,,bestFeat][,,candVar[2]]) 
            }
          }
        }
      }
      
    }
    
    if ((loss - bestLoss) < tolLoss) {print("loss not improving");return(YHat)}
    
    tmpSubset <- splitDataset(dataset, bestVal$ind, bestD_X)
    leftSet <- tmpSubset[[1]]; rightSet <- tmpSubset[[2]]
    leftInd <- tmpSubset[[3]]; rightInd <- tmpSubset[[4]]
    leftLen <- length(leftInd); rightLen <- length(rightInd)
    leftY <- leftSet$Y; rightY <- rightSet$Y
    splitstat1  <- PruningStat1(D,leftInd)
    if ((leftLen < tolN) || (rightLen < tolN)) {print("leaf size criterion met");return(YHat)}
    return(list(bestFeat, bestVal, bestLoss, splitstat1))
  }
  
}
buildTree <- function(dataset, leafFunc=regLeaf, method='BallImpurity',D, ops = c(0, 4), maxDepth=4, leafDepth=1, lastSplitFeat=NULL) {
 
  # input: 
  # input: dataset, leafFunc- the method to compute the leaf; method - loss function; ops - parameters specifying the tolerances on loss and N,  lastSplitFeat - the node for last split
  #   since this is non-Euclidean data, need to supply distance matrix D.
  #   maxDepth - max depth of the tree; leafDepth - depth of current node，lastSplitFeat - feature for last split
  # output:generate
  
  nodeset<-NULL
  X <- dataset$X; Y <- dataset$Y
  splitResult <- bestSplit(dataset, D, leafFunc, method, ops, lastSplitFeat)
  if ((sum(class(X) == 'array')==1) && (sum(class(X)=='matrix')!=1)) {n <- dim(X)[3]} else {n <- dim(X)[1]}
  
  n <- dim(D)[1]; d <-median(D)

  loss <-  BallImpurity(n,D,1:n,d)
  if (length(splitResult)>=3){if(length(splitResult[[2]])>1) {flag<-TRUE} else {flag<-FALSE}}
  if ((length(splitResult)==1) || (!((length(splitResult)>=3)&&flag) && (class(splitResult) != "numeric")) || (leafDepth >= maxDepth)) {
    yHat <- leafFunc(dataset$Y, method = method)
    resTree <- list()
    resTree$depth <- leafDepth
    resTree$nsample <- n
    resTree$leafValue <- yHat
    resTree$impur <- n * loss
    return(resTree)
  }
  
  if (length(class(X))==1&&class(X)=="array") {
    feat <- splitResult[[1]]; val <- splitResult[[2]]; impur <- splitResult[[3]]
  } else {
    feat <- splitResult[1]; val <- splitResult[2]; impur <- splitResult[3]; stat1<-splitResult[4]
  }
  resTree <- list()
  resTree$depth <- leafDepth
  resTree$nsample <- n
  resTree$splitInd <- feat
  resTree$splitVal <- val
  resTree$splitImpurity <- impur
  resTree$splitstat1<-stat1
  resTree$impur <- n * loss

  
  # Tree is split into left and right branches, and the splitting continues
  if (length(class(dataset$X))==1&&class(dataset$X)=="array") {
    tmpD_X <- shapeDist(dataset$X[,,,feat]); candVar <- val$ind
    subDataset <-splitDataset(dataset, candVar, tmpD_X)
  } else {subDataset <- splitDataset(dataset, feat, val)}
  lSet <- subDataset[[1]]; lD=D[subDataset$leftInd,subDataset$leftInd]; rSet <- subDataset[[2]]; rD=D[subDataset$rightInd,subDataset$rightInd]
  
  resTree$right <- buildTree(rSet,leafFunc, method,rD, ops, maxDepth, leafDepth+1, feat)
  resTree$left <- buildTree(lSet, leafFunc, method,lD, ops, maxDepth, leafDepth+1, feat)
  return(resTree)
}

getnodes<-function(tree){
  nodes<-NULL
  subtreelist<-list(tree)
  while(length(subtreelist[[1]])>0){
    nodes<-c(nodes,subtreelist[[1]]$splitInd)
    subtreelist<-c(subtreelist,list(subtreelist[[1]]$left,subtreelist[[1]]$right))
    subtreelist[[1]]<-NULL
  }
  return(nodes)
}
PruneTreeLeaf<-function(tree,thres){
  if(length(tree)==4){### leaf node reached
    print("No branches!")
    return(tree)
  }else if((length(tree$left)==4)&(length(tree$right)==4)){####the last branch
    flag=as.numeric(tree$splitstat1>thres)
    if(flag==1){### the last branch needs to be merged
      tree$splitInd<-NULL;tree$splitVal<-NULL;tree$splitImpurity<-NULL;tree$splitstat1<-NULL;
      tree$leafValue=(tree$right$leafValue*tree$right$nsample+tree$left$leafValue*tree$left$nsample)/tree$nsample
      tree$left<-NULL;tree$right<-NULL}
    return(tree)
  }else{ ### not the last branch, keep pruning. 
    tree$left=PruneTreeLeaf(tree$left,thres); tree$right=PruneTreeLeaf(tree$right,thres)
    return(tree)}
}
PruneTreeComplete<-function(tree,thres,depth){
  for(d in 1:depth){
    tree=PruneTreeLeaf(tree,thres)
  }
  return(tree)
}
