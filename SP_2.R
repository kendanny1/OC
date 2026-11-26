
rm(list=ls())
library(openxlsx)
library(expm)

for(KCnn in 2:6){
  result2=c()
  risum=c()
  result1=c()
  KCs1=c()
  cxse=c()
  
  # Define a function to calculate column sums
  cols2=function(x){
    if(length(x)==p){
      out=x
    }else{
      out=colSums(x)
    }
    return(out)
  }
  
  SSS_1128=c()
  cxk=c()
  KCs=c()
  ris=c()
  
  for(KC1 in 6:KCnn){
    
    # Check whether the corresponding result file exists
    if(file.exists(paste("nonzero_sep_",KC1,"_",KCnn,".RData",sep=""))==TRUE){
      
      # Load the result file from the current working directory
      load(paste("nonzero_sep_",KC1,"_",KCnn,".RData",sep=""))
      
      # Automatically determine the number of observations and variables
      n=nrow(df)
      p=ncol(df)
      
      X=scale(as.matrix(df))
      Csig3 <- matrix(0, p, p)
      bgamma <- matrix(0, 2*p, KCnn)
      
      # Calculate cluster-specific mean vectors
      for (d in 0:1) {
        for (k in 1:KCnn) {
          idx <- which(cxp == k & D == d)
          if (length(idx) == 0) {
            bgamma[(d*p+1):((d+1)*p), k] <- NA
          } else if (length(idx) == 1) {
            bgamma[(d*p+1):((d+1)*p), k] <- X[idx, ]
          } else {
            bgamma[(d*p+1):((d+1)*p), k] <- colMeans(X[idx, , drop = FALSE])
          }
        }
      }
      
      # Replace missing values with the global mean of each label
      for (d in 0:1) {
        idx_label <- which(D == d)
        global_mu <- colMeans(X[idx_label, , drop = FALSE])
        for (k in 1:KCnn) {
          mu_slot <- bgamma[(d*p+1):((d+1)*p), k]
          if (any(is.na(mu_slot))) bgamma[(d*p+1):((d+1)*p), k] <- global_mu
        }
      }
      
      # Calculate the initial covariance matrix Csig3
      Csig3 <- matrix(0, p, p)
      for (d in 0:1) {
        for (k in 1:KCnn) {
          idx <- which(cxp == k & D == d)
          if (length(idx) == 0) next
          mu_kd <- bgamma[(d*p+1):((d+1)*p), k]
          diffs <- sweep(X[idx, , drop = FALSE], 2, mu_kd)
          Csig3 <- Csig3 + t(diffs) %*% diffs
        }
      }
      Csig3 <- Csig3 / n
    }
  }
  
  for(KC1 in KCnn:6){
    
    # Check whether the corresponding result file exists
    if(file.exists(paste("nonzero_sep_",KC1,"_",KCnn,".RData",sep=""))==TRUE){
      
      # Load the result file from the current working directory
      load(paste("nonzero_sep_",KC1,"_",KCnn,".RData",sep=""))
      
      beta=gamma=c()
      kappa=1
      
      # Automatically determine the number of observations and variables
      n=nrow(df)
      p=ncol(df)
      
      X=scale(as.matrix(df))
      
      # Recalculate bgamma and Csig
      bgamma <- matrix(0, 2*p, KCnn)
      for (d in 0:1) {
        for (k in 1:KCnn) {
          idx <- which(cxp == k & D == d)
          if (length(idx) == 0) {
            bgamma[(d*p+1):((d+1)*p), k] <- NA
          } else if (length(idx) == 1) {
            bgamma[(d*p+1):((d+1)*p), k] <- X[idx, ]
          } else {
            bgamma[(d*p+1):((d+1)*p), k] <- colMeans(X[idx, , drop = FALSE])
          }
        }
      }
      
      # Replace missing values with the global mean of each label
      for (d in 0:1) {
        idx_label <- which(D == d)
        global_mu <- colMeans(X[idx_label, , drop = FALSE])
        for (k in 1:KCnn) {
          if (any(is.na(bgamma[(d*p+1):((d+1)*p), k])))
            bgamma[(d*p+1):((d+1)*p), k] <- global_mu
        }
      }
      
      for(i in 1:n) beta=c(beta,bgamma[1:p, cxp[i]])
      for(i in 1:n) beta=c(beta,bgamma[(p+1):(2*p), cxt[i]]-bgamma[1:p, cxp[i]])
      
      # Calculate the residual matrix and the objective function
      XXX=X2-A1%*%beta
      dim(XXX)=c(p,n)
      Sn=sum(XXX*(solve(Csig3)%*%XXX))/(2*n)
      
      SSS_1128=c(SSS_1128,Sn)
      cxk=rbind(cxk,cxp)
      KCs=c(KCs,KC1)
    }
  }
  
  # Select the clustering result with the minimum objective value
  if(length(SSS_1128)>0){
    cxs=cxk[which.min(SSS_1128),]
    cxse=rbind(cxse,cxk[which.min(SSS_1128),])
    rri=0
    risum=c(risum,rri)
    KCs1=c(KCs1,KCs[which.min(SSS_1128)])
    rm(Xlg)
  }
  
  # Save the workspace to the current working directory
  save.image(paste("nonzero_hsep_",KCnn,".RData",sep=""),version = 2)
}