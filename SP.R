rm(list=ls())
library(expm)

# Soft-thresholding function
ST = function(x, y) {
  out = numeric(length(x))
  nz = abs(x) > 0
  out[nz] = (1 - y / abs(x[nz])) * (abs(x[nz]) - y > 0) * x[nz]
  out[!nz] = 0
  return(out)
}

# Calculate column sums or return the original vector
cols2=function(x){
  if(length(x)==p){
    out=x
  }else{
    out=colSums(x)
  }
  return(out)
}

# Calculate row sums or return the original vector
rows2=function(x){
  if(length(x)==p){
    out=x
  }else{
    out=rowSums(x)
  }
  return(out)
}

for(KCn in 2:6){
  for(KC in KCn:6){
    go=c()
    gn=c()
    risp2=c()
    wm=c()
    
    # Load the corresponding clustering results from the current working directory
    load(paste("km_",KC,".RData",sep=""))
    
    # Automatically determine the number of observations (n) and variables (p)
    X=#Markers
    D=#outcome
    n=nrow(X)
    p=ncol(X)
    
    Cm=matrix(0,n,n*p)
    c_1210=gamma_1210=MSEkmsp_1210=MSEhesp_1210=MSEkmsig_1210=MSEhesig_1210=MSEgw_1210=c()
    
    A1=cbind(diag(n*p),bdiag(lapply(D, function(x) (x == 1) * diag(p))))
    A2=kronecker(diag(2),kronecker(kronecker(diag(n),rep(1,KCn)),diag(p)))
    A3=kronecker(diag(2),kronecker(rep(1,n),diag(KCn*p)))
    
    CC=rbind(
      cbind(
        bdiag(lapply(D, function(x) (x == 0) * ( KCn/(KCn^2+KCn) ) * diag(p) + (x == 1) * ( (KCn+1) /(KCn^2+2*KCn) ) * diag(p) )),
        bdiag(lapply(D, function(x) -(x == 1)/(KCn^2+2*KCn) * diag(p)))
      ),
      cbind(
        bdiag(lapply(D, function(x) -(x == 1)/(KCn^2+2*KCn) * diag(p))),
        bdiag(lapply(D, function(x) (x == 0) * ( 1/KCn ) * diag(p) + (x == 1) * ( (KCn+1) /(KCn^2+2*KCn) ) * diag(p) ))
      )
    )
    
    A4=solve(rbind(
      cbind(
        -kronecker(matrix(1,KCn,KCn),(sum(D==0)*KCn/(KCn^2+KCn)+sum(D==1)*((KCn+1)/(KCn^2+2*KCn)))*diag(p)) + n*diag(KCn*p),
        -kronecker(matrix(1,KCn,KCn),(-sum(D==1)/(KCn^2+2*KCn))*diag(p))
      ),
      cbind(
        -kronecker(matrix(1,KCn,KCn),(-sum(D==1)/(KCn^2+2*KCn))*diag(p)),
        -kronecker(matrix(1,KCn,KCn),(sum(D==0)*1/KCn+sum(D==1)*((KCn+1)/(KCn^2+2*KCn)))*diag(p)) + n*diag(KCn*p)
      )
    ))
    
    A5=t(A3) %*% A2 %*% CC %*% t(A1)
    A6=t(A3) %*% A2 %*% CC %*% t(A2)
    
    ril=c()
    msel=c()
    S2=c()
    SS=c()
    cxx=c()
    betat=c()
    gammat=c()
    cxt=Clek
    
    # Merge KC clusters into KCn clusters
    if(KCn<KC){
      X=#markers
      X[D==0,]=t(t(X[D==0,])-colSums(X[D==0,])/sum(D==0))
      X[D==1,]=t(t(X[D==1,])-colSums(X[D==1,])/sum(D==1))
      
      for(i in 1:(KC-KCn)){
        gammas=c()
        for(ij in 1:(KC-i+1)) gammas=cbind(gammas,cols2(X[cxt==ij,])/sum(cxt==ij))
        subs=c()
        for(ij in 1:(KC-i)){
          for(jk in (ij+1):(KC-i+1)) subs=cbind(subs,c(sum((gammas[,ij]-gammas[,jk])^2),ij,jk))
        }
        ij1=subs[2,which.min(subs[1,])]
        jk1=subs[3,which.min(subs[1,])]
        cxt[cxt==jk1]=ij1
        if(jk1<KC-i+1){
          cxt[cxt>jk1]=cxt[cxt>jk1]-1
        }
      } 
    }
    
    beta=gamma=c()
    kappa=1
    X=#markers %*%invCsig
    
    # Recalculate bgamma (cluster-specific mean vectors)
    bgamma <- matrix(0, 2*p, KCn)
    for (d in 0:1) {
      for (k in 1:KCn) {
        idx <- which(cxt == k & D == d)
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
      for (k in 1:KCn) {
        if (any(is.na(bgamma[(d*p+1):((d+1)*p), k])))
          bgamma[(d*p+1):((d+1)*p), k] <- global_mu
      }
    }
    
    for(i in 1:n) beta=c(beta,bgamma[1:p, cxt[i]])
    for(i in 1:n) beta=c(beta,bgamma[(p+1):(2*p), cxt[i]]-bgamma[1:p, cxt[i]])
    
    gamma=c(c(bgamma[1:p,]),c(bgamma[(p+1):(2*p),]))
    gamma[(KCn*p+1):(2*KCn*p)]=gamma[(KCn*p+1):(2*KCn*p)]-gamma[1:(KCn*p)]
    
    delta=A2%*%beta-A3%*%gamma
    X2=t(X)
    dim(X2)=c(n*p,1)
    Sn=sum(((X2-A1%*%beta))^2)/(2*n)
    SS=c(SS,Sn)
    cxx=cbind(cxx,cxt)
    
    lam_max=max(abs(delta))
    lam_seq=2^seq(log(lam_max,2)-10,log(lam_max,2),length=11)
    
    for(lambda in lam_seq){
      cxn=cxt
      cxo=rep(0,n)
      l0c=0
      
      while(sum(cxo==cxn)<n){
        l0c=l0c+1
        beta=gamma=c()
        cxo=cxn
        
        for(i in 1:KCn){
          if(sum(cxo==i)==0){
            cxo[cxo>i]=cxo[cxo>i]-1
          }
        }
        
        # Recalculate bgamma
        bgamma <- matrix(0, 2*p, KCn)
        for (d in 0:1) {
          for (k in 1:KCn) {
            idx <- which(cxo == k & D == d)
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
          for (k in 1:KCn) {
            if (any(is.na(bgamma[(d*p+1):((d+1)*p), k])))
              bgamma[(d*p+1):((d+1)*p), k] <- global_mu
          }
        }
        
        for(i in 1:n) beta=c(beta,bgamma[1:p, cxo[i]])
        for(i in 1:n) beta=c(beta,bgamma[(p+1):(2*p), cxo[i]]-bgamma[1:p, cxo[i]])
        
        gamma=c(c(bgamma[1:p,]),c(bgamma[(p+1):(2*p),]))
        gamma[(KCn*p+1):(2*KCn*p)]=gamma[(KCn*p+1):(2*KCn*p)]-gamma[1:(KCn*p)]
        
        delta=A2%*%beta-A3%*%gamma
        nu=rep(0,2*KCn*n*p)
        betan=beta
        gamman=gamma
        deltan=delta
        nun=nu
        met=10
        mett=c()
        obj_i=c()
        objo=Inf
        objn=1e+10
        
        # ADMM optimization
        while(met>n*1e-4*KCn && objo>objn ){
          met=0
          beta=betan
          gamma=gamman
          delta=deltan
          objo=objn
          nu=nun
          
          # Gamma update
          gamman=A4 %*% (-t(A3) %*% (delta-nu)+ A5%*% X2 + A6 %*% (delta-nu))
          
          # Beta update
          for (i in 1:n) {
            nui0=nu[((i-1)*KCn*p+1):((i*KCn*p))]
            nui1=nu[(n*KCn*p+(i-1)*KCn*p+1):(n*KCn*p+(i*KCn*p))]
            nusum=c(rowSums(matrix(nui0,p,KCn)),rowSums(matrix(nui1,p,KCn)))
            
            deltai0=delta[((i-1)*KCn*p+1):((i*KCn*p))]
            deltai1=delta[(n*KCn*p+(i-1)*KCn*p+1):(n*KCn*p+(i*KCn*p))]
            deltasum=c(rowSums(matrix(deltai0,p,KCn)),rowSums(matrix(deltai1,p,KCn)))
            
            gammasum = c(
              rowSums(matrix(gamman[1:(p*KCn)], p, KCn)),
              rowSums(matrix(gamman[(p*KCn+1):(2*p*KCn)], p, KCn))
            )
            
            betan[c(((i-1)*p+1):(i*p),((i-1)*p+1+n*p):(i*p+n*p))]=solve(
              t(cbind(diag(p),(D[i]==1)*diag(p))) %*% cbind(diag(p),(D[i]==1)*diag(p))+ KCn*diag(2*p)
            ) %*% (
              t(cbind(diag(p),(D[i]==1)*diag(p))) %*% X[i,] -nusum +gammasum+deltasum
            )
          }
          
          # Delta update
          deltan20=deltan21=c()
          obj_d=c()
          
          for(i in 1:n) {
            De=c()
            deltai20=deltai21=c()
            
            for(j in 1:KCn){
              deltail0=betan[((i-1)*p+1):(i*p)] - gamman[((j-1)*p+1):(j*p)] +nu[((i-1)*KCn*p+(j-1)*p+1):((i-1)*KCn*p+j*p)]
              deltail1=betan[((i-1)*p+1+n*p):(i*p+n*p)] - gamman[((j-1)*p+1+KCn*p):(j*p+KCn*p)] +nu[((i-1)*KCn*p+(j-1)*p+1+n*KCn*p):((i-1)*KCn*p+j*p+n*KCn*p)]
              
              De=c(De,lambda*(
                sum(abs(ST(deltail0,lambda))) +
                  sum(abs(ST(deltail1,lambda))) +
                  0.5 * sum((deltail0-ST(deltail0,lambda))^2) +
                  0.5 * sum((deltail1-ST(deltail1,lambda))^2)
              ))
              
              deltai20=c(deltai20,deltail0)
              deltai21=c(deltai21,deltail1)
            }
            
            ell=which.min(De)
            deltai20[((ell-1)*p+1):(ell*p)]=ST(deltai20[((ell-1)*p+1):(ell*p)],lambda)
            deltai21[((ell-1)*p+1):(ell*p)]=ST(deltai21[((ell-1)*p+1):(ell*p)],lambda)
            deltan20=c(deltan20,deltai20)
            deltan21=c(deltan21,deltai21)
            
            obj_d=c(obj_d,lambda*(
              sum(abs(deltai20[((ell-1)*p+1):(ell*p)]))+
                sum(abs(deltai21[((ell-1)*p+1):(ell*p)]))
            ))
          }
          
          deltan=c(deltan20,deltan21)
          
          # Nu update
          for(i in 1:n){
            for(j in 1:KCn){
              nun[((i-1)*KCn*p+(j-1)*p+1):((i-1)*KCn*p+j*p)] =
                nu[((i-1)*KCn*p+(j-1)*p+1):((i-1)*KCn*p+j*p)] +
                betan[((i-1)*p+1):(i*p)] -
                gamman[((j-1)*p+1):(j*p)] -
                deltan[((i-1)*KCn*p+(j-1)*p+1):((i-1)*KCn*p+j*p)]
              
              nun[((i-1)*KCn*p+(j-1)*p+1+n*KCn*p):((i-1)*KCn*p+j*p+n*KCn*p)] =
                nu[((i-1)*KCn*p+(j-1)*p+1+n*KCn*p):((i-1)*KCn*p+j*p+n*KCn*p)] +
                betan[((i-1)*p+1+n*p):(i*p+n*p)] -
                gamman[((j-1)*p+1+KCn*p):(j*p+KCn*p)] -
                deltan[((i-1)*KCn*p+(j-1)*p+1+n*KCn*p):((i-1)*KCn*p+j*p+n*KCn*p)]
              
              met=met+
                sum((betan[((i-1)*p+1):(i*p)] -
                       gamman[((j-1)*p+1):(j*p)] -
                       deltan[((i-1)*KCn*p+(j-1)*p+1):((i-1)*KCn*p+j*p)])^2)+
                sum((betan[((i-1)*p+1+n*p):(i*p+n*p)] -
                       gamman[((j-1)*p+1+KCn*p):(j*p+KCn*p)] -
                       deltan[((i-1)*KCn*p+(j-1)*p+1+n*KCn*p):((i-1)*KCn*p+j*p+n*KCn*p)])^2)
            }
          }
          
          mett=c(mett,met)
          objn=sum(((X2-A1%*%beta))^2)/2+
            sum(nun*(A2%*% betan-A3%*% gamman-deltan))+
            sum(((A2%*% betan-A3%*% gamman-deltan))^2)/2+
            sum(obj_d)
          
          obj_i=c(obj_i,objn)
        }
        
        # Reclassify observations based on the updated parameters
        R=A2%*%beta-A3%*%gamma
        RR=c()
        RR2=c()
        
        for(i in 1:n){
          RRR=c()
          for(j in 1:KCn){
            RRR=c(RRR,
                  sum(abs(R[((i-1)*KCn*p+(j-1)*p+1):((i-1)*KCn*p+j*p)]))+
                    sum(abs(R[(n*KCn*p+(i-1)*KCn*p+(j-1)*p+1):(n*KCn*p+(i-1)*KCn*p+j*p)]))
            )
          }
          RR=c(RR,which.min(RRR))
          RR2=c(RR2,min(RRR))
        }
        
        cxn=RR
        if(l0c>=5) break
      }
      
      beta=gamma=c()
      for(i in 1:n) beta=c(beta,bgamma[1:p, cxo[i]])
      for(i in 1:n) beta=c(beta,bgamma[(p+1):(2*p), cxo[i]]-bgamma[1:p, cxo[i]])
      
      SS=c(SS,sum(((X2-A1%*%beta))^2)/(2*n))
      cxx=cbind(cxx,cxo)
    }
    
    wm=c(wm,which.min(SS))
    cxp=cxx[,which.min(SS)]
    c_1210=rbind(c_1210,cxp)
    
    rm(Xlg)
    
    # Save the workspace to the current working directory
    save.image(paste("nonzero_sep_",KC,"_",KCn,".RData",sep=""),version = 2)
  }
}