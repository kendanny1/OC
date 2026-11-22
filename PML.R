rm(list=ls())



#parameters
nl=388L #number of labeled data
nu=0L #number of unlabeled data
nd=nl+nu
KCe=4
p=7L    #number of attributes
mu1=rep(0,p) #mean vector of group 2
mu2=rep(1.5,p) #mean vector of group 1
mu3=rep(c(1.5,0),p/2)
mu123=c(mu1,mu2,mu3)
sig=matrix(0.3,p,p)+diag(0.7,p,p) #covariance matrix
si=rep(sig[upper.tri(sig,diag=TRUE)],3)
'|' <- function(a,b)ifelse(a==0 & b==0, 0, a/b)
K=function(u){
  out=0.75*(1-u^2)*(abs(u)<=1)
  return(out)
}
K4=function(u){
  (1-3*u^2)*(1-u^2)^2*(abs(u)<=1)*105/64
}
cols2=function(x){
  if(length(x)==p){
    out=x
  }else{
    out=colSums(x)
  }
  return(out)
}
cols3=function(x){
  if(length(x)==p){
    out=sum(x)
  }else{
    out=colSums(x)
  }
  return(out)
}
mo=function(x){
  if(length(x)==p){
    out=matrix(0,p,p)
  }else{
    out=x%*%t(x)
  }
  return(out)
}
t2=function(x){
  if(length(x)==p){
    out=x
  }else{
    out=t(x)
  }
  return(out)
}
simt=1000
mm=sm=pmm=psm=pgm=pom=conv=rep(0,simt)
mmm=mpmm=rep(0,KCe*p)
msm=mpsm=mpgm=mpom=rep(0,p*(p+1)/2)
bd=matrix(0,simt,4)
ncp=rep(0,simt)
est_sigma=matrix(0,simt,p*(p+1)/2)
est_mu=matrix(0,simt,p*KCe)
A=rep(0,simt)
iii=388
k2=5
k=1
y_plot=seq(0,30,length=3000)
fy_plot=matrix(0,simt,3000)
tim=50
P=matrix(0,simt,KCe)
lfne=c()
j1=1
load(paste("/Users/kendanny/Downloads/Optimal_Classification/application/nonzero_hsep_",KCe/2,".RData",sep=""))
KCe=4
#load data
Cl=cxse[1,]+diabetes_data[,9]*KCe/2
Xl=as.matrix(scale(df))
p=7L    #number of attributes

#order data
D=cbind(Xl,Cl)
D=D[order(D[,p+1],decreasing=FALSE),]
pl=outer(D[,p+1],1:KCe,FUN="==")
dim(pl)=c(nl,KCe)
nc=colSums(pl)
Nc=rep(0,KCe+1)
Nc[2:(KCe+1)]=cumsum(nc)

#initial labeled estimator
Cmu=c()
for (i in 1:KCe){
  Dc=D[D[,p+1]==i,1:p]
  Cmu=c(Cmu,cols2(Dc)/nc[i])
}
Csig=matrix(0,p,p)
for(i in 1:KCe){
  Dc=D[D[,p+1]==i,1:p]
  mc=cols2(Dc)/nc[i]
  Csig=Csig+mo(t2(Dc)-mc)
}
Csig=Csig/nl



Y=c()
for(i in 1:KCe){ 
  Dc=D[D[,p+1]==i,1:p]
  ylc=colSums((solve(Csig)%*%(t2(Dc)-Cmu[((i-1)*p+1):(i*p)]))*(t2(Dc)-Cmu[((i-1)*p+1):(i*p)]))
  Y=c(Y,(ylc^(p/2)+1)^(2/p)-1)
}


#AMISE criterion
y_seq=seq(min(Y),max(Y),length=5000)
l=function(h){
  K_1=K(outer(Y,Y,FUN="-")/h)/h+K(outer(Y,Y,FUN="+")/h)/h
  dim(K_1)=c(nl,nl)
  diag(K_1)=0
  A_2=2*sum(K_1)/(nl*(nl-1))
  K_2=K(outer(y_seq,Y,FUN="-")/h)/h+K(outer(y_seq,Y,FUN="+")/h)/h
  f_hat_h_y_seq=rowSums(K_2)/nl
  K_3=nl/(nl-1)*f_hat_h_y_seq-K_2/(nl-1)
  out=sum(K_3^2)*(max(Y)-min(Y))/(nl*5000)-A_2
  return(out)
}
h_seq=seq(0.1,sd(Y),length=100)
l_seq=c()
for(i in 1:100) l_seq=c(l_seq,l(h_seq[i]))
h_hat0=nlminb(h_seq[which.min(l_seq)],l)$par
bd[j1,1]=h_hat0

#initial estimate of density generator
g_hat0=function(y){
  K_4=K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="-")/h_hat0)/h_hat0+K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="+")/h_hat0)/h_hat0
  out=rowSums(K_4)*(y^(p/2)+1)^(2/p-1)*gamma(p/2)/(nl*pi^(p/2))
  out
}



cho=chol(Csig)
ini=c(Cmu,cho[upper.tri(cho,diag=TRUE)],nc/nl)
#groupwise pseudo likelihood


#pseudo likelihood estimator
V2=matrix(0,p,p)
V2[upper.tri(V2,diag=TRUE)]=ini[(p*KCe+1):(p*KCe+p*(p+1)/2)]
pssi=t(V2)%*%V2
detp=det(pssi)^(-0.5)

pmu=Y2=c()
for(i in 1:KCe){
  pmu=c(pmu,ini[((i-1)*p+1):(i*p)])
}

for(i in 1:KCe){ 
  Dc=D[D[,p+1]==i,1:p]
  ylc=colSums((solve(pssi)%*%(t2(Dc)-pmu[(p*(i-1)+1):(p*i)])*(t2(Dc)-pmu[(p*(i-1)+1):(p*i)])))
  Y2=c(Y2,(ylc^(p/2)+1)^(2/p)-1)
}




ri3=0
Yt=c()
for(i in 1:KCe){  #Yic for unlabeled
  Yt=cbind(Yt,cols3((solve(pssi)%*%(t(Xl)-pmu[(p*(i-1)+1):(p*i)]))*(t(Xl)-pmu[(p*(i-1)+1):(p*i)])))
}
pp=matrix(0,nd,KCe)
for (i in 1:KCe){
  pp[,i]=detp*g_hat0(Yt[,i])*ini[KCe*p+p*(p+1)/2+i]
}
Clee=c()
for(i in 1:nl){
  cc=(diabetes_data[i,9]==1)*KCe/2
  a=which(pp[i,]==max(pp[i,(cc+1):(cc+KCe/2)]))
  if(length(a)>1){
    b=sample(a,1,replace=FALSE,prob=rep(1,length(a)))
  }
  if(length(a)==1){
    b=a
  }
  Clee=c(Clee,b)
}
l1c=0
while (sum(Cl==Clee)<nl) {
  Cl=Clee
  Xl=as.matrix(scale(df))
  p=7L    #number of attributes
  
  #order data
  D=cbind(Xl,Cl)
  D=D[order(D[,p+1],decreasing=FALSE),]
  pl=outer(D[,p+1],1:KCe,FUN="==")
  dim(pl)=c(nl,KCe)
  nc=colSums(pl)
  Nc=rep(0,KCe+1)
  Nc[2:(KCe+1)]=cumsum(nc)
  
  #initial labeled estimator
  Cmu=c()
  for (i in 1:KCe){
    Dc=D[D[,p+1]==i,1:p]
    Cmu=c(Cmu,cols2(Dc)/nc[i])
  }
  Csig=matrix(0,p,p)
  for(i in 1:KCe){
    Dc=D[D[,p+1]==i,1:p]
    mc=cols2(Dc)/nc[i]
    Csig=Csig+mo(t2(Dc)-mc)
  }
  Csig=Csig/nl
  
  
  
  Y=c()
  for(i in 1:KCe){ 
    Dc=D[D[,p+1]==i,1:p]
    ylc=colSums((solve(Csig)%*%(t2(Dc)-Cmu[((i-1)*p+1):(i*p)]))*(t2(Dc)-Cmu[((i-1)*p+1):(i*p)]))
    Y=c(Y,(ylc^(p/2)+1)^(2/p)-1)
  }
  
  
  h_hat0=h_hat0*nl^(3/80)
  
  #initial estimate of density generator
  g_hat0=function(y){
    K_4=K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="-")/h_hat0)/h_hat0+K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="+")/h_hat0)/h_hat0
    out=rowSums(K_4)*(y^(p/2)+1)^(2/p-1)*gamma(p/2)/(nl*pi^(p/2))
    out
  }
  
  
  
  cho=chol(Csig)
  ini=c(Cmu,cho[upper.tri(cho,diag=TRUE)],nc/nl)
  #groupwise pseudo likelihood
  
  Ls=function(th){
    g_hat0=function(y){
      K_4=K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="-")/h_hat0)/h_hat0+K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="+")/h_hat0)/h_hat0
      out=rowSums(K_4)*(y^(p/2)+1)^(2/p-1)*gamma(p/2)/(nl*pi^(p/2))
      out
    }
    inir=ini
    inir[(p*KCe+2):(p*KCe+p*(p+1)/2)]=th
    
    for(i in 1:KCe){
      mu=inir[((i-1)*p+1):(i*p)]
      siv=inir[(p*KCe+1):(p*KCe+p*(p+1)/2)]
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      U=t(solve(V))
      Dc=D[D[,p+1]==i,1:p]
      Y[(Nc[i]+1):(Nc[i+1])]=(cols3((U%*%(t2(Dc)-mu))^2)^(p/2)+1)^(2/p)-1
    }
    lfyu=0
    for(i in 1:(KCe-1)){
      mu=ini[((i-1)*p+1):(i*p)]
      siv=inir[(p*KCe+1):(p*KCe+p*(p+1)/2)]
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      U=t(solve(V))
      Dc=D[D[,p+1]==i,1:p]
      Yl=cols3((U%*%(t2(Dc)-mu))^2)
      lf=sum(log(g_hat0(Yl)))+nc[i]*(sum(0.5*log(diag(U)^2)))
      lfyu=lfyu+lf
    }  
    mu=ini[((KCe-1)*p+1):(KCe*p)]
    siv=inir[(p*KCe+1):(p*KCe+p*(p+1)/2)]
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    U=t(solve(V))
    Dc=D[D[,p+1]==KCe,1:p]
    Yl=cols3((U%*%(t2(Dc)-mu))^2)
    lf=sum(log(g_hat0(Yl)))+nc[KCe]*(sum(0.5*log(diag(U)^2)))
    lfyu=lfyu+lf
    -lfyu/nl
  }
  
  L1=function(th){
    g_hat0=function(y){
      K_4=K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="-")/h_hat0)/h_hat0+K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="+")/h_hat0)/h_hat0
      out=rowSums(K_4)*(y^(p/2)+1)^(2/p-1)*gamma(p/2)/(nl*pi^(p/2))
      out
    }
    inir=ini
    inir[c(1:p)]=th
    vvv=ini[p*KCe+p*(p+1)/2+KCe+1]
    for(i in 1:KCe){
      mu=inir[((i-1)*p+1):(i*p)]
      siv=inir[(p*KCe+1):(p*KCe+p*(p+1)/2)]
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      U=t(solve(V))
      Dc=D[D[,p+1]==i,1:p]
      Y[(Nc[i]+1):(Nc[i+1])]=(cols3((U%*%(t2(Dc)-mu))^2)^(p/2)+1)^(2/p)-1
    }
    lfyu=0
    mu=th[1:p]
    siv=ini[(p*KCe+1):(p*KCe+p*(p+1)/2)]
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    U=t(solve(V))
    Dc=D[D[,p+1]==1,1:p]
    Yl=cols3((U%*%(t2(Dc)-mu))^2)
    lf=sum(log(g_hat0(Yl)))+nc[1]*(sum(0.5*log(diag(U)^2)))
    lfyu=lfyu+lf
    mu=ini[((KCe-1)*p+1):(KCe*p)]
    siv=ini[(p*KCe+1):(p*KCe+p*(p+1)/2)]
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    U=t(solve(V))
    Dc=D[D[,p+1]==KCe,1:p]
    Yl=cols3((U%*%(t2(Dc)-mu))^2)
    lf=sum(log(g_hat0(Yl)))+nc[KCe]*(sum(0.5*log(diag(U)^2)))
    lfyu=lfyu+lf
    -lfyu/nl
  }
  
  L2=function(th){
    g_hat0=function(y){
      K_4=K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="-")/h_hat0)/h_hat0+K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="+")/h_hat0)/h_hat0
      out=rowSums(K_4)*(y^(p/2)+1)^(2/p-1)*gamma(p/2)/(nl*pi^(p/2))
      out
    }
    inir=ini
    inir[(p+1):(2*p)]=th
    for(i in 1:KCe){
      mu=inir[((i-1)*p+1):(i*p)]
      siv=inir[(p*KCe+1):(p*KCe+p*(p+1)/2)]
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      U=t(solve(V))
      Dc=D[D[,p+1]==i,1:p]
      Y[(Nc[i]+1):(Nc[i+1])]=(cols3((U%*%(t2(Dc)-mu))^2)^(p/2)+1)^(2/p)-1
    }
    lfyu=0
    mu=th[1:p]
    siv=ini[(p*KCe+1):(p*KCe+p*(p+1)/2)]
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    U=t(solve(V))
    Dc=D[D[,p+1]==2,1:p]
    Yl=cols3((U%*%(t2(Dc)-mu))^2)
    lf=sum(log(g_hat0(Yl)))+nc[2]*(sum(0.5*log(diag(U)^2)))
    lfyu=lfyu+lf
    -lfyu/nl
  }
  
  
  L3=function(th){
    g_hat0=function(y){
      K_4=K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="-")/h_hat0)/h_hat0+K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="+")/h_hat0)/h_hat0
      out=rowSums(K_4)*(y^(p/2)+1)^(2/p-1)*gamma(p/2)/(nl*pi^(p/2))
      out
    }
    inir=ini
    inir[(2*p+1):(3*p)]=th
    for(i in 1:KCe){
      mu=inir[((i-1)*p+1):(i*p)]
      siv=inir[(p*KCe+1):(p*KCe+p*(p+1)/2)]
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      U=t(solve(V))
      Dc=D[D[,p+1]==i,1:p]
      Y[(Nc[i]+1):(Nc[i+1])]=(cols3((U%*%(t2(Dc)-mu))^2)^(p/2)+1)^(2/p)-1
    }
    lfyu=0
    mu=th[1:p]
    siv=ini[(p*KCe+1):(p*KCe+p*(p+1)/2)]
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    U=t(solve(V))
    Dc=D[D[,p+1]==3,1:p]
    Yl=cols3((U%*%(t2(Dc)-mu))^2)
    lf=sum(log(g_hat0(Yl)))+nc[3]*(sum(0.5*log(diag(U)^2)))
    lfyu=lfyu+lf
    -lfyu/nl
  }
  
  
  
  L4=function(th){
    g_hat0=function(y){
      K_4=K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="-")/h_hat0)/h_hat0+K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="+")/h_hat0)/h_hat0
      out=rowSums(K_4)*(y^(p/2)+1)^(2/p-1)*gamma(p/2)/(nl*pi^(p/2))
      out
    }
    inir=ini
    inir[(3*p+1):(4*p)]=th
    for(i in 1:KCe){
      mu=inir[((i-1)*p+1):(i*p)]
      siv=inir[(p*KCe+1):(p*KCe+p*(p+1)/2)]
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      U=t(solve(V))
      Dc=D[D[,p+1]==i,1:p]
      Y[(Nc[i]+1):(Nc[i+1])]=(cols3((U%*%(t2(Dc)-mu))^2)^(p/2)+1)^(2/p)-1
    }
    lfyu=0
    mu=th[1:p]
    siv=ini[(p*KCe+1):(p*KCe+p*(p+1)/2)]
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    U=t(solve(V))
    Dc=D[D[,p+1]==4,1:p]
    Yl=cols3((U%*%(t2(Dc)-mu))^2)
    lf=sum(log(g_hat0(Yl)))+nc[4]*(sum(0.5*log(diag(U)^2)))
    lfyu=lfyu+lf
    -lfyu/nl
  }
  lfo=0
  lfn=10
  l4c=0
  while(abs(lfo-lfn)>0.001){
    lfo=lfn
    ini[(p*KCe+2):(p*KCe+p*(p+1)/2)]=nlminb(ini[(p*KCe+2):(p*KCe+p*(p+1)/2)],Ls,control = list(trace=1,rel.tol=1e-3,x.tol=5e-4,iter.max=7))$par
    ini[p*KCe+p*(p+1)/2+KCe]=1-sum(ini[(p*KCe+p*(p+1)/2+1):(p*KCe+p*(p+1)/2+KCe-1)])
    
    ini[c(1:p)]=nlminb(ini[c(1:p)],L1,control = list(trace=1,rel.tol=1e-3,x.tol=5e-4,iter.max=7))$par
    ini[p*KCe+p*(p+1)/2+KCe]=1-sum(ini[(p*KCe+p*(p+1)/2+1):(p*KCe+p*(p+1)/2+KCe-1)])
    ini[c((p+1):(2*p))]=nlminb(ini[c((p+1):(2*p))],L2,control = list(trace=1,rel.tol=1e-3,x.tol=5e-4,iter.max=7))$par
    ini[p*KCe+p*(p+1)/2+KCe]=1-sum(ini[(p*KCe+p*(p+1)/2+1):(p*KCe+p*(p+1)/2+KCe-1)])
    ini[c((2*p+1):(3*p))]=nlminb(ini[c((2*p+1):(3*p))],L3,control = list(trace=1,rel.tol=1e-3,x.tol=5e-4,iter.max=7))$par
    ini[p*KCe+p*(p+1)/2+KCe]=1-sum(ini[(p*KCe+p*(p+1)/2+1):(p*KCe+p*(p+1)/2+KCe-1)])
    ini[(3*p+1):(4*p)]=nlminb(ini[(3*p+1):(4*p)],L4,control = list(trace=1,rel.tol=1e-3,x.tol=5e-4,iter.max=7))$par
    ini[p*KCe+p*(p+1)/2+KCe]=1-sum(ini[(p*KCe+p*(p+1)/2+1):(p*KCe+p*(p+1)/2+KCe-1)])
    inir=ini
    for(i in 1:KCe){
      mu=inir[((i-1)*p+1):(i*p)]
      siv=inir[(p*KCe+1):(p*KCe+p*(p+1)/2)]
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      U=t(solve(V))
      Dc=D[D[,p+1]==i,1:p]
      Y[(Nc[i]+1):(Nc[i+1])]=(cols3((U%*%(t2(Dc)-mu))^2)^(p/2)+1)^(2/p)-1
    }
    lfn=0
    lfnl=0
    Yl=Yl2=c()
    for (i in 1:KCe){
      mu=ini[((i-1)*p+1):(i*p)]
      siv=ini[(p*KCe+1):(p*KCe+p*(p+1)/2)]
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      U=t(solve(V))
      Dc=D[D[,p+1]==i,1:p]
      Yl=c(Yl,cols3((U%*%(t2(Dc)-mu))^2))
      Yl2=c(Yl2,(cols3((U%*%(t2(Dc)-mu))^2)^(p/2)+1)^(2/p)-1)
    }
    K_1=K(outer(Yl2,Yl2,FUN="-")/h_hat0)/h_hat0+K(outer(Yl2,Yl2,FUN="+")/h_hat0)/h_hat0
    dim(K_1)=c(length(Yl2),length(Yl2))
    diag(K_1)=0
    out=rowSums(K_1)*(Yl^(p/2)+1)^(2/p-1)*gamma(p/2)/(nl*pi^(p/2))
    lfn=lfn+sum(log(out[out!=0]))+sum(out!=0)*(sum(0.5*log(diag(U)^2)))+sum(out==0)*(log(min(out[out!=0]))+sum(0.5*log(diag(U)^2)))
    lfnl=lfnl+sum(out!=0)
    lfn=-lfn/nl
    l4c=l4c+1
    if(l4c>=15) break
  }
  h_hat0=h_hat0*nl^(-3/80)
  lfni=lfn
  
  #pseudo likelihood estimator
  V2=matrix(0,p,p)
  V2[upper.tri(V2,diag=TRUE)]=ini[(p*KCe+1):(p*KCe+p*(p+1)/2)]
  pssi=t(V2)%*%V2
  detp=det(pssi)^(-0.5)
  
  pmu=Y2=c()
  for(i in 1:KCe){
    pmu=c(pmu,ini[((i-1)*p+1):(i*p)])
  }
  
  for(i in 1:KCe){ 
    Dc=D[D[,p+1]==i,1:p]
    ylc=colSums((solve(pssi)%*%(t2(Dc)-pmu[(p*(i-1)+1):(p*i)])*(t2(Dc)-pmu[(p*(i-1)+1):(p*i)])))
    Y2=c(Y2,(ylc^(p/2)+1)^(2/p)-1)
  }
  ri3=0
  Yt=c()
  for(i in 1:KCe){  #Yic for unlabeled
    Yt=cbind(Yt,cols3((solve(pssi)%*%(t(Xl)-pmu[(p*(i-1)+1):(p*i)]))*(t(Xl)-pmu[(p*(i-1)+1):(p*i)])))
  }
  pp=matrix(0,nd,KCe)
  for (i in 1:KCe){
    pp[,i]=detp*g_hat0(Yt[,i])*ini[KCe*p+p*(p+1)/2+i]
  }
  Clee=c()
  for(i in 1:nl){
    cc=(diabetes_data[i,9]==1)*KCe/2
    a=which(pp[i,]==max(pp[i,(cc+1):(cc+KCe/2)]))
    if(length(a)>1){
      b=sample(a,1,replace=FALSE,prob=rep(1,length(a)))
    }
    if(length(a)==1){
      b=a
    }
    Clee=c(Clee,b)
  } 
  l1c=l1c+1
  if(l1c>=20) break
}







r_seq=seq(0,((max(Y)+1)^(p/2)-1)^(2/p),length=5000)
A=sum(r_seq^(p+1)*g_hat0(r_seq^2)*2*pi^(p/2)/gamma(p/2)*((max(Y)+1)^(p/2)-1)^(2/p)/5000/p)
save.image("/Users/kendanny/Downloads/Optimal_Classification/application/dependent_2.RData")
load("/Users/kendanny/Downloads/Optimal_Classification/application/dependent_boot_2.RData")
load("/Users/kendanny/Downloads/Optimal_Classification/application/dependent_2.RData")
Cl=Clee
#order data
D=cbind(Xl,Cl)
D=D[order(D[,p+1],decreasing=FALSE),]
pl=outer(D[,p+1],1:KCe,FUN="==")
dim(pl)=c(nl,KCe)
nc=colSums(pl)
Nc=rep(0,KCe+1)
Nc[2:(KCe+1)]=cumsum(nc)
AUC=c()
par(mfrow = c(1, 3),mai=c(1.02,0.82,0.82,0.42))
for(i in 1:(KCe/2)){
  D1=Xl[Clee==i,]
  D2=Xl[Clee==i+KCe/2,]
  Y11=cols3((solve(pssi)%*%(t(D1)-pmu[(p*(i-1)+1):(p*i)]))*(t(D1)-pmu[(p*(i-1)+1):(p*i)]))
  Y12=cols3((solve(pssi)%*%(t(D1)-pmu[(KCe/2*p+p*(i-1)+1):(KCe/2*p+p*i)]))*(t(D1)-pmu[(KCe/2*p+p*(i-1)+1):(KCe/2*p+p*i)]))
  Y21=cols3((solve(pssi)%*%(t(D2)-pmu[(p*(i-1)+1):(p*i)]))*(t(D2)-pmu[(p*(i-1)+1):(p*i)]))
  Y22=cols3((solve(pssi)%*%(t(D2)-pmu[(KCe/2*p+p*(i-1)+1):(KCe/2*p+p*i)]))*(t(D2)-pmu[(KCe/2*p+p*(i-1)+1):(KCe/2*p+p*i)]))
  S1=g_hat0(Y12)/g_hat0(Y11)
  S2=g_hat0(Y22)/g_hat0(Y21)
  S2[S2==Inf]=max(c(S2[S2!=Inf],S1))+1
  S_seq=sort(c(-0.1,S1,S2,max(S1,S2)+1))
  
  TPR=colSums(outer(S2,S_seq,FUN=">"))/nc[i+KCe/2]
  FPR=colSums(outer(S1,S_seq,FUN=">"))/nc[i]
  In_s=outer(S1,S2,FUN="<")
  ROC_cl2=c()
  for(j in 1:1001){
    ROC_cl2=c(ROC_cl2,TPR[min(which(FPR<=ROC_seq[j]))])
  }
  ROC_cl=c(ROC_cl,ROC_cl2)
  AUC=c(AUC,sum(In_s)/( (nc[i]/nl) * (nc[i+KCe/2]/nl) * nl^2))
  plot(0, 0, type = "l", col = "blue", lwd = 2,
       xlim = c(0, 1), ylim = c(0, 1),
       xaxs = "i", yaxs = "i", 
       axes = FALSE,
       xlab = "FPR",
       ylab = "TPR",
       main = paste("Cluster",i,sep=" "),
       font.lab=2,
       cex.main = 2,  # 放大標題
       cex.lab = 1.45)   # 放大軸標籤
  
  # 自訂軸刻度，每 0.1 一格 + 放大字體
  axis(1, at = seq(0, 1, by = 0.1),
       labels = ifelse(seq(0, 1, by = 0.1) == 0, "", seq(0, 1, by = 0.1)),
       cex.axis = 2)
  axis(2, at = seq(0, 1, by = 0.1), labels = ifelse(seq(0, 1, by = 0.1) == 0, "", seq(0, 1, by = 0.1)),
       cex.axis = 2)
  text(0, 0, "0", pos = 1, offset = 0.6, xpd = TRUE, cex = 2)
  
  ROC_cl4=ROC_cl2+R_0.95[i]*apply(ROC_TPR[,((i-1)*1001+1):(i*1001)],MARGIN=2,FUN="sd")
  ROC_cl4[ROC_cl4>1]=1
  ROC_cl5=ROC_cl2-R_0.95[i]*apply(ROC_TPR[,((i-1)*1001+1):(i*1001)],MARGIN=2,FUN="sd")
  ROC_cl5[ROC_cl5<0]=0
  # 加上外框
  lines(FPR,TPR,col="black",lwd=2,lty=1)
  lines(ROC_seq,ROC_cl4,col="black",lwd=2,lty=2)
  lines(ROC_seq,ROC_cl5,col="black",lwd=2,lty=2)
}

AUC2=c()
S11=S12=S21=S22=0
for(i in 1:(KCe/2)){
  D1=Xl[diabetes_data[,9]==0,]
  D2=Xl[diabetes_data[,9]==1,]
  Y11=cols3((solve(pssi)%*%(t(D1)-pmu[(p*(i-1)+1):(p*i)]))*(t(D1)-pmu[(p*(i-1)+1):(p*i)]))
  Y12=cols3((solve(pssi)%*%(t(D1)-pmu[(KCe/2*p+p*(i-1)+1):(KCe/2*p+p*i)]))*(t(D1)-pmu[(KCe/2*p+p*(i-1)+1):(KCe/2*p+p*i)]))
  Y21=cols3((solve(pssi)%*%(t(D2)-pmu[(p*(i-1)+1):(p*i)]))*(t(D2)-pmu[(p*(i-1)+1):(p*i)]))
  Y22=cols3((solve(pssi)%*%(t(D2)-pmu[(KCe/2*p+p*(i-1)+1):(KCe/2*p+p*i)]))*(t(D2)-pmu[(KCe/2*p+p*(i-1)+1):(KCe/2*p+p*i)]))
  S11=S11+g_hat0(Y11)*nc[i]/nl
  S12=S12+g_hat0(Y12)*nc[i+KCe/2]/nl
  S21=S21+g_hat0(Y21)*nc[i]/nl
  S22=S22+g_hat0(Y22)*nc[i+KCe/2]/nl
  
}
S1=S12/S11
S2=S22/S21
S2[S2==Inf]=max(S2[S2!=Inf])
S_seq=sort(c(-0.1,S1,S2,max(S1,S2)+1))
TPR=colSums(outer(S2,S_seq,FUN=">"))/sum(nc[3:4])
FPR=colSums(outer(S1,S_seq,FUN=">"))/sum(nc[1:2])
In_s=outer(S1,S2,FUN="<")
ROC_cl2=c()
for(j in 1:1001){
  ROC_cl2=c(ROC_cl2,TPR[min(which(FPR<=ROC_seq[j]))])
}
ROC_cl4=ROC_cl2+R_0.95[3]*apply(ROC_TPR[,((3-1)*1001+1):(3*1001)],MARGIN=2,FUN="sd")
ROC_cl4[ROC_cl4>1]=1
ROC_cl5=ROC_cl2-R_0.95[3]*apply(ROC_TPR[,((3-1)*1001+1):(3*1001)],MARGIN=2,FUN="sd")
ROC_cl5[ROC_cl5<0]=0
plot(0, 0, type = "l", col = "blue", lwd = 2,
     xlim = c(0, 1), ylim = c(0, 1),
     xaxs = "i", yaxs = "i", 
     axes = FALSE,
     xlab = "FPR",
     ylab = "TPR",
     main = "All Participants",
     font.lab=2,
     cex.main = 2,  # 放大標題
     cex.lab = 1.45)   # 放大軸標籤

# 自訂軸刻度，每 0.1 一格 + 放大字體
axis(1, at = seq(0, 1, by = 0.1),
     labels = ifelse(seq(0, 1, by = 0.1) == 0, "", seq(0, 1, by = 0.1)),
     cex.axis = 2)
axis(2, at = seq(0, 1, by = 0.1), labels = ifelse(seq(0, 1, by = 0.1) == 0, "", seq(0, 1, by = 0.1)),
     cex.axis = 2)
text(0, 0, "0", pos = 1, offset = 0.6, xpd = TRUE, cex = 2)

# 加上外框

# 再把 ROC 曲線疊上去（畫在軸上層）
lines(FPR,TPR,col="black",lwd=2,lty=1)
lines(ROC_seq,ROC_cl4,col="black",lwd=2,lty=2)
lines(ROC_seq,ROC_cl5,col="black",lwd=2,lty=2)
AUC2=sum(In_s)/( (sum(nc[1:2])/nl) * (sum(nc[3:4])/nl) * nl^2 )


save.image(paste("/Users/kendanny/Downloads/Optimal_Classification/application/dependent_",KCe/2,".RData",sep=""),version = 2)

