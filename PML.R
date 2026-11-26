# Clear the current R environment before starting the analysis.
rm(list=ls())


# --------------------------------------------------
# Parameters
# --------------------------------------------------

# Set the number of unlabeled observations.
nu=0L

# Set the number of clusters.
KCe=4


# --------------------------------------------------
# Kernel and matrix utility functions
# --------------------------------------------------

# Define the Epanechnikov kernel function.
K=function(u){
  out=0.75*(1-u^2)*(abs(u)<=1)
  return(out)
}

# Calculate column sums for matrix input.
# If the input is a vector, return it directly.
cols2=function(x){
  if(length(x)==p){
    out=x
  }else{
    out=colSums(x)
  }
  return(out)
}

# Calculate the total sum for vector input.
# For matrix input, calculate column sums.
cols3=function(x){
  if(length(x)==p){
    out=sum(x)
  }else{
    out=colSums(x)
  }
  return(out)
}

# Construct a zero matrix for vector input.
# Otherwise, calculate the cross-product matrix.
mo=function(x){
  if(length(x)==p){
    out=matrix(0,p,p)
  }else{
    out=x%*%t(x)
  }
  return(out)
}

# Return a vector unchanged or transpose a matrix.
t2=function(x){
  if(length(x)==p){
    out=x
  }else{
    out=t(x)
  }
  return(out)
}


# --------------------------------------------------
# Load clustering results and data
# --------------------------------------------------

# Load the clustering results from the current working directory.
load(paste("nonzero_hsep_",KCe/2,".RData",sep=""))

# Determine the sample size and number of attributes from the loaded data.
nl=nrow(df) # Number of labeled observations
nd=nl+nu    # Total number of observations
p=ncol(df)  # Number of attributes

KCe=4

# Initialize the cluster labels using the previous clustering results.
Cl=cxse[1,]+D*KCe/2

# Standardize the feature matrix.
Xl=as.matrix(scale(df))


# --------------------------------------------------
# Order observations by cluster labels
# --------------------------------------------------

# Combine the standardized data and cluster labels.
D2=cbind(Xl,Cl)

# Sort observations according to their cluster labels.
D2=D2[order(D2[,p+1],decreasing=FALSE),]

# Construct the cluster membership indicator matrix.
pl=outer(D2[,p+1],1:KCe,FUN="==")
dim(pl)=c(nl,KCe)

# Calculate the number of observations in each cluster.
nc=colSums(pl)

# Calculate the cumulative number of observations in each cluster.
Nc=rep(0,KCe+1)
Nc[2:(KCe+1)]=cumsum(nc)


# --------------------------------------------------
# Initial estimation of cluster means and covariance
# --------------------------------------------------

# Initialize the vector of cluster means.
Cmu=c()

# Estimate the mean vector for each cluster.
for (i in 1:KCe){
  Dc=D2[D2[,p+1]==i,1:p]
  Cmu=c(Cmu,cols2(Dc)/nc[i])
}

# Initialize the pooled covariance matrix.
Csig=matrix(0,p,p)

# Calculate the pooled covariance matrix using all clusters.
for(i in 1:KCe){
  Dc=D2[D2[,p+1]==i,1:p]
  mc=cols2(Dc)/nc[i]
  Csig=Csig+mo(t2(Dc)-mc)
}

# Normalize the pooled covariance matrix by the sample size.
Csig=Csig/nl


# --------------------------------------------------
# Calculate transformed Mahalanobis distances
# --------------------------------------------------

# Initialize the vector of transformed distances.
Y=c()

# Calculate the transformed squared Mahalanobis distances
# for observations in each cluster.
for(i in 1:KCe){ 
  Dc=D2[D2[,p+1]==i,1:p]
  ylc=colSums((solve(Csig)%*%(t2(Dc)-Cmu[((i-1)*p+1):(i*p)]))*(t2(Dc)-Cmu[((i-1)*p+1):(i*p)]))
  Y=c(Y,(ylc^(p/2)+1)^(2/p)-1)
}


# --------------------------------------------------
# Bandwidth selection using the AMISE criterion
# --------------------------------------------------

# Construct a sequence of evaluation points.
y_seq=seq(min(Y),max(Y),length=5000)

# Define the AMISE-based objective function.
l=function(h){
  
  # Calculate the kernel matrix for the observed distances.
  K_1=K(outer(Y,Y,FUN="-")/h)/h+K(outer(Y,Y,FUN="+")/h)/h
  dim(K_1)=c(nl,nl)
  
  # Remove diagonal elements.
  diag(K_1)=0
  
  # Calculate the first component of the criterion.
  A_2=2*sum(K_1)/(nl*(nl-1))
  
  # Calculate the kernel estimates on the evaluation grid.
  K_2=K(outer(y_seq,Y,FUN="-")/h)/h+K(outer(y_seq,Y,FUN="+")/h)/h
  f_hat_h_y_seq=rowSums(K_2)/nl
  
  # Calculate the adjusted kernel estimate.
  K_3=nl/(nl-1)*f_hat_h_y_seq-K_2/(nl-1)
  
  # Evaluate the objective function.
  out=sum(K_3^2)*(max(Y)-min(Y))/(nl*5000)-A_2
  
  return(out)
}

# Define the candidate bandwidth sequence.
h_seq=seq(0.1,sd(Y),length=100)

# Evaluate the AMISE criterion over the candidate bandwidths.
l_seq=c()
for(i in 1:100) l_seq=c(l_seq,l(h_seq[i]))

# Obtain the initial bandwidth estimate.
h_hat0=nlminb(h_seq[which.min(l_seq)],l)$par


# --------------------------------------------------
# Initial estimation of the density generator
# --------------------------------------------------

# Define the initial kernel estimator of the density generator.
g_hat0=function(y){
  
  # Calculate the kernel contributions.
  K_4=K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="-")/h_hat0)/h_hat0+K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="+")/h_hat0)/h_hat0
  
  # Calculate the density estimate.
  out=rowSums(K_4)*(y^(p/2)+1)^(2/p-1)*gamma(p/2)/(nl*pi^(p/2))
  
  out
}


# --------------------------------------------------
# Initialize model parameters
# --------------------------------------------------

# Obtain the Cholesky decomposition of the pooled covariance matrix.
cho=chol(Csig)

# Construct the initial parameter vector.
# It contains cluster means, covariance parameters, and cluster proportions.
ini=c(Cmu,cho[upper.tri(cho,diag=TRUE)],nc/nl)


# --------------------------------------------------
# Initial pseudo-likelihood estimation
# --------------------------------------------------

# Initialize the upper-triangular covariance parameter matrix.
V2=matrix(0,p,p)
V2[upper.tri(V2,diag=TRUE)]=ini[(p*KCe+1):(p*KCe+p*(p+1)/2)]

# Calculate the covariance matrix used in the pseudo-likelihood.
pssi=t(V2)%*%V2

# Calculate the determinant-related constant.
detp=det(pssi)^(-0.5)

# Initialize the cluster mean vector and distance vector.
pmu=Y2=c()

# Extract the initial mean vector for each cluster.
for(i in 1:KCe){
  pmu=c(pmu,ini[((i-1)*p+1):(i*p)])
}

# Calculate the transformed distances under the initial covariance estimate.
for(i in 1:KCe){ 
  Dc=D2[D2[,p+1]==i,1:p]
  ylc=colSums((solve(pssi)%*%(t2(Dc)-pmu[(p*(i-1)+1):(p*i)])*(t2(Dc)-pmu[(p*(i-1)+1):(p*i)])))
  Y2=c(Y2,(ylc^(p/2)+1)^(2/p)-1)
}


# --------------------------------------------------
# Initial cluster assignment
# --------------------------------------------------

# Calculate the distances of all observations from each cluster.
ri3=0
Yt=c()

for(i in 1:KCe){  # Calculate the distance for each cluster.
  Yt=cbind(Yt,cols3((solve(pssi)%*%(t(Xl)-pmu[(p*(i-1)+1):(p*i)]))*(t(Xl)-pmu[(p*(i-1)+1):(p*i)])))
}

# Initialize the pseudo-likelihood matrix.
pp=matrix(0,nd,KCe)

# Calculate the pseudo-likelihood for each observation and cluster.
for (i in 1:KCe){
  pp[,i]=detp*g_hat0(Yt[,i])*ini[KCe*p+p*(p+1)/2+i]
}

# Initialize the updated cluster assignments.
Clee=c()

# Assign each labeled observation to the most likely eligible cluster.
for(i in 1:nl){
  
  # Determine the eligible cluster range.
  cc=(D[i]==1)*KCe/2
  
  # Find the cluster with the largest pseudo-likelihood.
  a=which(pp[i,]==max(pp[i,(cc+1):(cc+KCe/2)]))
  
  # Randomly select a cluster if multiple clusters share the maximum value.
  if(length(a)>1){
    b=sample(a,1,replace=FALSE,prob=rep(1,length(a)))
  }
  
  # Assign the cluster directly if there is only one maximum.
  if(length(a)==1){
    b=a
  }
  
  Clee=c(Clee,b)
}


# --------------------------------------------------
# Iterative estimation procedure
# --------------------------------------------------

# Initialize the iteration counter.
l1c=0

# Repeat the estimation procedure until the cluster assignments stabilize.
while (sum(Cl==Clee)<nl) {
  
  # Update the cluster labels.
  Cl=Clee
  
  # Standardize the feature matrix.
  Xl=as.matrix(scale(df))
  
  # Combine the data and updated cluster labels.
  D2=cbind(Xl,Cl)
  
  # Sort observations according to the updated cluster labels.
  D2=D2[order(D2[,p+1],decreasing=FALSE),]
  
  # Construct the cluster membership indicator matrix.
  pl=outer(D2[,p+1],1:KCe,FUN="==")
  dim(pl)=c(nl,KCe)
  
  # Calculate the number of observations in each cluster.
  nc=colSums(pl)
  
  # Calculate cumulative cluster sample sizes.
  Nc=rep(0,KCe+1)
  Nc[2:(KCe+1)]=cumsum(nc)
  
  
  # ------------------------------------------------
  # Re-estimate cluster means and covariance matrix
  # ------------------------------------------------
  
  # Initialize the cluster mean vector.
  Cmu=c()
  
  # Estimate the mean vector for each cluster.
  for (i in 1:KCe){
    Dc=D2[D2[,p+1]==i,1:p]
    Cmu=c(Cmu,cols2(Dc)/nc[i])
  }
  
  # Initialize the pooled covariance matrix.
  Csig=matrix(0,p,p)
  
  # Calculate the pooled covariance matrix.
  for(i in 1:KCe){
    Dc=D2[D2[,p+1]==i,1:p]
    mc=cols2(Dc)/nc[i]
    Csig=Csig+mo(t2(Dc)-mc)
  }
  
  # Normalize the pooled covariance matrix.
  Csig=Csig/nl
  
  
  # ------------------------------------------------
  # Recalculate transformed distances
  # ------------------------------------------------
  
  # Initialize the transformed distance vector.
  Y=c()
  
  # Calculate transformed distances for each cluster.
  for(i in 1:KCe){ 
    Dc=D2[D2[,p+1]==i,1:p]
    ylc=colSums((solve(Csig)%*%(t2(Dc)-Cmu[((i-1)*p+1):(i*p)]))*(t2(Dc)-Cmu[((i-1)*p+1):(i*p)]))
    Y=c(Y,(ylc^(p/2)+1)^(2/p)-1)
  }
  
  
  # Adjust the bandwidth for the next iteration.
  h_hat0=h_hat0*nl^(3/80)
  
  
  # ------------------------------------------------
  # Update the density generator
  # ------------------------------------------------
  
  # Define the updated kernel estimator of the density generator.
  g_hat0=function(y){
    K_4=K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="-")/h_hat0)/h_hat0+K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="+")/h_hat0)/h_hat0
    out=rowSums(K_4)*(y^(p/2)+1)^(2/p-1)*gamma(p/2)/(nl*pi^(p/2))
    out
  }
  
  
  # ------------------------------------------------
  # Initialize parameters for pseudo-likelihood optimization
  # ------------------------------------------------
  
  # Calculate the Cholesky decomposition of the updated covariance matrix.
  cho=chol(Csig)
  
  # Construct the updated initial parameter vector.
  ini=c(Cmu,cho[upper.tri(cho,diag=TRUE)],nc/nl)
  
  
  # ------------------------------------------------
  # Groupwise pseudo-likelihood objective function
  # ------------------------------------------------
  
  # Define the objective function for covariance parameter estimation.
  Ls=function(th){
    
    # Update the density generator using the current bandwidth.
    g_hat0=function(y){
      K_4=K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="-")/h_hat0)/h_hat0+K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="+")/h_hat0)/h_hat0
      out=rowSums(K_4)*(y^(p/2)+1)^(2/p-1)*gamma(p/2)/(nl*pi^(p/2))
      out
    }
    
    # Copy the current parameter vector.
    inir=ini
    
    # Replace the covariance parameters with the candidate values.
    inir[(p*KCe+2):(p*KCe+p*(p+1)/2)]=th
    
    # Calculate transformed distances for all clusters.
    for(i in 1:KCe){
      mu=inir[((i-1)*p+1):(i*p)]
      siv=inir[(p*KCe+1):(p*KCe+p*(p+1)/2)]
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      U=t(solve(V))
      Dc=D2[D2[,p+1]==i,1:p]
      Y[(Nc[i]+1):(Nc[i+1])]=(cols3((U%*%(t2(Dc)-mu))^2)^(p/2)+1)^(2/p)-1
    }
    
    # Initialize the groupwise log pseudo-likelihood.
    lfyu=0
    
    # Calculate the log pseudo-likelihood for the first KCe-1 clusters.
    for(i in 1:(KCe-1)){
      mu=ini[((i-1)*p+1):(i*p)]
      siv=inir[(p*KCe+1):(p*KCe+p*(p+1)/2)]
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      U=t(solve(V))
      Dc=D2[D2[,p+1]==i,1:p]
      Yl=cols3((U%*%(t2(Dc)-mu))^2)
      lf=sum(log(g_hat0(Yl)))+nc[i]*(sum(0.5*log(diag(U)^2)))
      lfyu=lfyu+lf
    }
    
    # Calculate the contribution from the final cluster.
    mu=ini[((KCe-1)*p+1):(KCe*p)]
    siv=inir[(p*KCe+1):(p*KCe+p*(p+1)/2)]
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    U=t(solve(V))
    Dc=D2[D2[,p+1]==KCe,1:p]
    Yl=cols3((U%*%(t2(Dc)-mu))^2)
    lf=sum(log(g_hat0(Yl)))+nc[KCe]*(sum(0.5*log(diag(U)^2)))
    lfyu=lfyu+lf
    
    # Return the negative average log pseudo-likelihood.
    -lfyu/nl
  }
  
  
  # ------------------------------------------------
  # Objective function for updating the first cluster mean
  # ------------------------------------------------
  
  L1=function(th){
    
    # Define the density generator for the current bandwidth.
    g_hat0=function(y){
      K_4=K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="-")/h_hat0)/h_hat0+K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="+")/h_hat0)/h_hat0
      out=rowSums(K_4)*(y^(p/2)+1)^(2/p-1)*gamma(p/2)/(nl*pi^(p/2))
      out
    }
    
    # Copy the current parameter vector and update the first mean.
    inir=ini
    inir[c(1:p)]=th
    
    vvv=ini[p*KCe+p*(p+1)/2+KCe+1]
    
    # Recalculate transformed distances under the candidate parameters.
    for(i in 1:KCe){
      mu=inir[((i-1)*p+1):(i*p)]
      siv=inir[(p*KCe+1):(p*KCe+p*(p+1)/2)]
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      U=t(solve(V))
      Dc=D2[D2[,p+1]==i,1:p]
      Y[(Nc[i]+1):(Nc[i+1])]=(cols3((U%*%(t2(Dc)-mu))^2)^(p/2)+1)^(2/p)-1
    }
    
    # Initialize the log pseudo-likelihood.
    lfyu=0
    
    # Calculate the contribution from the first cluster.
    mu=th[1:p]
    siv=ini[(p*KCe+1):(p*KCe+p*(p+1)/2)]
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    U=t(solve(V))
    Dc=D2[D2[,p+1]==1,1:p]
    Yl=cols3((U%*%(t2(Dc)-mu))^2)
    lf=sum(log(g_hat0(Yl)))+nc[1]*(sum(0.5*log(diag(U)^2)))
    lfyu=lfyu+lf
    
    # Calculate the contribution from the final cluster.
    mu=ini[((KCe-1)*p+1):(KCe*p)]
    siv=ini[(p*KCe+1):(p*KCe+p*(p+1)/2)]
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    U=t(solve(V))
    Dc=D2[D2[,p+1]==KCe,1:p]
    Yl=cols3((U%*%(t2(Dc)-mu))^2)
    lf=sum(log(g_hat0(Yl)))+nc[KCe]*(sum(0.5*log(diag(U)^2)))
    lfyu=lfyu+lf
    
    # Return the negative average log pseudo-likelihood.
    -lfyu/nl
  }
  
  
  # ------------------------------------------------
  # Objective function for updating the second cluster mean
  # ------------------------------------------------
  
  L2=function(th){
    
    # Define the density generator for the current bandwidth.
    g_hat0=function(y){
      K_4=K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="-")/h_hat0)/h_hat0+K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="+")/h_hat0)/h_hat0
      out=rowSums(K_4)*(y^(p/2)+1)^(2/p-1)*gamma(p/2)/(nl*pi^(p/2))
      out
    }
    
    # Copy the current parameter vector and update the second mean.
    inir=ini
    inir[(p+1):(2*p)]=th
    
    # Recalculate transformed distances under the candidate parameters.
    for(i in 1:KCe){
      mu=inir[((i-1)*p+1):(i*p)]
      siv=inir[(p*KCe+1):(p*KCe+p*(p+1)/2)]
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      U=t(solve(V))
      Dc=D2[D2[,p+1]==i,1:p]
      Y[(Nc[i]+1):(Nc[i+1])]=(cols3((U%*%(t2(Dc)-mu))^2)^(p/2)+1)^(2/p)-1
    }
    
    # Initialize the log pseudo-likelihood.
    lfyu=0
    
    # Calculate the contribution from the second cluster.
    mu=th[1:p]
    siv=ini[(p*KCe+1):(p*KCe+p*(p+1)/2)]
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    U=t(solve(V))
    Dc=D2[D2[,p+1]==2,1:p]
    Yl=cols3((U%*%(t2(Dc)-mu))^2)
    lf=sum(log(g_hat0(Yl)))+nc[2]*(sum(0.5*log(diag(U)^2)))
    lfyu=lfyu+lf
    
    # Return the negative average log pseudo-likelihood.
    -lfyu/nl
  }
  
  
  # ------------------------------------------------
  # Objective function for updating the third cluster mean
  # ------------------------------------------------
  
  L3=function(th){
    
    # Define the density generator for the current bandwidth.
    g_hat0=function(y){
      K_4=K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="-")/h_hat0)/h_hat0+K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="+")/h_hat0)/h_hat0
      out=rowSums(K_4)*(y^(p/2)+1)^(2/p-1)*gamma(p/2)/(nl*pi^(p/2))
      out
    }
    
    # Copy the current parameter vector and update the third mean.
    inir=ini
    inir[(2*p+1):(3*p)]=th
    
    # Recalculate transformed distances under the candidate parameters.
    for(i in 1:KCe){
      mu=inir[((i-1)*p+1):(i*p)]
      siv=inir[(p*KCe+1):(p*KCe+p*(p+1)/2)]
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      U=t(solve(V))
      Dc=D2[D2[,p+1]==i,1:p]
      Y[(Nc[i]+1):(Nc[i+1])]=(cols3((U%*%(t2(Dc)-mu))^2)^(p/2)+1)^(2/p)-1
    }
    
    # Initialize the log pseudo-likelihood.
    lfyu=0
    
    # Calculate the contribution from the third cluster.
    mu=th[1:p]
    siv=ini[(p*KCe+1):(p*KCe+p*(p+1)/2)]
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    U=t(solve(V))
    Dc=D2[D2[,p+1]==3,1:p]
    Yl=cols3((U%*%(t2(Dc)-mu))^2)
    lf=sum(log(g_hat0(Yl)))+nc[3]*(sum(0.5*log(diag(U)^2)))
    lfyu=lfyu+lf
    
    # Return the negative average log pseudo-likelihood.
    -lfyu/nl
  }
  
  
  # ------------------------------------------------
  # Objective function for updating the fourth cluster mean
  # ------------------------------------------------
  
  L4=function(th){
    
    # Define the density generator for the current bandwidth.
    g_hat0=function(y){
      K_4=K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="-")/h_hat0)/h_hat0+K(outer((y^(p/2)+1)^(2/p)-1,Y,FUN="+")/h_hat0)/h_hat0
      out=rowSums(K_4)*(y^(p/2)+1)^(2/p-1)*gamma(p/2)/(nl*pi^(p/2))
      out
    }
    
    # Copy the current parameter vector and update the fourth mean.
    inir=ini
    inir[(3*p+1):(4*p)]=th
    
    # Recalculate transformed distances under the candidate parameters.
    for(i in 1:KCe){
      mu=inir[((i-1)*p+1):(i*p)]
      siv=inir[(p*KCe+1):(p*KCe+p*(p+1)/2)]
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      U=t(solve(V))
      Dc=D2[D2[,p+1]==i,1:p]
      Y[(Nc[i]+1):(Nc[i+1])]=(cols3((U%*%(t2(Dc)-mu))^2)^(p/2)+1)^(2/p)-1
    }
    
    # Initialize the log pseudo-likelihood.
    lfyu=0
    
    # Calculate the contribution from the fourth cluster.
    mu=th[1:p]
    siv=ini[(p*KCe+1):(p*KCe+p*(p+1)/2)]
    V=matrix(0,p,p)
    V[upper.tri(V,diag=TRUE)]=siv
    U=t(solve(V))
    Dc=D2[D2[,p+1]==4,1:p]
    Yl=cols3((U%*%(t2(Dc)-mu))^2)
    lf=sum(log(g_hat0(Yl)))+nc[4]*(sum(0.5*log(diag(U)^2)))
    lfyu=lfyu+lf
    
    # Return the negative average log pseudo-likelihood.
    -lfyu/nl
  }
  
  
  # ------------------------------------------------
  # Iterative optimization of model parameters
  # ------------------------------------------------
  
  # Initialize the previous and current objective values.
  lfo=0
  lfn=10
  l4c=0
  
  # Continue optimization until the objective value converges.
  while(abs(lfo-lfn)>0.001){
    
    # Store the current objective value.
    lfo=lfn
    
    # Update the covariance parameters by minimizing the groupwise objective.
    ini[(p*KCe+2):(p*KCe+p*(p+1)/2)]=nlminb(ini[(p*KCe+2):(p*KCe+p*(p+1)/2)],Ls,control = list(trace=1,rel.tol=1e-3,x.tol=5e-4,iter.max=7))$par
    
    # Update the cluster proportions.
    ini[p*KCe+p*(p+1)/2+KCe]=1-sum(ini[(p*KCe+p*(p+1)/2+1):(p*KCe+p*(p+1)/2+KCe-1)])
    
    # Update the first cluster mean.
    ini[c(1:p)]=nlminb(ini[c(1:p)],L1,control = list(trace=1,rel.tol=1e-3,x.tol=5e-4,iter.max=7))$par
    
    # Update the cluster proportions.
    ini[p*KCe+p*(p+1)/2+KCe]=1-sum(ini[(p*KCe+p*(p+1)/2+1):(p*KCe+p*(p+1)/2+KCe-1)])
    
    # Update the second cluster mean.
    ini[c((p+1):(2*p))]=nlminb(ini[c((p+1):(2*p))],L2,control = list(trace=1,rel.tol=1e-3,x.tol=5e-4,iter.max=7))$par
    
    # Update the cluster proportions.
    ini[p*KCe+p*(p+1)/2+KCe]=1-sum(ini[(p*KCe+p*(p+1)/2+1):(p*KCe+p*(p+1)/2+KCe-1)])
    
    # Update the third cluster mean.
    ini[c((2*p+1):(3*p))]=nlminb(ini[c((2*p+1):(3*p))],L3,control = list(trace=1,rel.tol=1e-3,x.tol=5e-4,iter.max=7))$par
    
    # Update the cluster proportions.
    ini[p*KCe+p*(p+1)/2+KCe]=1-sum(ini[(p*KCe+p*(p+1)/2+1):(p*KCe+p*(p+1)/2+KCe-1)])
    
    # Update the fourth cluster mean.
    ini[(3*p+1):(4*p)]=nlminb(ini[(3*p+1):(4*p)],L4,control = list(trace=1,rel.tol=1e-3,x.tol=5e-4,iter.max=7))$par
    
    # Update the cluster proportions.
    ini[p*KCe+p*(p+1)/2+KCe]=1-sum(ini[(p*KCe+p*(p+1)/2+1):(p*KCe+p*(p+1)/2+KCe-1)])
    
    # Store the updated parameter vector.
    inir=ini
    
    # Recalculate transformed distances under the updated parameters.
    for(i in 1:KCe){
      mu=inir[((i-1)*p+1):(i*p)]
      siv=inir[(p*KCe+1):(p*KCe+p*(p+1)/2)]
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      U=t(solve(V))
      Dc=D2[D2[,p+1]==i,1:p]
      Y[(Nc[i]+1):(Nc[i+1])]=(cols3((U%*%(t2(Dc)-mu))^2)^(p/2)+1)^(2/p)-1
    }
    
    # Initialize the objective value and temporary vectors.
    lfn=0
    lfnl=0
    Yl=Yl2=c()
    
    # Calculate transformed distances for each cluster.
    for (i in 1:KCe){
      mu=ini[((i-1)*p+1):(i*p)]
      siv=ini[(p*KCe+1):(p*KCe+p*(p+1)/2)]
      V=matrix(0,p,p)
      V[upper.tri(V,diag=TRUE)]=siv
      U=t(solve(V))
      Dc=D2[D2[,p+1]==i,1:p]
      Yl=c(Yl,cols3((U%*%(t2(Dc)-mu))^2))
      Yl2=c(Yl2,(cols3((U%*%(t2(Dc)-mu))^2)^(p/2)+1)^(2/p)-1)
    }
    
    # Calculate the kernel matrix for the transformed distances.
    K_1=K(outer(Yl2,Yl2,FUN="-")/h_hat0)/h_hat0+K(outer(Yl2,Yl2,FUN="+")/h_hat0)/h_hat0
    dim(K_1)=c(length(Yl2),length(Yl2))
    
    # Remove diagonal elements from the kernel matrix.
    diag(K_1)=0
    
    # Calculate the density-related objective contribution.
    out=rowSums(K_1)*(Yl^(p/2)+1)^(2/p-1)*gamma(p/2)/(nl*pi^(p/2))
    
    # Update the objective value.
    lfn=lfn+sum(log(out[out!=0]))+sum(out!=0)*(sum(0.5*log(diag(U)^2)))+sum(out==0)*(log(min(out[out!=0]))+sum(0.5*log(diag(U)^2)))
    
    # Count the number of nonzero density estimates.
    lfnl=lfnl+sum(out!=0)
    
    # Normalize the objective value.
    lfn=-lfn/nl
    
    # Increase the inner iteration counter.
    l4c=l4c+1
    
    # Stop the inner optimization after 15 iterations.
    if(l4c>=15) break
  }
  
  
  # ------------------------------------------------
  # Update bandwidth and calculate pseudo-likelihood
  # ------------------------------------------------
  
  # Reverse the bandwidth adjustment.
  h_hat0=h_hat0*nl^(-3/80)
  
  # Store the final objective value from the inner optimization.
  lfni=lfn
  
  # Recalculate the covariance matrix used in the pseudo-likelihood.
  V2=matrix(0,p,p)
  V2[upper.tri(V2,diag=TRUE)]=ini[(p*KCe+1):(p*KCe+p*(p+1)/2)]
  pssi=t(V2)%*%V2
  
  # Calculate the determinant-related constant.
  detp=det(pssi)^(-0.5)
  
  # Initialize the cluster mean vector and transformed distance vector.
  pmu=Y2=c()
  
  # Extract the updated cluster means.
  for(i in 1:KCe){
    pmu=c(pmu,ini[((i-1)*p+1):(i*p)])
  }
  
  # Calculate transformed distances for each cluster.
  for(i in 1:KCe){ 
    Dc=D2[D2[,p+1]==i,1:p]
    ylc=colSums((solve(pssi)%*%(t2(Dc)-pmu[(p*(i-1)+1):(p*i)])*(t2(Dc)-pmu[(p*(i-1)+1):(p*i)])))
    Y2=c(Y2,(ylc^(p/2)+1)^(2/p)-1)
  }
  
  
  # ------------------------------------------------
  # Update cluster assignments
  # ------------------------------------------------
  
  # Calculate the distances of observations from each cluster.
  ri3=0
  Yt=c()
  
  for(i in 1:KCe){  # Calculate the distance for each cluster.
    Yt=cbind(Yt,cols3((solve(pssi)%*%(t(Xl)-pmu[(p*(i-1)+1):(p*i)]))*(t(Xl)-pmu[(p*(i-1)+1):(p*i)])))
  }
  
  # Initialize the pseudo-likelihood matrix.
  pp=matrix(0,nd,KCe)
  
  # Calculate the pseudo-likelihood for each observation and cluster.
  for (i in 1:KCe){
    pp[,i]=detp*g_hat0(Yt[,i])*ini[KCe*p+p*(p+1)/2+i]
  }
  
  # Initialize the updated cluster assignments.
  Clee=c()
  
  # Reassign each labeled observation to the most likely eligible cluster.
  for(i in 1:nl){
    cc=(D[i]==1)*KCe/2
    a=which(pp[i,]==max(pp[i,(cc+1):(cc+KCe/2)]))
    
    # Randomly select a cluster if multiple clusters share the maximum value.
    if(length(a)>1){
      b=sample(a,1,replace=FALSE,prob=rep(1,length(a)))
    }
    
    # Use the unique maximum if only one cluster is selected.
    if(length(a)==1){
      b=a
    }
    
    Clee=c(Clee,b)
  } 
  
  # Increase the outer iteration counter.
  l1c=l1c+1
  
  # Stop the outer iteration after 20 iterations.
  if(l1c>=20) break
}


# --------------------------------------------------
# Final calculation and output
# --------------------------------------------------

# Construct the sequence of integration points.
r_seq=seq(0,((max(Y)+1)^(p/2)-1)^(2/p),length=5000)

# Calculate the final integral-based quantity.
A=sum(r_seq^(p+1)*g_hat0(r_seq^2)*2*pi^(p/2)/gamma(p/2)*((max(Y)+1)^(p/2)-1)^(2/p)/5000/p)

# Save the complete R workspace in the current working directory.
save.image("dependent_2.RData")