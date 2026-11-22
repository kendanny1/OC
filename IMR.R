library(expm)
library(openxlsx)
p=7
tr=c()
for(KC in 2:6){
  result=c()
  
  n=388
  load("/Users/kendanny/Downloads/Optimal_Classification/application/pima_age.RData")
  n=388
  Cle=c()
  Clek=c()
  rii=rii2=c()
  X=as.matrix(scale(df))
  D=diabetes_data[,9]
  X2=X
  X2[D==0,]=t(t(X2[D==0,])-colSums(X2[D==0,])/sum(D==0))
  X2[D==1,]=t(t(X2[D==1,])-colSums(X2[D==1,])/sum(D==1))
  kx_km=kmeans(X2,KC,algorithm ="Lloyd",iter.max=100)
  cx_km=kx_km$cluster
  WSS_kmeans <- kx_km$tot.withinss
  dx=dist(X2,method="euclidean")
  hx=hclust(dx,method="ward.D2")
  cx_hc=cutree(hx,k=KC)
  # 計算 hierarchical clustering 的群內平方和 WSS
  WSS_hc <- 0
  for (i in 1:KC) {
    Xi <- X2[cx_hc == i, , drop = FALSE]
    mu_hc <- colMeans(Xi)
    WSS_hc <- WSS_hc + sum(rowSums((Xi - mu_hc)^2))
  }
  if (WSS_kmeans < WSS_hc) {
    cx=cx_km
  }else{
    cx=cx_hc
  }
  Cle=c(Cle,cx)
  # ---- 在 kmeans 之後，建立初始 bgamma 與 Csig（替換原本該段）
  Csig <- matrix(0, p, p)
  bgamma <- matrix(0, 2*p, KC)
  
  for (d in 0:1) {
    for (k in 1:KC) {
      idx <- which(cx == k & D == d)
      if (length(idx) == 0) {
        bgamma[(d*p+1):((d+1)*p), k] <- NA
      } else if (length(idx) == 1) {
        bgamma[(d*p+1):((d+1)*p), k] <- X[idx, ]
      } else {
        bgamma[(d*p+1):((d+1)*p), k] <- colMeans(X[idx, , drop = FALSE])
      }
    }
  }
  
  # 用全 label 平均填補 NA
  for (d in 0:1) {
    idx_label <- which(D == d)
    global_mu <- colMeans(X[idx_label, , drop = FALSE])
    for (k in 1:KC) {
      mu_slot <- bgamma[(d*p+1):((d+1)*p), k]
      if (any(is.na(mu_slot))) bgamma[(d*p+1):((d+1)*p), k] <- global_mu
    }
  }
  
  # 計算初始 Csig
  Csig <- matrix(0, p, p)
  for (d in 0:1) {
    for (k in 1:KC) {
      idx <- which(cx == k & D == d)
      if (length(idx) == 0) next
      mu_kd <- bgamma[(d*p+1):((d+1)*p), k]
      diffs <- sweep(X[idx, , drop = FALSE], 2, mu_kd)
      Csig <- Csig + t(diffs) %*% diffs
    }
  }
  Csig <- Csig / n
  
  # ---- EM-like 更新迴圈（替換原本 while 區塊）
  cx2 <- cx
  cx <- rep(0, n)
  max_iter <- 100
  iter <- 0
  eps_reg <- 1e-8
  
  while (!all(cx == cx2) && iter < max_iter) {
    iter <- iter + 1
    cx <- cx2
    
    # invert with regularization
    invCsig <- solve(Csig)
    
    # 計算 Yu (n x KC) 向量化
    Yu <- matrix(0, n, KC)
    for (k in 1:KC) {
      mu0 <- bgamma[1:p, k]; mu1 <- bgamma[(p+1):(2*p), k]
      Mu <- matrix(0, n, p)
      if (sum(D == 0) > 0) Mu[D == 0, ] <- matrix(mu0, sum(D==0), p, byrow = TRUE)
      if (sum(D == 1) > 0) Mu[D == 1, ] <- matrix(mu1, sum(D==1), p, byrow = TRUE)
      Diff <- X - Mu
      tmp <- Diff %*% invCsig
      Yu[, k] <- rowSums(tmp * Diff)
    }
    
    # 更新 cx2（向量化）
    cx2 <- apply(Yu, 1, which.min)
    
    # 重新計算 bgamma 與 Csig（同上）
    bgamma <- matrix(0, 2*p, KC)
    for (d in 0:1) {
      for (k in 1:KC) {
        idx <- which(cx2 == k & D == d)
        if (length(idx) == 0) {
          bgamma[(d*p+1):((d+1)*p), k] <- NA
        } else if (length(idx) == 1) {
          bgamma[(d*p+1):((d+1)*p), k] <- X[idx, ]
        } else {
          bgamma[(d*p+1):((d+1)*p), k] <- colMeans(X[idx, , drop = FALSE])
        }
      }
    }
    
    # 用 label 全體平均補 NA
    for (d in 0:1) {
      idx_label <- which(D == d)
      global_mu <- colMeans(X[idx_label, , drop = FALSE])
      for (k in 1:KC) {
        if (any(is.na(bgamma[(d*p+1):((d+1)*p), k])))
          bgamma[(d*p+1):((d+1)*p), k] <- global_mu
      }
    }
    
    # 計算 Csig
    Csig <- matrix(0, p, p)
    for (d in 0:1) {
      for (k in 1:KC) {
        idx <- which(cx2 == k & D == d)
        if (length(idx) == 0) next
        mu_kd <- bgamma[(d*p+1):((d+1)*p), k]
        diffs <- sweep(X[idx, , drop = FALSE], 2, mu_kd)
        Csig <- Csig + t(diffs) %*% diffs
      }
    }
    Csig <- Csig / n
  } # end while
  
  if (iter == max_iter) warning("達到最大迭代次數，可能未收斂")
  
  tr=c(tr,sum(diag(Csig)))
  Clek=c(Clek,cx2)
  save.image(paste("/Users/kendanny/Downloads/Optimal_Classification/application/km_",KC,".RData",sep=""),version = 2)
}
Csig=(t(X)-colSums(X)/n)%*%t(t(X)-colSums(X)/n)/n
tr=c(sum(diag(Csig)),tr)


