# Replication Code

This repository contains the computational code used in the paper:

**Teng, J.-C., Fan, S.-H., Chiang, C.-T., and Lim, A. (2026). Optimal Classification with Latent Clusters**

The R scripts in this repository implement the estimation, clustering, and cluster-selection procedures described in the paper. The role of each script is summarized below.

## Code Description

### 1. `IMR.R`

This script implements the procedure described in **Section A.2.1** of the paper. It combines residual-based clustering with the proposed iterative Mahalanobis reassignment to obtain initial estimates of the model parameters.

### 2. `SP.R`

This script implements the algorithm described in **Section A.2.2** of the paper.

For each specified number of clusters, the estimates of $\beta$ and $\gamma$ obtained from `IMR.R` are used as the initial values for the algorithm. The script then performs the proposed estimation procedure and obtains the corresponding parameter estimates and objective function values.

### 3. `SP_2.R`

This script compares the results obtained from `SP.R` across different numbers of clusters and selects the estimates of the model parameters corresponding to the minimum objective function value.

### 4. `PML.R`

This script implements the **pseudo-maximum likelihood (PML) estimation**.


## Summary of the Computational Procedure

The overall computational workflow can be summarized as follows:

1. `IMR.R` obtains initial estimates using the procedure in Section A.2.1.
2. `SP.R` uses the estimated $\beta$ and $\gamma$ as initial values and implements the algorithm in Section A.2.2.
3. `SP_2.R` compares the results across different numbers of clusters and selects the estimates of the model parameters based on the minimum objective function value.
4. `PML.R` implements the pseudo-maximum likelihood estimation procedure.
