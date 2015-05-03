# ==============================================================================
# Experimenting with coordinate descent library provided by Tong Tong and Lange
# ==============================================================================
set.seed(1001)
rm(list = ls())
install.packages("CDLasso")
library("CDLasso")
# ------------------------------------------------------------------------------
# Simulate the data with more predictors than individuals
# ------------------------------------------------------------------------------
n <- 500
p <- 2000
nz <- c(1:5)
true.beta     <- rep(0, p)
true.beta[nz] <- c(1, 1, 1, 1, 1)
X <- matrix(rnorm(n * p), n, p)
Y <- X %*% true.beta
rownames(X) <- 1:nrow(X)
colnames(X) <- 1:ncol(X)
# ------------------------------------------------------------------------------
# Perform the L2 regression. Requires the predictors to be the rows
# ------------------------------------------------------------------------------
outL2     <- l2.reg(X, Y, 2)
outL2est  <- l2.reg(X[outL2$selected, ], Y, lambda = 0)
crossval2 <- cv.l2.reg(X, Y, 10, seq(0, 1, 0.01))
plot(crossval2)
out2      <- l2.reg(X, Y, crossval2$lam.opt)
# ------------------------------------------------------------------------------
# Let's try it on our data
# ------------------------------------------------------------------------------
setwd("~/Dropbox/Post_Doc_QBI/Fin_Mix_Reg/")
# ------------------------------------------------------------------------------
# Read in the data phenotype and genotypes
# ------------------------------------------------------------------------------
X     <- read.table("Data/cape_verde_r10k.txt")
pheno <- read.table("Data/cape_verde_pheno.txt", header = T)
rownames(pheno) <- pheno[, 1]
Y <- pheno[, 3] # Let's try for skin first as there are less mising values
X <- X[-which(Y == -9), ]
Y <- Y[-which(Y == -9)]
n <- length(Y)
library(MASS)
# ------------------------------------------------------------------------------
# Get rid of the NAs for now. Deal with them properly later 
# ------------------------------------------------------------------------------
X[is.na(X)] <- sample( c(0, 1, 2), length(which(is.na(X))), replace = TRUE)
X  <- as.matrix(X)
X2 <- X
X  <- X[, 1:5000]
P  <- dim(X)[2]
# ------------------------------------------------------------------------------
# Correct geno for allele frequencies - puts each column on mu = 0, var = 1
# ------------------------------------------------------------------------------
for (j in seq(1, P)) {
 q <- sum(X[, j] / (2 * n))			    #Calculate the reference allele freq
 X[, j] <- (X[, j] - 2 * q) / sqrt(2 * q * (1 - q))	#Centre and then scale 
}
# ------------------------------------------------------------------------------
# Perform the L2 regression. Requires the predictors to be the rows
# ------------------------------------------------------------------------------
crossval2 <- cv.l2.reg(t(X), Y, 10, seq(10, 40, 1))
plot(crossval2)
out2      <- l2.reg(t(X), Y, 21)
out2$selected

