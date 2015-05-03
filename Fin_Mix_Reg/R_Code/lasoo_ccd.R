# ==============================================================================
# R script to prototype cyclic coordinate descent for lasso penalty
# Algorithm is outline on Tong Tong and Lange 2008
# ==============================================================================
rm(list = ls())
setwd("~/Dropbox/Post_Doc_QBI/Fin_Mix_Reg/")
source('~/Dropbox/Post_Doc_QBI/Fin_Mix_Reg/R_Code/scad_scad_deriv_func.R')
# -------------------------------------------------------------------------------
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
# Initialise the betas to be all 0 as discussed by Tong Tong and Lange
# ------------------------------------------------------------------------------
beta       <- matrix(0, nrow = dim(X)[2], ncol = 1)
bwd.dds    <- matrix(0, nrow = dim(X)[2], ncol = 1)
fwd.dds    <- matrix(0, nrow = dim(X)[2], ncol = 1)
dds        <- matrix(0, nrow = dim(X)[2], ncol = 2)
dg.dbeta   <- matrix(0, nrow = dim(X)[2], ncol = 1)
lambda     <- 1
mu         <- 0
for (i in seq(1, 100))
{
# ------------------------------------------------------------------------------
# For each beta we calculate the forward and backward directional derivative
# ------------------------------------------------------------------------------
for (k in seq(1, length(beta)))
{
  # Calculate the partial derivative of the sum of squares
  dg.dbeta[k]  <- -sum(((Y - mu - X%*%beta) * X[, k]))
  # Find the backward directional derivatives
  if (beta[k] > 0) 
  {
	bwd.dds[k] <- -dg.dbeta[k] - lambda
  } else 
  {
    bwd.dds[k] <- -dg.dbeta[k] + lambda
  }
  # Find the forward directional derivatives
  if (beta[k] >= 0) 
  {
	fwd.dds[k] <- dg.dbeta[k] + lambda
  } else 
  {
    fwd.dds[k] <- dg.dbeta[k] - lambda
  }
  dds[k, ] <- c(bwd.dds[k], fwd.dds[k])
}
# ------------------------------------------------------------------------------
# If both directional derivatives are non-negative then we don't update the
# beta value for that cycle. If either is negative then we solve for the min
# in that direction. Due to convexity it should be impossible for both DDs to
# be negative.
# ------------------------------------------------------------------------------
beta.old <- beta
up.bwd <- which((dds[, 1] < 0))
up.fwd <- which((dds[, 2] < 0))
print(length(up.bwd))
print(length(up.fwd))
stopifnot(intersect(up.bwd, up.fwd) == 0)
for (k in up.bwd)
{
  up.pot.bwd  <- beta[k] - ((dg.dbeta[k] - lambda) / sum(X[, k]^2))
  beta[k]     <- min(0, up.pot.bwd)
}
for (k in up.fwd)
{
  up.pot.fwd  <- beta[k] - ((dg.dbeta[k] + lambda) / sum(X[, k]^2))
  beta[k]     <- max(0, up.pot.fwd)
}
# ------------------------------------------------------------------------------
# Update the intercept
# ------------------------------------------------------------------------------
mu <- (1 / n) * sum(Y - X%*%beta)
# ------------------------------------------------------------------------------
# Report the summed difference in betas
# ------------------------------------------------------------------------------
#sum((beta - beta.old)^2)
print(sum((Y - mu - X%*%beta)^2))
plot(density(X %*% beta), lwd = 2,  col = "blue")
lines(density(Y), lwd = 2, col = "green")
}