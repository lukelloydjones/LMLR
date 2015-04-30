# ==============================================================================
# R script to prototype the finite mixture of regression model for GWAS in 
# admixed populations
# ==============================================================================
rm(list = ls())
setwd("~/Dropbox/Post_Doc_QBI/Fin_Mix_Reg/")
source('~/Dropbox/Post_Doc_QBI/Fin_Mix_Reg/R_Code/scad_scad_deriv_func.R')
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
X  <- X[, 1:100]
P  <- dim(X)[2]
# ------------------------------------------------------------------------------
# Correct geno for allele frequencies - puts each column on mu = 0, var = 1
# ------------------------------------------------------------------------------
for (j in seq(1, P)) {
 q <- sum(X[, j] / (2 * n))			    #Calculate the reference allele freq
 X[, j] <- (X[, j] - 2 * q) / sqrt(2 * q * (1 - q))	#Centre and then scale 
}
# ------------------------------------------------------------------------------
# Initialise the parameters - start with k = 2 
# ------------------------------------------------------------------------------
beta.1 <- matrix(1, nrow = dim(X)[2], ncol = 1)
beta.2 <- matrix(1, nrow = dim(X)[2], ncol = 1)
pi.1   <- 0.4
pi.2   <- 1 - pi.1
var.1  <- 1.5
var.2  <- 1.6
eps    <- 10e-5 # The perturbation factor from Hunter and Li
# ------------------------------------------------------------------------------
# Initialise the weights
# ------------------------------------------------------------------------------
mu.1  <- X %*% beta.1[, 1]
mu.2  <- X %*% beta.2[, 1]
w.1   <- ((pi.1)  * dnorm(Y, mu.1, sqrt(var.1))) / 
         (((pi.1) * dnorm(Y, mu.1, sqrt(var.1))) + 
         ((pi.2)  * dnorm(Y, mu.2, sqrt(var.2))))
w.2   <- ((pi.2)  * dnorm(Y, mu.2, sqrt(var.2))) / 
         (((pi.1) * dnorm(Y, mu.1, sqrt(var.1))) + 
         ((pi.2)  * dnorm(Y, mu.2, sqrt(var.2))))
# ==============================================================================
# Begin the EM Algorithm 
# ==============================================================================
for (r in seq(1, 5000)) {
# ------------------------------------------------------------------------------
# Update the pis and variances 
# ------------------------------------------------------------------------------
pi.1  <- max(sum(w.1) / n, 10e-6)
pi.2  <- max(sum(w.2) / n, 10e-6)
var.1 <- max(sum(w.1 * (Y - X %*% beta.1[, 1]) ^ 2) / sum(w.1), 10e-6)
var.2 <- max(sum(w.2 * (Y - X %*% beta.2[, 1]) ^ 2) / sum(w.2), 10e-6)
# ------------------------------------------------------------------------------
# Update the the weights given the pis and variances
# ------------------------------------------------------------------------------
mu.1  <- X %*% beta.1[, 1]
mu.2  <- X %*% beta.2[, 1]
w.1   <- ((pi.1)  * dnorm(Y, mu.1, sqrt(var.1))) / 
         (((pi.1) * dnorm(Y, mu.1, sqrt(var.1))) + 
         ((pi.2)  * dnorm(Y, mu.2, sqrt(var.2))))
w.2   <- ((pi.2)  * dnorm(Y, mu.2, sqrt(var.2))) / 
         (((pi.1) * dnorm(Y, mu.1, sqrt(var.1))) + 
         ((pi.2)  * dnorm(Y, mu.2, sqrt(var.2))))
# ------------------------------------------------------------------------------
# Update the betas given the pis, variances and weights
# ------------------------------------------------------------------------------
# Generate the penalty matrices with the perturbed denominator
beta.1.pen <- unlist(lapply(beta.1, p.scad.deriv, lambda = 1, a = 3.7, n = n))
beta.2.pen <- unlist(lapply(beta.2, p.scad.deriv, lambda = 1, a = 3.7, n = n))
diag.1     <- c(beta.1.pen / (eps + beta.1))
diag.2     <- c(beta.2.pen / (eps + beta.2))
W.mat.1    <- diag(diag.1)
W.mat.2    <- diag(diag.2)
# ------------------------------------------------------------------------------
# Update beta 1
# ------------------------------------------------------------------------------
# Build the denominator of the first group
X.sqr.pi      <- matrix(0, ncol = P, nrow = P)
for (i in seq(1, n)) {
  X.sqr.temp <- X[i, ] %*% t(X[i, ])
  X.sqr.pi      <- X.sqr.pi + w.1[i] *  X.sqr.temp
}
beta.1.denom     <- pi.1 * W.mat.1 + (1 / var.1) * X.sqr.pi
beta.1.denom.inv <- solve(beta.1.denom)
#beta.1.denom.inv <- ginv(beta.1.denom)
# Build the numerator of the first group
X.y      <- matrix(0, ncol = P, nrow = 1)
for (i in seq(1, n)){
  X.y.temp <- X[i, ] * Y[i]
  X.y      <- X.y + w.1[i] *  X.y.temp
}
beta.1.numer <- as.numeric((1 / var.1) * X.y) 
# Update beta.1
beta.1.m <- beta.1
beta.1   <- beta.1.denom.inv %*% beta.1.numer
# ------------------------------------------------------------------------------
# Update beta 2
# ------------------------------------------------------------------------------
# Build the denominator of the second group
X.sqr.pi      <- matrix(0, ncol = P, nrow = P)
for (i in seq(1, n)) {
  X.sqr.temp <- X[i, ] %*% t(X[i, ])
  X.sqr.pi   <- X.sqr.pi  + w.2[i] *  X.sqr.temp
}
beta.2.denom     <- pi.2 * W.mat.2 + (1 / var.2) * X.sqr.pi 
beta.2.denom.inv <- solve(beta.2.denom)
#beta.2.denom.inv <- ginv(beta.2.denom)
# Build the numerator of the first group
X.y      <- matrix(0, ncol = P, nrow = 1)
for (i in seq(1, n)){
  X.y.temp <- X[i, ] * Y[i]
  X.y      <- X.y + w.2[i] *  X.y.temp
}
beta.2.numer <- as.numeric((1 / var.2) * X.y) 
# Update beta.2
beta.2.m <- beta.2
beta.2   <- beta.2.denom.inv %*% beta.2.numer
# ------------------------------------------------------------------------------
# Update the weights given the betas
# ------------------------------------------------------------------------------
w.1   <- ((pi.1)  * dnorm(Y, mu.1, sqrt(var.1))) / 
         (((pi.1) * dnorm(Y, mu.1, sqrt(var.1))) + 
         ((pi.2)  * dnorm(Y, mu.2, sqrt(var.2))))
w.2   <- ((pi.2)  * dnorm(Y, mu.2, sqrt(var.2))) / 
         (((pi.1) * dnorm(Y, mu.1, sqrt(var.1))) + 
         ((pi.2)  * dnorm(Y, mu.2, sqrt(var.2))))
# ------------------------------------------------------------------------------
# Go back to step 1 
# ------------------------------------------------------------------------------
plot(density(X %*% beta.1), lwd = 2,  col = "blue")
lines(density(X %*% beta.2), lwd = 2, col = "pink")
lines(density(pi.1 * (X %*% beta.1) + pi.2 * (X %*% beta.2)), lwd = 2, col = "yellow")
lines(density(Y), lwd = 2, col = "green")
#plot(density(Y), lwd = 2, col = "green")
print(c(pi.1, pi.2, var.1, var.2))
}