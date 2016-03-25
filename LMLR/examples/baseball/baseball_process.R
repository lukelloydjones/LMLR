# ==============================================================================
# Script to do initial processing of Baseball data
# Author: Luke Lloyd-Jones
# Date started: 30/01/2016
# Data last updated 08/02/2016
# ==============================================================================
rm(list = ls())
library(ggplot2)
setwd("~/Dropbox/post_doc_qbi/Fin_Mix_Reg/baseball/")
bball <- read.csv("baseball_dat.csv", header = F)
head(bball)
Y <- log(bball[, 1])
hist(log(bball[, 1]), 13)
plot(density(Y))
# -----------------------------------------------------------------------------
# Do a nice histogram density plot of the salaries
# -----------------------------------------------------------------------------
m <- ggplot(bball, aes(x = log(bball[, 1])))
m + geom_histogram(aes(y = ..density..), binwidth = 0.3) + 
    geom_line(lwd = 2, adjust = 0.9, stat = "density") +
    xlab("Log of salary in 1000s") +
    ylab("Density") +
    theme(axis.title=element_text(size=16,face="bold"))
ggsave(filename = "salaries_1990s.pdf", 
       plot = last_plot(), 
       path = "~/Desktop", 
       scale = 1, width = 45, height = 20, 
       units = "cm", dpi = 300)
# -----------------------------------------------------------------------------
# Process the covariates
# -----------------------------------------------------------------------------
X <- bball[, 2:17]
dim(X)
X2 <- cbind(X[, 1] * X[, 13], X[, 1] * X[, 14], X[, 1] * X[, 15], X[, 1] * X[, 16],
            X[, 3] * X[, 13], X[, 3] * X[, 14], X[, 3] * X[, 15], X[, 3] * X[, 16],
            X[, 7] * X[, 13], X[, 7] * X[, 14], X[, 7] * X[, 15], X[, 7] * X[, 16],
            X[, 8] * X[, 13], X[, 8] * X[, 14], X[, 8] * X[, 15], X[, 8] * X[, 16])
X3 <- cbind(X, X2)
X4 <- scale(X3)
#write.table(Y,  "bball_resp.csv", col.names = F, row.names = F, quote = F, sep = ",")
#write.table(X4, "bball_pred.csv", col.names = F, row.names = F, quote = F, sep = ",")
# Write out the betas to take in as initial guesses - lasso set
beta.str.1  <- c(0, 0, 0, 0.2, 0, 0, -0.19, 0.26, 
                 0, 0, 0, 0, 0.79, 0.72, 0.15, 0,
             -0.21, 0.63, 0.34, 0, 0, 0.14, 0, -0.18,
                 0, 0, 0, 0, 0.29, -0.14, 0, 0) 
beta.str.2  <- c(-0.32, 0.29, -0.70, 0.96, 0, 0, 0, 0, 
                 0, 0, 0, 0, 0.70, 0, 0.50, -0.36,
                 0, 0, 0, 0, 0, -0.38, 0, 0.74,
                 0, 0, 0.34, 0, -0.46, 0, 0, 0) 
beta.str <- cbind(beta.str.1, beta.str.2)
# Write out the betas to take in as initial guesses - SCAD set
# beta.str.1  <- c(0, 0, 0, 0.4, 0, 0, 0, 0, 
                 # 0, 0, 0, 0, 1.49, -1.24, 0.68, 0,
                 # -0.93, 1.17, 0, 0, 0, 0, 0, 0,
                 # 0, 0, 0, 0, 0.39, 0, 0, 0) 
# beta.str.2  <- c(0.58, 0, 0, 0, 0, 0, 0, 0, 
                 # 0, 0, 0, 0, 0.95, 0, 0, 0,
                 # 0, 0, 0, 0, 0, 0, 0, 0,
                 # 0.71, 0, 0, 0, -1.12, 0, 0, 0) 
# beta.str.1  <- array(0.01, 32)
# beta.str.2  <- array(0.01, 32)
# beta.str <- cbind(beta.str.1, beta.str.2)
# write.table(beta.str, "betas_str.txt", col.names = F, 
            # row.names = F, quote = F, sep = ",")
# The predictions from Chen and Khalili
X = X4
alphas_mat   <- matrix(0, nrow = 1, ncol = 2)
pis_mat      <- matrix(0, nrow = 1, ncol = 2)
sigma_mat    <- matrix(0, nrow = 1, ncol = 2)
alphas_mat[1, ] <- c(6.41, 7.00)
pis_mat[1, ]    <- c(0.72, 0.28)
sigma_mat[1, ]  <- c(0.25, 0.25)
n = dim(X)[1]
g <- 2
tau <- matrix(0, nrow = g, ncol = n)
for (k in 1:n) 
{
inner <- 0
for (j in 1:g) 
{
	inner <- inner + 
	         pis_mat[, j] * 
	         dnorm(Y[k],
	         alphas_mat[, j] + sum(beta.str[, j] * X[k, ]),
	         sqrt(sigma_mat[, j]))
}
for (j in 1:g) 
{
	 tau[j, k] <- pis_mat[, j] * 
	              dnorm(Y[k], 
	              alphas_mat[, j] + sum(beta.str[, j] * X[k, ]),
	              sqrt(sigma_mat[, j])) / 
	              inner
}
}
X <- as.matrix(X)
c1.kc <- round(tau[1, ]) * (alphas_mat[, 1] + X%*%(beta.str[, 1]))
c2.kc <- round(tau[2, ]) * (alphas_mat[, 1] + X%*%(beta.str[, 1]))
Y.pred.kc <- tau[1, ] * (alphas_mat[, 1] + X%*%(beta.str[, 1])) + 
	         tau[2, ] * (alphas_mat[, 2] + X%*%(beta.str[, 2]))
# ==============================================================================
# Process the results from our model
# ==============================================================================
X = X4
alphas_mat   <- matrix(0, nrow = 1, ncol = 2)
pis_mat      <- matrix(0, nrow = 1, ncol = 2)
sigma_mat    <- matrix(0, nrow = 1, ncol = 2)
rsqrs        <- matrix(0, nrow = 1, ncol = 2)
setwd("~/Dropbox/Post_Doc_QBI/Fin_Mix_Reg/baseball/")
library("GenABEL")
# Read in the beta results for this file			    
beta.res   <- read.csv("beta_estimates.txt", header = F, fill = T)
# Read in the estimates of pi
pi.res     <- read.table("pi_estimates.txt", header = F, fill = T)
pis_mat    <- t(pi.res)
# Read in the estimates of means
mn.res     <- read.table("alpha_estimates.txt", header = F, fill = T)
alphas_mat <- t(mn.res)
# Read in the estimates of sigmas
sigma.res  <- read.table("sigma_estimates.txt", header = F, fill = T)
sigma_mat  <- t(sigma.res)
n = dim(X)[1]
g <- 2
tau <- matrix(0, nrow = g, ncol = n)
for (k in 1:n) 
{
inner <- 0
for (j in 1:g) 
{
	inner <- inner + 
	         pis_mat[, j] * 
	         dnorm(Y[k],
	         alphas_mat[, j] + sum(beta.res[, j] * X[k, ]),
	         sqrt(sigma_mat[, j]))
}
for (j in 1:g) 
{
	 tau[j, k] <- pis_mat[, j] * 
	              dnorm(Y[k], 
	              alphas_mat[, j] + sum(beta.res[, j] * X[k, ]),
	              sqrt(sigma_mat[, j])) / 
	              inner
}
}
names <-     c("(Intercept)",     			
     "AVG",		 		
     "OBP",          			
     "R",               			
     "H",               			
     "2B",             			
     "3B",            			
     "HR",            		
     "RBI",           			
     "BB" ,           			
     "SO"  ,          			
     "SB"   ,         			
     "ERS"	,				
     "FAE"	,			
     "FA"	,			
     "AE"	,			
     "ARB"	,			
     "AVG_FAE",	
     "AVG_FA",		
     "AVG_AE",		
     "AVG_ARB",	 
     "R_FAE",		
     "R_FA",		
     "R_AE",		
     "R_ARB",		
     "HR_FAE",		    
     "HR_FA"	,	
     "HR_AE",		
     "HR_ARB",		
     "RBI_FAE",		
     "RBI_FA"	,	
     "RBI_AE"	,	
     "RBI_ARB")
X <- as.matrix(X)
colnames(X) <- names[-1]
Y.pred <- tau[1, ] * (alphas_mat[, 1] + X%*%(beta.res[, 1])) + 
	      tau[2, ] * (alphas_mat[, 2] + X%*%(beta.res[, 2]))
p = dim(X)[2]
beta.marg <- matrix(0, nrow = p, ncol = 1)
dat <- as.data.frame(cbind(Y, X[, 1:dim(X)[2]]))
y.lm <- lm(Y ~ ., data = dat)
library(MASS)
step <- stepAIC(y.lm, k = log(n), direction = "both")
summary(step)
y.lm.pred <- step$fitted.values
# Do a nice plot of tau scores versus
# ----------------------------------- 
x.dat <- data.frame(cbind(Y, tau[1, ]))
colnames(x.dat) <- c("Y", "tau_1")
p <- ggplot(x.dat, aes(tau_1, Y))
p + geom_point() + xlab(expression(tau[1])) + theme(text = element_text(size=20)) + 
    ylab("Log of salary in 1000s") +
    theme(axis.title=element_text(size=16,face="bold"))
ggsave(filename = "tau_salaries_1990s.pdf", 
       plot = last_plot(), 
       path = "~/Desktop", 
       scale = 1, width = 45, height = 20, 
       units = "cm", dpi = 300)
# --------------------------------------------
# Bind up for table
# Make a column of zeros and a column of names
# --------------------------------------------
stp.res <- data.frame(names, 0)
colnames(stp.res) <- c("Cov", "Estimate_LM")
rownames(stp.res) <- stp.res[, 1]
stp.res.2 <- data.frame(names(summary(step)$coefficients[, 1]), 
                              summary(step)$coefficients[, 1])
colnames(stp.res.2) <- c("Cov", "Estimate_LM")
stp.res[rownames(stp.res.2), 2] <- stp.res.2[, 2]
# Our results
beta.res.2 <- data.frame(colnames(X), beta.res)
colnames(beta.res.2) <- c("Cov", "Estimate_C1", "Estimate_C2")
int <- data.frame("(Intercept)", alphas_mat)
colnames(int) <-  c("Cov", "Estimate_C1", "Estimate_C2")
beta.res.bnd <- rbind(int, beta.res.2)
# Chen and Khalili's results
alphas_mat.kc <- c(6.41, 7.00)
beta.str.2 <- data.frame(colnames(X), beta.str)
colnames(beta.str.2) <- c("Cov", "Estimate_C1_KC", "Estimate_C2_KC")
int <- data.frame("(Intercept)", alphas_mat.kc[1], alphas_mat.kc[2] )
colnames(int) <-  c("Cov", "Estimate_C1_KC", "Estimate_C2_KC")
beta.res.bnd.kc <- rbind(int, beta.str.2)
# Bind up the final table
all.bnd <- inner_join(stp.res, beta.res.bnd.kc, beta.res.bnd, by = "Cov")
all.bnd.2 <- inner_join(all.bnd, beta.res.bnd, by = "Cov")
library(knitr)
kable(cbind(all.bnd.2[, 1], round(all.bnd.2[, -1], 2)), format = "latex")
# -------------------------------------------------------------------------
# Plot the predicted and expected with the linear model and Chen's fit
# -------------------------------------------------------------------------
library(reshape)
library(ggplot2)
# Density plot
x <- cbind(Y, y.lm.pred, Y.pred.kc, Y.pred)
colnames(x) <- c("Observed", "Linear model step BIC", "MIXLASSO", "LMLR")
x.m <- melt(x)
colnames(x.m) <- c("X1", "Model", "X3")
m <- ggplot()
m + geom_line(data = x.m, aes(x = X3, group = Model, 
                 col = Model, linetype = Model), 
                 alpha = 1, stat="density",
                 binwidth = 0.3, lwd = 1.1) +   
    scale_linetype_manual(values = c(rep("solid", 4))) +
    xlab("Log of salary in 1000s") +
    xlim(c(3, 13)) +
    ylab("Density") +
    theme(axis.title=element_text(size = 16, face = "bold"), 
          legend.position = "none")
ggsave(filename = "pred_salaries_1990s.pdf", 
       plot = last_plot(), 
       path = "~/Desktop", 
       scale = 1, width = 45, height = 20, 
       units = "cm", dpi = 300)
# R^2 of prediction values
plot(Y, Y.pred, xlim = c(4, 13), ylim = c(4, 13))
points(Y, y.lm$fitted.values, col = "blue")
points(Y, Y.pred.kc, col = "green")
# Our model
r1      <- lm(Y.pred ~ Y)
mse.lmr <- sum((Y.pred - Y) ^ 2) / length(Y.pred)
mse.lmr
summary(r1)$adj.r.squared
# KC
r2 <- lm(y.lm.pred ~ Y)
mse.lm <- sum((y.lm.pred - Y) ^ 2) / length(Y.pred)
mse.lm
summary(r2)$adj.r.square
# Linear model
r3 <- lm(Y.pred.kc ~ Y)
mse.kc <- sum((Y.pred.kc - Y) ^ 2) / length(Y.pred)
mse.kc
summary(r3)$adj.r.square


