# ===============================================================
# Try and simulate some data for mixture of regression 
# ===============================================================
# Set the number of individuals to be n and p
n <- 100
p <- 500
# Create an observed Z matrix drawn from N(0, 1)
X <- matrix(rnorm(n * p), n, p)
# Add in the intercept column 
X <- cbind(array(1, n), X)
# Create two sets of betas each with 10 non-zero elements
beta.1 <- c(rep(1, 5), rep(0, 495))
beta.2 <- c(rep(0, 495), rep(1, 5))
# Add the intercepts for each group
beta.1 <- c(5, beta.1)
beta.2 <- c(-5, beta.2)
# Check the dimensions
dim(X)
length(beta.1)
length(beta.2)
# Cycle over the rows and sample a value from each mean
mu.1 <- array(0, n)
mu.2 <- array(0, n)
for (i in seq(1, n)) {	
  mu.1[i] <- rnorm(1, sum(X[i, ] * beta.1), 1)
  mu.2[i] <- rnorm(1, sum(X[i, ] * beta.2), 2)
}
# Create the phenotype
Y <- c(sample(mu.1, 40), sample(mu.2, 60))
plot(density(Y), col = "green", lwd = 2.5)
# Run our lasso on it and get the results
#x <- CDLasso(X[, -1], Y, 15)
# Plot and compare
#lines(density(X%*%x), col = "blue", lwd = 2.5)
# Lasso seems to work but chooses the wrong set of betas
# ===============================================================
# Try and fix mixture of regressions code
# ===============================================================
# ------------------------------------------------------------------------------
# Initialise the parameters - start with k = 2 
# ------------------------------------------------------------------------------
beta.1 <- matrix(0, nrow = dim(X)[2] , ncol = 1)
beta.2 <- matrix(0, nrow = dim(X)[2] , ncol = 1)
beta.1[1:6] <- c(5, 1, 1, 1, 1, 1)
beta.2[c(1, 497:501)] <- c(-5, 1, 1, 1, 1, 1)
# Remember that the first element is the mean element of each group
pi.1   <- 0.4
pi.2   <- 1 - pi.1
var.1  <- 1
var.2  <- 2
eps    <- 10e-5 # The perturbation factor from Hunter and Li
lambda <- 0.5
# ------------------------------------------------------------------------------
# Initialise the weights
# ------------------------------------------------------------------------------
mu.1  <- X %*% beta.1
mu.2  <- X %*% beta.2
w.1   <- ((pi.1)  * dnorm(Y, mu.1, sqrt(var.1))) / 
         (((pi.1) * dnorm(Y, mu.1, sqrt(var.1))) + 
         ((pi.2)  * dnorm(Y, mu.2, sqrt(var.2))))
w.2   <- ((pi.2)  * dnorm(Y, mu.2, sqrt(var.2))) / 
         (((pi.1) * dnorm(Y, mu.1, sqrt(var.1))) + 
         ((pi.2)  * dnorm(Y, mu.2, sqrt(var.2))))
# ==============================================================================
# Begin the EM Algorithm 
# ==============================================================================
for (r in seq(1, 500)) {
	# ------------------------------------------------------------------------------
	# Update the the weights given the pis and variances
	# ------------------------------------------------------------------------------
	mu.1  <- X %*% beta.1
	mu.2  <- X %*% beta.2
	w.1 <- array(0, n)
	w.2 <- array(0, n)
	for (i in seq(1, n))
	{
	 w.1[i]   <- (pi.1  * dnorm(Y[i], beta.1[1], sqrt(var.1)) ) / 
	          (pi.1  * dnorm(Y[i], beta.1[1], sqrt(var.1)) + 
	           pi.2  * dnorm(Y[i], beta.2[1], sqrt(var.2)) )
	 w.2[i]   <- (pi.2  * dnorm(Y[i], beta.2[1], sqrt(var.2)) ) / 
	          (pi.1  * dnorm(Y[i], beta.1[1], sqrt(var.1)) + 
	           pi.2  * dnorm(Y[i], beta.2[1], sqrt(var.2)) )
	}
	# ------------------------------------------------------------------------------
	# Update the pis and variances 
	# ------------------------------------------------------------------------------
	#pi.1  <- max(sum(w.1) / n, 10e-6)
	#pi.2  <- max(sum(w.2) / n, 10e-6)
	pi.1  <- 0.4
	pi.2  <- 0.6
	# var.1 <- max(sum(w.1 * (Y - X %*% beta.1) ^ 2) / sum(w.1), 10e-6)
	# var.2 <- max(sum(w.2 * (Y - X %*% beta.2) ^ 2) / sum(w.2), 10e-6)
	# ------------------------------------------------------------------------------
	# Generate the penalisation matrices
	# ------------------------------------------------------------------------------
	# Generate the penalty matrices with the perturbed denominator
	#beta.1.pen <- unlist(lapply(beta.1, p.scad.deriv, lambda = 1, a = 3.7, n = n))
	#beta.2.pen <- unlist(lapply(beta.2, p.scad.deriv, lambda = 1, a = 3.7, n = n))
	beta.1.pen <- pi.1 * sqrt(n) * lambda
	beta.2.pen <- pi.2 * sqrt(n) * lambda
	diag.pen.1 <- c(beta.1.pen / (eps + abs(beta.1[-1])))
	diag.pen.2 <- c(beta.2.pen / (eps + abs(beta.2[-1])))
	E.mat.1    <- diag(diag.pen.1)
	E.mat.2    <- diag(diag.pen.2)
	# ------------------------------------------------------------------------------
	# Update beta 1 - Generate the weighted matrices - do the update 
	# ------------------------------------------------------------------------------
	W.mat.1 <- diag(w.1) * (1 / var.1)
	W.mat.2 <- diag(w.2) * (1 / var.2)
	beta.1  <- solve(t(X[, -1]) %*% W.mat.1 %*% X[, -1] + E.mat.1) %*% (t(X[, -1]) %*% W.mat.1 %*% Y)
	beta.2  <- solve(t(X[, -1]) %*% W.mat.2 %*% X[, -1] + E.mat.2) %*% (t(X[, -1]) %*% W.mat.2 %*% Y)
	beta.1 <- c(5,  beta.1)
	beta.2 <- c(-5, beta.2)
	# ------------------------------------------------------------------------------
	# Update beta 1
	# ------------------------------------------------------------------------------
	# # Build the denominator of the first group
	# X.sqr.pi      <- matrix(0, ncol = p, nrow = p)
	# for (i in seq(1, n)) {
	  # X.sqr.temp    <- X[i, -1] %*% t(X[i, -1])
	  # X.sqr.pi      <- X.sqr.pi + w.1[i] *  X.sqr.temp
	# }
	# beta.1.denom     <- pi.1 * W.mat.1 + (1 / var.1) * X.sqr.pi
	# beta.1.denom.inv <- solve(beta.1.denom)
	# #beta.1.denom.inv <- ginv(beta.1.denom)
	# # Build the numerator of the first group
	# X.y      <- matrix(0, ncol = p, nrow = 1)
	# for (i in seq(1, n)){
	  # X.y.temp <- X[i, -1] * Y[i]
	  # X.y      <- X.y + w.1[i] *  X.y.temp
	# }
	# beta.1.numer <- as.numeric((1 / var.1) * X.y) 
	# # Update beta.1
	# beta.1.m <- beta.1
	# beta.1   <- beta.1.denom.inv %*% beta.1.numer
	# # ------------------------------------------------------------------------------
	# # Update beta 2
	# # ------------------------------------------------------------------------------
	# # Build the denominator of the second group
	# X.sqr.pi      <- matrix(0, ncol = p, nrow = p)
	# for (i in seq(1, n)) {
	  # X.sqr.temp <- X[i, -1] %*% t(X[i, -1])
	  # X.sqr.pi   <- X.sqr.pi  + w.2[i] *  X.sqr.temp
	# }
	# beta.2.denom     <- pi.2 * W.mat.2 + (1 / var.2) * X.sqr.pi 
	# beta.2.denom.inv <- solve(beta.2.denom)
	# #beta.2.denom.inv <- ginv(beta.2.denom)
	# # Build the numerator of the first group
	# X.y      <- matrix(0, ncol = p, nrow = 1)
	# for (i in seq(1, n)){
	  # X.y.temp <- X[i, -1] * Y[i]
	  # X.y      <- X.y + w.2[i] *  X.y.temp
	# }
	# beta.2.numer <- as.numeric((1 / var.2) * X.y) 
	# # Update beta.2
	# beta.2.m <- beta.2
	# beta.2   <- beta.2.denom.inv %*% beta.2.numer
	# # ------------------------------------------------------------------------------
	# # Update the weights given the betas
	# # ------------------------------------------------------------------------------
	# w.1   <- ((pi.1)  * dnorm(Y, mu.1, sqrt(var.1))) / 
	         # (((pi.1) * dnorm(Y, mu.1, sqrt(var.1))) + 
	         # ((pi.2)  * dnorm(Y, mu.2, sqrt(var.2))))
	# w.2   <- ((pi.2)  * dnorm(Y, mu.2, sqrt(var.2))) / 
	         # (((pi.1) * dnorm(Y, mu.1, sqrt(var.1))) + 
	         # ((pi.2)  * dnorm(Y, mu.2, sqrt(var.2))))
	#beta.1 <- c(5, beta.1)
	#beta.2 <- c(-5, beta.2)
	# ------------------------------------------------------------------------------
	# Go back to step 1 
	# ------------------------------------------------------------------------------
	#plot(density(X %*% beta.1), lwd = 2,  col = "blue")
	#lines(density(X %*% beta.2), lwd = 2, col = "pink")
	#lines(density(pi.1 * (X %*% beta.1) + pi.2 * (X %*% beta.2)), lwd = 2, col = "yellow")
	#lines(density(Y), lwd = 2, col = "green")
	#plot(density(Y), lwd = 2, col = "green")
	#print(c(pi.1, pi.2, var.1, var.2))
}