# ==============================================================================
# Try and simulate some data for mixture of regression 
# ==============================================================================
# b1 <- rnorm(300, 0.5)
# b2 <- rnorm(200,  -0.5)
# beta.1 <- rnorm(500, 1)
# beta.2 <- rnorm(500, -1)
beta.1.t <- c(rep(0.5,10), rep(0,90))
beta.2.t <- c(rep(0,90), rep(0.5,10))

#plot(density(beta.1), lwd = 2, col = "blue")
#lines(density(beta.2), lwd = 2, col = "green")
#samp.spc <- seq(1, 500)
#sm.1 <- sample(samp.spc, 300)
#sm.2 <- sample(samp.spc, 200)
#beta.1 <- matrix(0, nrow = dim(X)[2], ncol = 1)
#beta.2 <- matrix(0, nrow = dim(X)[2], ncol = 1)
#beta.1[sm.1] <- b1
#beta.2[sm.2] <- b2
mu.1 <- array(0, 684)
mu.2 <- array(0, 684)
for (i in seq(1, 684)) {
	
  mu.1[i] <- rnorm(1, sum(X[i, ] * beta.1.t), 1)
  mu.2[i] <- rnorm(1, sum(X[i, ] * beta.2.t), 1)
}
#plot(density(mu.1), lwd = 2, col = "blue")
#lines(density(mu.2), lwd = 2, col = "green")
Y <- c(sample(mu.1, 284), sample(mu.2, 400))
hist(Y)

