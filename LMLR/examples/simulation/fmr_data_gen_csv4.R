# ===========================================================
# Author: Luke Lloyd-Jones
# Start date: 09/09/15 
# Date last updated: 05/02/15
# Code to generate SNP data and phenotypes in plink format
# ===========================================================
	
FMRDataGen4Csv <- function(n, p, ext, sim.name, 
                           sd, mu, g, pis, beta) {
  # Args:
  #  n: number of individuals 
  #  p: number of SNPs or covariates 
  #  ext: place and extension to write plink files to 
  #  sim.name: the name of the simulation for reference
  #  no.not.zero: number of non zero SNPs
  #  eff.size: the effect size of the non zero SNPs
  #  sd: standard deviation 
  #  mu.1: mean of the first mixture
  #  mu.2: meand of the second mixture
  # 
  # Returns:
  # Writes plink files of size n and p to ext file
  #
  # Specify a path
  path <- ext
  # Simulate a SNP like matrix 
  snp.mat  <- matrix(NA, ncol = p, nrow = n)
  gen.freq <- runif(p, 0.05, 0.5)
  # Generate the 0, 1, 2 from binomial with the potential 
  # successes  to be 2
  for(i in 1:p) 
  {
    q <- gen.freq[i]
    snp.mat[,i] <- rbinom(n, 2, q)
  }
  for (j in seq(1, p)) 
  {
    q <- sum(snp.mat[, j] / (2 * n))			    
    snp.mat[, j] <- (snp.mat[, j] - 2 * q) / 
                     sqrt(2 * q * (1 - q))	
  }
  X <- snp.mat
  # Set the means and the empty outcome vector
  Y <- array(0, n)
  # ---------------------
  # Generate Y new school
  # ---------------------
  n.grp <- pis * n
  a <- cumsum(n.grp)
  b <- cumsum(n.grp) + 1
  inds <- rbind(c(1, b[-length(b)]), a)
  for (j in seq(1, g))
  {
  	for (i in seq(inds[1, j], inds[2, j]))
  	{
  	  print(i)
  	  Y[i] <- mu[j] + sum(beta[, j] * X[i, ])  + rnorm(1, 0, sd)	
  	}
  }
  # print(n.grp)
  print(beta)
  print(plot(density(Y)))
  # # ---------------------
  # # Generate Y old school
  # # ---------------------
  # for (i in seq(1, (n * (5 / 10)))) 
  # {
    # Y[i] <- mu.1 + 
            # sum(beta[, 1] * X[i, ])  #+ 
            # rnorm(1, 0, sd)
  # }
  # for (i in seq((n * (5 / 10) + 1), n)) {
    # Y[i] <- mu.2 + 
            # sum(beta[, 2] * X[i, ]) #+ 
            # rnorm(1, 0, sd)
  # }
  write.table(X, paste0(ext, sim.name, "_geno.csv"),  
              row.names = F, col.names = F,
              sep = ",",
              quote = F)
  write.table(Y, paste0(ext, sim.name, "_pheno.csv"),
              row.names = F, col.names = F,
              sep = ",", 
              quote = F)
}

