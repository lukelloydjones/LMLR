# ===========================================================
# Authors : Luke Lloyd-Jones and Hien Nguyen 
# Data started : 10/05/2015
# Date last updated: 11/07/2017
# Make 50 replicates of 2 component test data
# ===========================================================
# -----------------------------------------------------------
# Generate the data for sim 3 - scenario 1
# -----------------------------------------------------------
rm(list = ls())
setwd("")
source('fmr_data_gen_csv4.R', chdir = F)
# ---------------------------------------
# Simulation 1
# ---------------------------------------
ext      <- ""
sim.name <- "sim"
for (i in seq(1, 1))
{
  sim.name.i <- paste0(ext, sim.name, i)
  n = 200
  p = 20
  sd = 1
  mu <- c(-20, 20)
  g  <- 2
  pis <- c(0.4, 0.6)
  # Beta vector
  beta.1 <- c(2, -1, 2)
  beta.2 <- c(-2, 2, -1)
  beta.tru   <- matrix(0, nrow = p, ncol = g)
  beta.tru[1:length(beta.1), 1]  <- beta.1
  beta.tru[(length(beta.1) + 2):(length(beta.1) + 1 + length(beta.2)), 2]  <- beta.2
  # Simulate the name
  FMRDataGen4Csv(n, p, "", sim.name.i, sd, mu, g, pis, beta.tru)
  beta <- matrix(0, nrow = p, ncol = g)
  beta <- beta + 0.1
  betas.out <- paste0("betas_str_", i, ".txt")
  write.table(beta, betas.out, col.names = F, row.names = F, sep = ",")
}
