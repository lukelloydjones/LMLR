#!/bin/bash
# File path to X matrix
# File path to Y matrix
# Lasso constraint parameter lower bound 20 good for 500, 1000
# Lasso constraint parameter upper bound 20 good for 500, 1000
# Convergence criterion for MM algorithm within each run
# Maximim iterations for MM algorithm within each run
# Perturbation parameter for screening zero effects
# Path to write output
# Path to starting beta parameter matrix
# Set of starting sigma parameters. 2 mixture components in this example
# Set of starting mean  parameters. 2 mixture components in this example
# Set of starting pi    parameters. 2 mixture components in this example
# Path to tee out standard output. Adjust "out/baseball_3g_mreg.log" comp only
../../lmlr sim1_geno.csv \
           sim1_pheno.csv \
		   1 \
		   80 \
           1e-4 \
           500000 \
           1e-16 \
           out/sim_ \
           betas_str_1.txt \
           1.5 2.2\
           -10 10\
           0.48 0.52 \
           2>&1 | tee out/sim.log
