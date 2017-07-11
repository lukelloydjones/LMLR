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
# Set of starting sigma parameters. 3 mixture components in this example
# Set of starting mean  parameters. 3 mixture components in this example
# Set of starting pi    parameters. 3 mixture components in this example
# Path to tee out standard output. Adjust "out/baseball_3g_mreg.log" comp only
../../lmlr bball_pred.csv \
           bball_resp.csv \
		   1 \
		   100 \
           1e-4 \
           500000 \
           1e-16 \
           out/baseball_3g_mreg_ \
           kac_3g_betas.txt \
           1.5 2.2 2.2\
           6.0 6.5 7.0\
           0.4 0.2 0.4 \
           2>&1 | tee out/baseball_3g_mreg.log
