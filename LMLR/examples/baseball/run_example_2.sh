#!/bin/bash
# File path to X matrix (requires full path)
# File path to Y matrix (requires full path)
# Lasso constraint parameter lower bound 20 good for 500, 1000
# Lasso constraint parameter upper bound 20 good for 500, 1000
# Convergence criterion
# Maximim iterations
# Perturbation paramter
# First sigma parameter
# Second sigma parameter
# First mu parameter
# Second mu parameter
# First pi parameter
# Second pi parameter
# Path to write output (requires full path). Also the path to the starting betas
# The starting betas must be named "betas_str.txt"
../../lmlr /Users/hayleywise/Dropbox/Git_Repos/LMLR/examples/baseball/bball_pred.csv \
           /Users/hayleywise/Dropbox/Git_Repos/LMLR/examples/baseball/bball_resp.csv \
           15 \
           60 \
           1e-4 \
           500000 \
           1e-16 \
           1.2 \
           1.1 \
           6.5 \
           6.6 \
           0.6 \
           0.4 \
           /Users/hayleywise/Dropbox/Git_Repos/LMLR/examples/baseball/ \
           2>&1 | tee /Users/hayleywise/Dropbox/Git_Repos/LMLR/examples/baseball/baseball.log
