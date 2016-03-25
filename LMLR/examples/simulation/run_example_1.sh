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
../../lmlr /Users/hayleywise/Dropbox/Git_Repos/LMLR/examples/simulation/sim_x.csv \
           /Users/hayleywise/Dropbox/Git_Repos/LMLR/examples/simulation/sim_y.csv \
           10 \
           60 \
           1e-4 \
           500000 \
           1e-16 \
           2.2 \
           2.1 \
           -10 \
           10 \
           0.6 \
           0.4 \
           /Users/hayleywise/Dropbox/Git_Repos/LMLR/examples/simulation/ \
           2>&1 | tee /Users/hayleywise/Dropbox/Git_Repos/LMLR/examples/simulation/simulation.log
