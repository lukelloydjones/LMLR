#!/bin/bash
# File path to X matrix
# File path to Y matrix
# Lasso constraint parameter lower bound 20 good for 500, 1000
# Lasso constraint parameter upper bound 20 good for 500, 1000
# Convergence criterion
# Maximim iterations
# Perturbation paramter
# Path to write output
# Set of starting sigma parameters
# Set of starting mean  parameters
# Set of starting pi    parameters
# Bath to starting matrix of betas for initialising
./fmr_lasso       /Users/l.lloydjones/Dropbox/Post_Doc_QBI/Fin_Mix_Reg/cpp_6/test/sim_1_1_geno.csv  \
                  /Users/l.lloydjones/Dropbox/Post_Doc_QBI/Fin_Mix_Reg/cpp_6/test/sim_1_1_pheno.csv \
                  1 \
                  100 \
                  1e-4 \
                  500000 \
                  1e-16 \
                  /Users/l.lloydjones/Dropbox/Post_Doc_QBI/Fin_Mix_Reg/cpp_6/test/out/ \
                  /Users/l.lloydjones/Dropbox/Post_Doc_QBI/Fin_Mix_Reg/cpp_6/test/out/betas_str_1.txt \
                  1.5 2.2\
                  -15 15\
                  0.32 0.68  \
                  2>&1 | tee /Users/l.lloydjones/Dropbox/Post_Doc_QBI/Fin_Mix_Reg/cpp_6/test/sim_test_brent.log