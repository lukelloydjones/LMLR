# ==============================================================================
# Attempts to translate Lange's CD code from Fortran
# ==============================================================================
# ------------------------------------------------------------------------------
# Coordinate descent module
# ------------------------------------------------------------------------------
n <- 500
p <- 2000
nz <- c(1:5)
true.beta     <- rep(0, p)
true.beta[nz] <- c(1, 1, 1, 1, 1)
X <- matrix(rnorm(n * p), n, p)
Y <- X %*% true.beta
rownames(X) <- 1:nrow(X)
colnames(X) <- 1:ncol(X)
# ------------------------------------------------------------------------------
# End - Coordinate descent module
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Subroutine L1GREEDY - Skip this one
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Subroutine LASSO_PENALISED_L2_REGRESSION
# ------------------------------------------------------------------------------
# Function of 
#LassoPenalisedRegressionL2 <- function(X_INPUT, Y, LAMBDA, CASES,PREDICTORS, 
#                                    L2, R,
 #                                   OBJECTIVE, PENALTY, ESTIMATE)
  # X_in - real matrix with nrows = predictors, ncols = cases, but then transposed
  # Y - real vector of size Cases
  # lambda - double
  # Cases - int
  # Predictors - int
  # L2
  # R
  # Objective
  # Penalty
  # Estimate 
# Uses module Coordinate Descent from above
# Defines variables and memory
# Integers
# Doubles
CRITERION = 10e-5
EPSILON   = 10e-8
A <- 0.0 
B <- 0.0 
C <- 0.0 
DL2 <- 0.0
LAMBDA <- 2
L2 <- 0.0 
OBJECTIVE <- 0.0
PENALTY <- 0.0
NEW_OBJECTIVE <- 0.0
LEFT_L2 <- 0.0
LEFT_OBJECTIVE <- 0.0
LEFT_PENALTY <- 0.0
LEFT_ROOT <- 0.0
RIGHT_L2 <- 0.0
RIGHT_OBJECTIVE <- 0.0 
RIGHT_PENALTY <- 0.0
RIGHT_ROOT <- 0.0
M <- length(Y)
X <- cbind(array(1, M), X)
N <- dim(X)[2]
ESTIMATE <- array(0,   N)
#X_INPUT  <- matrix(,  nrow = M, ncol = N )
SUM_X_SQUARES <-  matrix(0, nrow =  N, ncol = 1)
LEFT_R  <- matrix(, nrow = M, ncol = 1)
R       <- matrix(, nrow = M, ncol = 1)
RIGHT_R <- matrix(, nrow = M, ncol = 1)
# Initialise the number of cases M and the number of predictors N
# ------------------------------------------------------------------------------
# Initialise the residual vector and the penalty
# ------------------------------------------------------------------------------
R <- Y
if (abs(ESTIMATE[1] > 0)) { R = R - ESTIMATE[1]}
PENALTY = 0
for (i in seq(2, N))
{
	A = ESTIMATE[i]
	B = abs(A)
	if (B > 0) {
		R = R - A * X[, i]
		PENALTY = PENALTY + B
	}
}
# ------------------------------------------------------------------------------
# Initialise the objective funtion and penalty. This is for the check on the
# objective function
# ------------------------------------------------------------------------------
L2 <- sum(R ^ 2)
PENALTY = LAMBDA * PENALTY
OBJECTIVE = L2 / 2 + PENALTY
# ------------------------------------------------------------------------------
# Start the main loop. When we read in X the first column must be a vecotr of 1s
# ------------------------------------------------------------------------------
for (ITERATION in seq(1, 1000))
{
  # Update the intercept
  # --------------------
  A = ESTIMATE[1]
  ESTIMATE[1] = A + sum(R) / M
  R = R + A - ESTIMATE[1]
  
  # Update the other regression coefficients
  # ----------------------------------------
  
  for (i in seq(2, N))
  {
  	#i = 1
  	DL2 = -sum(R * X[, i])
  	A = ESTIMATE[i]
  	B = abs(A)
  	if (B < EPSILON)
  	{
  	  if (DL2 + LAMBDA >= 0 & -DL2 + LAMBDA >= 0)
  	  {
  	    next	  
  	  }
  	}
  	# Find the root to the right of 0
  	# -------------------------------
  	if (SUM_X_SQUARES[i] <= 0) {SUM_X_SQUARES[i] = sum(X[, i] ^ 2)}
  	RIGHT_ROOT = max(A - (DL2 + LAMBDA) / SUM_X_SQUARES[i], 0)
  	RIGHT_L2   = 0.0
  	C = A - RIGHT_ROOT
  	for (j in seq(1, M))
  	{
        RIGHT_R[j] <- R[j] + C * X[j, i]
        RIGHT_L2   <- RIGHT_L2 + RIGHT_R[j] ^ 2
  	}
  	RIGHT_PENALTY   = PENALTY + LAMBDA * (RIGHT_ROOT - B)
  	RIGHT_OBJECTIVE = RIGHT_L2 / 2 + RIGHT_PENALTY
  	# Find the root to the left of 0
  	# -------------------------------
  	LEFT_ROOT = min(A - (DL2 - LAMBDA)/SUM_X_SQUARES[i], 0)
  	LEFT_L2   = 0.0
  	C = A - LEFT_ROOT
  	for (j in seq(1, M))
  	{
        LEFT_R[j] <- R[j] + C * X[j, i]
        LEFT_L2   <- LEFT_L2 + LEFT_R[j] ^ 2
  	}
  	LEFT_PENALTY   = PENALTY + LAMBDA * (abs(LEFT_ROOT) - B)
  	LEFT_OBJECTIVE = LEFT_L2 / 2 + LEFT_PENALTY
  	# Choose between the two roots
  	# ----------------------------
  	if (RIGHT_OBJECTIVE <= LEFT_OBJECTIVE)
  	{
  	  R           = RIGHT_R
  	  ESTIMATE[i] = RIGHT_ROOT
  	  L2          = RIGHT_L2
  	  PENALTY     = RIGHT_PENALTY
  	} else
  	{
  	  R           = LEFT_R
  	  ESTIMATE[i] = LEFT_ROOT
  	  L2          = LEFT_L2
  	  PENALTY     = LEFT_PENALTY	
  	}
  }
  
  NEW_OBJECTIVE = L2 / 2 + PENALTY
  
  # Check for descent failure or convergence. If neither occurs,
  # record the new value of the objective function
  
  if (NEW_OBJECTIVE > OBJECTIVE)
  {
  	print("*** ERROR *** OBJECTIVE FUNCTION INCREASE")
  	stopifnot(NEW_OBJECTIVE <= OBJECTIVE)
  }
  if (OBJECTIVE - NEW_OBJECTIVE < CRITERION)
  {
  	print("***We Have convergence***")
  	break
  } else 
  {
  	OBJECTIVE = NEW_OBJECTIVE
  	print(NEW_OBJECTIVE)
  }
  
}
ESTIMATE  
  
  
  
  
  
  
  
  







