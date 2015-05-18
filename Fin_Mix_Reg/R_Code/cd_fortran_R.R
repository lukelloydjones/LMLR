# ==============================================================================
# Attempts to translate Lange's CD code from Fortran
# ==============================================================================
# ------------------------------------------------------------------------------
# Simulated data with many more predictors than individuals
# ------------------------------------------------------------------------------
n <- 100
p <- 500
nz <- c(1:5)
true.beta     <- rep(0, p)
true.beta[nz] <- array(1, length(nz))
X <- matrix(rnorm(n * p), n, p)
Y <- X %*% true.beta
rownames(X) <- 1:nrow(X)
colnames(X) <- 1:ncol(X)
# ------------------------------------------------------------------------------
# Function LASSO_PENALISED_L2_REGRESSION translated
# ------------------------------------------------------------------------------
# Uses module Coordinate Descent from above
# Defines variables and memory
# Integers
# Doubles
CDLasso <- function(X_mat, Y_mat, lambda)
{
CRITERION = 10e-5
EPSILON   = 10e-8
A <- 0.0 
B <- 0.0 
C <- 0.0 
DL2 <- 0.0
LAMBDA <- lambda
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
Y <- Y_mat
M <- length(Y)
# Add a columns of ones fro the intercept
X <- cbind(array(1, M), X_mat)
N <- dim(X)[2]
# Array to store the beta estimates
ESTIMATE <- array(0,   N)
# Residual vectors for left and right derivatives
SUM_X_SQUARES <- matrix(0, nrow =  N, ncol = 1)
LEFT_R        <- matrix(, nrow = M, ncol = 1)
R             <- matrix(, nrow = M, ncol = 1)
RIGHT_R       <- matrix(, nrow = M, ncol = 1)
# Initialise the number of cases M and the number of predictors N
# ------------------------------------------------------------------------------
# Initialise the residual vector and the penalty
# ------------------------------------------------------------------------------
R <- Y
if (abs(ESTIMATE[1] > 0)) { R = R - ESTIMATE[1]} # Adjust for the mean
PENALTY = 0
for (i in seq(2, N))
{
	A = ESTIMATE[i]
	B = abs(A)
	if (B > 0) {
		R = R - A * X[, i]     # Update the residual for each estimate
		PENALTY = PENALTY + B  # Update the penalty
	}
}
# ------------------------------------------------------------------------------
# Initialise the objective funtion and penalty. This is for the check on the
# objective function. 
# ------------------------------------------------------------------------------
L2 <- sum(R ^ 2)             # Sum of the residuals
PENALTY = LAMBDA * PENALTY   # Lambda times the sum of the penalty
OBJECTIVE = L2 / 2 + PENALTY # Objective is the sum of squares plut the penalty
# ------------------------------------------------------------------------------
# Start the main loop. When we read in X the first column must be a vector of 1s
# ------------------------------------------------------------------------------
for (ITERATION in seq(1, 1000))
{
  # Update the intercept
  # --------------------
  A = ESTIMATE[1]
  ESTIMATE[1] = A + sum(R) / M  # There is a double negative accounted for here
  R = R + A - ESTIMATE[1]       # Residual update of the mean
  
  # Update the other regression coefficients
  # ----------------------------------------
  
  for (i in seq(2, N))
  {
  	#i = 2
  	DL2 = -sum(R * X[, i])
  	A = ESTIMATE[i]
  	B = abs(A)
  	# Go to the next update if the directional derivatives are both positive
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
  	# Update the residuals to the right
  	for (j in seq(1, M))
  	{
        RIGHT_R[j] <- R[j] + C * X[j, i]
        RIGHT_L2   <- RIGHT_L2 + RIGHT_R[j] ^ 2
  	}
  	RIGHT_PENALTY   = PENALTY + LAMBDA * (RIGHT_ROOT - B)
  	RIGHT_OBJECTIVE = RIGHT_L2 / 2 + RIGHT_PENALTY
  	# Find the root to the left of 0
  	# -------------------------------
  	LEFT_ROOT = min(A - (DL2 - LAMBDA) / SUM_X_SQUARES[i], 0)
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
  	stop("*** ERROR *** OBJECTIVE FUNCTION INCREASE")
  	break
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

return(ESTIMATE)
} 
# ------------------------------------------------------------------------------
# Cross validation
# ------------------------------------------------------------------------------  
CrossVal <- function(X, Y, lambda, k)
{
	# Split X and Y over the folds
	n <- dim(X)[1]
	no.each.fold <- n / k
	mse <- array(0, k)
	for (i in seq(0, k - 1))
	{
	  k.elem <- seq(i * no.each.fold + 1, (i + 1) * no.each.fold)
	  # The training sets
	  X.train <- X[-k.elem, ]
	  Y.train <- Y[-k.elem]
	  # Calculate the lasso parameters from the training
	  est.k <- cd_lasso(X.train, Y.train, lambda)
	  # The prediction sets
	  X.pred  <- X[k.elem, ]
	  dim(X.pred)
	  length(est.k[-1])
	  Y.est   <- est.k[1] + X.pred%*%est.k[-1]
	  Y.pred  <- Y[k.elem]
	  mse[i + 1] <- sum(abs(Y.est - Y.pred) ) / length(Y.pred)
	}
	return(mean(mse))
}

lassoParam  <- function(lam)
{	x <- cross_val(X, Y, lam, 10) 
	print(lam)
	print(x)
	x 
} 

 
optimise(lassoParam, c(0, 20))
  







