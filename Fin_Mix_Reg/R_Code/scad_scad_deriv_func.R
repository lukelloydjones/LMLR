# =====================================================================
# Functions for SCAD penalty function and its derivative
# =====================================================================
# ---------------------------------------------------------------------
# SCAD penalty function
# ---------------------------------------------------------------------
p.scad <- function(beta, lambda, a, n) {
  # Calculates the SCAD penalty - corresponds to a quadratic spline
  # with knots at lambda and a * lambda
  # Args:
  #  beta   - the regression parameter
  #  lambda - > 0 parameter
  #  a      - > 2 parameter   
  #  n      - number of individuals
  #
  # Returns:
  #  The resultant value of the scad penalty
  
	if (sqrt(n) * abs(beta) <= lambda)
	{
		p <- lambda * sqrt(n) * abs(beta)
		
	} else if (lambda < sqrt(n) * abs(beta) & 
	           sqrt(n) * abs(beta) <= a * lambda)
	{
		p <- - sqrt(n) * ((sqrt(n) * beta ^ 2 - 2 * a * lambda * 
		       abs(beta) + lambda ^ 2) / (2 * (a - 1)))
	} else
	{
		# If abs(beta) > a * lambda
		
		p <-  ((a + 1) * lambda ^ 2) / 2
	}
	stopifnot(p >= 0)
	return(p)
}
# ---------------------------------------------------------------------
# SCAD penalty function derivative
# ---------------------------------------------------------------------
p.scad.deriv <- function(beta, lambda, a, n) {
  # Calculates the SCAD penalty first derivative-   
  #
  # Args:
  #  beta   - the regression parameter
  #  lambda - > 0 parameter
  #  a      - > 2 parameter   
  #  n      - number of individuals
  #
  # Returns:
  #  The resultant value of the scad penalty
  
	if (sqrt(n) * abs(beta) <= lambda)
	{
		p.dash <- lambda * sqrt(n)
		
	} else if (lambda < sqrt(n) * abs(beta) )
	{
		# May want to check that max is the best thing to do here
		p.dash <- sqrt(n) * max((a * lambda - sqrt(n) * abs(beta)), 10e-3) / (a - 1)
	} 
	stopifnot(p.dash >= 0)
	return(p.dash)
}