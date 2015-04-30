# Investigating SCAD penalty function

SCAD <- function(beta, lambda, a)
{
	if (abs(beta) <= lambda)
	{
		p <- lambda*abs(beta)
		
	} else if (lambda < abs(beta) & abs(beta) <= a * lambda)
	{
		p <- -((beta ^ 2 - 2 * a * lambda * abs(beta) + lambda ^ 2) / (2 * (a - 1)))
	} else
	{
		# If abs(beta) > a * lambda
		
		p <- ((a + 1) * lambda ^ 2) / 2
	}
	
	return(p)
}


beta <- seq(-5, 5, 0.01)
p <- unlist(lapply(beta, SCAD, lambda = 1, a = 3.7))
plot(beta, p)

SCADThreshold <- function(beta, lambda, a)
{
	if (abs(beta) <= 2 * lambda)
	{
		p <- max((abs(beta) - lambda), 0) * sign(beta)
		
	} else if (2 * lambda < abs(beta) & abs(beta) <= a * lambda)
	{
		p <- ((a - 1) * beta - sign(beta) * a * lambda) / (a - 2)
	} else
	{
		# If abs(beta) > a * lambda
		
		p <- beta
	}
	
	return(p)
}

beta <- seq(-5, 5, 0.01)
p <- unlist(lapply(beta, SCADThreshold, lambda = 1, a = 3.7))
plot(beta, p)