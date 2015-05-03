\name{print.l2.reg}
\alias{print.l2.reg}

\title{
Print results of Cyclic Coordinate Descent for L2 Regression
}
\description{
Print short summary of results of Cyclic Coordinate Descent for L2 Regression. Includes number of cases and predictors, lambda used, estimate of coeffcients produced, the number of selected predictors, and the names of selected predictors. 
}
\usage{
\method{print}{l2.reg}(x, ...)
}

\arguments{
  \item{x}{
Output of l2.reg. Must be of class "l2.reg"
}
  \item{...}{
N/A
}
}
\details{
\code{print.L1_REG} produces selected output from l2.reg. For more output, see \link{summary.l2.reg}.
}

\references{
	Wu, T.T. and Lange, K. (2008). Coordinate Descent Algorithms for Lasso Penalized Regression. Annals of Applied Statistics, Volume 2, No 1, 224-244.

}

\author{
Edward Grant, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
}

\seealso{
\code{\link{summary.l2.reg}}

\code{\link{l2.reg}}
}
\examples{
set.seed(100)
n=500
p=2000
nzfixed = c(1:5)
true.beta<-rep(0,p)
true.beta[nzfixed] = c(1,1,1,1,1)

x=matrix(rnorm(n*p),p,n)
y = t(x) \%*\% true.beta

rownames(x)<-1:nrow(x)
colnames(x)<-1:ncol(x)

#Lasso penalized L2 regression
out<-l2.reg(x,y,lambda=2)

#Re-estimate parameters without penalization
out2<-l2.reg(x[out$selected,],y,lambda=0)
print(out2)
}