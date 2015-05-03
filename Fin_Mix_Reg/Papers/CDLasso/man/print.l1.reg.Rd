\name{print.l1.reg}
\alias{print.l1.reg}

\title{
Print results of Greedy Coordinate Descent for L1 Regression
}
\description{
Print short summary of results of Greedy Coordinate Descent for L1 Regression. Includes number of cases and predictors, lambda used, estimate of coeffcients produced, the number of selected predictors, and the names of selected predictors. 
}
\usage{
\method{print}{l1.reg}(x, ...)
}

\arguments{
  \item{x}{
Output of l1.reg. Must be of class "l1.reg"
}
  \item{...}{
N/A
}
}
\details{
\code{print.l1.reg} produces selected output from \code{l1.reg}. For more output, see \link{summary.l1.reg}.
}

\references{
	Wu, T.T. and Lange, K. (2008). Coordinate Descent Algorithms for Lasso Penalized Regression. Annals of Applied Statistics, Volume 2, No 1, 224-244.

}

\author{
Edward Grant, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
}

\seealso{
\code{\link{summary.l1.reg}}

\code{\link{l1.reg}}
}
\examples{
set.seed(100)
n=500
p=2000
nz = c(1:5)
true.beta<-rep(0,p)
true.beta[nz] = c(1,1,1,1,1)

x=matrix(rnorm(n*p),p,n)
y = t(x) \%*\% true.beta

rownames(x)<-1:nrow(x)
colnames(x)<-1:ncol(x)

#Lasso penalized L1 regression
out<-l1.reg(x,y,lambda=50)

#Re-estimate parameters without penalization
out2<-l1.reg(x[out$selected,],y,lambda=0)
print(out2)
}