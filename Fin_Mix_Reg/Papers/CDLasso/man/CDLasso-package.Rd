\name{CDLasso-package}
\alias{CDLasso-package}
\alias{CDLasso}
\docType{package}
\title{
Coordinate descent algorithms for L1 and L2 regression
}
\description{
Greedy coordinate descent for L1 regression and cyclic coordinate descent for L2 regression with \emph{p} predictors and \emph{n} cases
}
\details{
\tabular{ll}{
Package: \tab CDLasso \cr
Title: \tab Coordinate Descent for L1 and L2 Regression \cr
Version: \tab 1.1 \cr
Date: \tab 2013-13-03 \cr
Author: \tab Edward Grant, Kenneth Lange, Tong Tong Wu \cr
Maintainer: \tab Edward Grant <edward.m.grant@gmail.com> \cr
Description: \tab Coordinate Descent for L1, L2, and Logistic Regression \cr
License: \tab GPL-2 \cr
}

}
\author{
Edward Grant, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
}

\references{
	Wu, T.T. and Lange, K. (2008). Coordinate Descent Algorithms for Lasso Penalized Regression. Annals of Applied Statistics, Volume 2, No 1, 224-244.

}

\seealso{
\code{\link{l1.reg}}

\code{\link{l2.reg}}

\code{\link{logit.reg}}
}
\examples{
set.seed(1001)
n=500
p=2000
nz = c(1:5)
true.beta<-rep(0,p)
true.beta[nz] = c(1,1,1,1,1)

x=matrix(rnorm(n*p),p,n)
y = t(x) \%*\% true.beta

logity=exp(y)/(1+exp(y))
ylog=rbinom(n=length(logity),prob=logity,size=1)

rownames(x)<-1:nrow(x)
colnames(x)<-1:ncol(x)

#L1
outL1<-l1.reg(x,y,lambda=50)
outL1est<-l1.reg(x[outL1$selected,],y,lambda=0)

#L2
outL2<-l2.reg(x,y,2)
outL2est<-l2.reg(x[outL2$selected,],y,lambda=0)

#Logistic
outLOGIT<-logit.reg(x,ylog,lambda=50)
outLOGITest<-logit.reg(x[outLOGIT$selected,],ylog,lambda=0)

}
