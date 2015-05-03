\name{plot.cv.l1.reg}
\alias{plot.cv.l1.reg}

\title{
Cross validation plot
}
\description{
Plot cross validation results across lambdas for greedy coordinate descent for L1 regression
}
\usage{
\method{plot}{cv.l1.reg}(x, ...)
}

\arguments{
  \item{x}{
Output of \link{cv.l1.reg}. Must be of class "cv.l1.reg"
}
  \item{...}{
N/A
}
}
\details{
\code{plot.cv.l1.reg} plots the prediction error of k fold cross validation across lambda values. 
}

\references{
	Wu, T.T. and Lange, K. (2008). Coordinate Descent Algorithms for Lasso Penalized Regression. Annals of Applied Statistics, Volume 2, No 1, 224-244.

}

\author{
Edward Grant, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
}

\seealso{
\code{\link{cv.l1.reg}}

\code{\link{l1.reg}}
}
\examples{
set.seed(100)
n=50 
p=200
nzfixed = c(1:5)
true.beta<-rep(0,p)
true.beta[nzfixed] = c(1,1,1,1,1)

x=matrix(rnorm(n*p),p,n)
y = t(x) \%*\% true.beta

rownames(x)<-1:nrow(x)
colnames(x)<-1:ncol(x)

crossval<-cv.l1.reg(x,y,10,(0:20)*2)
plot(crossval)
out<-l1.reg(x,y,lambda=crossval$lam.opt)
out
}