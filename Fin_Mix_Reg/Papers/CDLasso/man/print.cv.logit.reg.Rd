\name{print.cv.logit.reg}
\alias{print.cv.logit.reg}

\title{
Print results of Greedy Coordinate Descent for Logistic Regression
}
\description{
Print short summary of results of cross validation for Greedy Coordinate Descent for Logistric Regression. 
}
\usage{
\method{print}{cv.logit.reg}(x, ...)
}

\arguments{
  \item{x}{
Output of cv.logit.reg. Must be of class "cv.logit.reg"
}
  \item{...}{
N/A
}
}
\details{
\code{print.cv.logit.reg} produces output from \code{cv.logit.reg}.
}

\references{
	Wu, T.T., Chen, Y.F., Hastie, T., Sobel E. and Lange, K. (2009). Genome-wide association analysis by lasso penalized logistic regression. Bioinformatics, Volume 25, No 6, 714-721.

}

\author{
Edward Grant, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
}

\seealso{
\code{\link{cv.logit.reg}}

}
\examples{
set.seed(101)
n=250;p=50
beta=c(1,1,1,1,1,rep(0,p-5))
x=matrix(rnorm(n*p),p,n)
xb = t(x) \%*\% beta
logity=exp(xb)/(1+exp(xb))
y=rbinom(n=length(logity),prob=logity,size=1)

rownames(x)<-1:nrow(x)
colnames(x)<-1:ncol(x)
lam.vec = (0:15)*2

#K-fold cross validation
cv <- cv.logit.reg(x,y,5,lam.vec)
plot(cv)

#Lasso penalized logistic regression using optimal lambda
out<-logit.reg(x,y,cv$lam.opt)

#Re-estimate parameters without penalization
out2<-logit.reg(x[out$selected,],y,0)
out2$estimate
}