\name{print.logit.reg}
\alias{print.logit.reg}

\title{
Print results of Cyclic Coordinate Descent for Logistic Regression
}
\description{
Print short summary of results of Cyclic Coordinate Descent for Logistic Regression. Includes number of cases and predictors, lambda used, estimate of coeffcients produced, the number of selected predictors, and the names of selected predictors. 
}
\usage{
\method{print}{logit.reg}(x, ...)
}

\arguments{
  \item{x}{
Output of logit.reg. Must be of class "logit.reg"
}
  \item{...}{
N/A
}
}
\details{
\code{print.logit.reg} produces selected output from logit.reg. For more output, see \link{summary.logit.reg}.
}

\references{
	Wu, T.T., Chen, Y.F., Hastie, T., Sobel E. and Lange, K. (2009). Genome-wide association analysis by lasso penalized logistic regression. Bioinformatics, Volume 25, No 6, 714-721.

}

\author{
Edward Grant, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
}

\seealso{
\code{\link{summary.logit.reg}}

\code{\link{logit.reg}}
}
\examples{
set.seed(1001)
n=500;p=5000
beta=c(1,1,1,1,1,rep(0,p-5))
x=matrix(rnorm(n*p),p,n)
xb = t(x) \%*\% beta
logity=exp(xb)/(1+exp(xb))
y=rbinom(n=length(logity),prob=logity,size=1)

rownames(x)<-1:nrow(x)
colnames(x)<-1:ncol(x)

#Lasso penalized logistic regression using optimal lambda
out<-logit.reg(x,y,50)
print(out)

#Re-estimate parameters without penalization
out2<-logit.reg(x[out$selected,],y,0)
print(out2)
}