\name{cv.logit.reg}
\alias{cv.logit.reg}
\alias{cv.logit.reg.default}
\title{
k-fold Cross Validation
}
\description{
k-fold Cross Validation to find optimal \code{lambda} for Cyclic Coordinate Descent for logistic regression
}
\usage{
cv.logit.reg(x, y, k, lam.vec)
}
\arguments{
  \item{x}{
\code{p} x \code{n} design matrix - Note that the rows of X correspond to predictors and the columns to cases.
}
  \item{y}{
Outcome of length \code{n}. Outcome must be 0 and 1.
}
  \item{k}{
Number of folds for k-fold cross validation
}
  \item{lam.vec}{
Vector of penalization parameters
}
}
\details{
K-fold cross validation to select optimal lambda for use in cyclic coordinate descent for logistic regression \code{\link{logit.reg}}. The optimal value is considered the \code{lambda} value that retuns the lowest testing error over the cross validation. If more than one \code{lambda} value give the minumum testing error, the largest \code{lambda} is selected.
Plot of the cross validation can be viewed through \code{\link{plot.cv.logit.reg}}
}

\value{
  \item{k}{
  The value of \code{K} used for the K-fold cross validation.
  }
  \item{lam.vec}{
  The values of lambda tested.
  }
  \item{lam.opt}{
  The determined lambda value among \code{lam.vec} that returns the smallest prediction error. This value is the optimal \code{lambda} value for use in \code{\link{logit.reg}}.
  }
  \item{error.cv}{
  The prediction error matrix returned by cross validation method.
  }
  \item{num.pred}{
  The number of selected predictors when using the corresponding lambda value.
  }
}

\references{
	Wu, T.T., Chen, Y.F., Hastie, T., Sobel E. and Lange, K. (2009). Genome-wide association analysis by lasso penalized logistic regression. Bioinformatics, Volume 25, No 6, 714-721.

}
\author{
Edward Grant, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
}

\seealso{
\code{\link{logit.reg}}

\code{\link{plot.cv.logit.reg}}

}

\examples{
set.seed(1001)
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
cv

#Lasso penalized logistic regression using optimal lambda
out<-logit.reg(x,y,cv$lam.opt)

#Re-estimate parameters without penalization
out2<-logit.reg(x[out$selected,],y,0)
out2
}