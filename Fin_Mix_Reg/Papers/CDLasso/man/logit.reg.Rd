\name{logit.reg}
\alias{logit.reg}
\alias{logit.reg.default}

\title{
Cyclic Coordinate Descent for Logistic regression
}
\description{
Cyclic Coordinate Descent for Logistic regression with \emph{p} predictors and \emph{n} cases
}
\usage{
logit.reg(X, Y, lambda = 1)
}
\arguments{
  \item{X}{
\code{p} x \code{n} design matrix - Note that the rows of X correspond to predictors and the columns to cases.
}
  \item{Y}{
Outcome of length \code{n}
}
  \item{lambda}{
Penalization Parameter. For optimal \code{lambda}, use \link{cv.logit.reg}.
}
}

\details{
\code{logit.reg} performs an algorithm for estimating regression coefficients in a penalized logistic regression model. The algorithm is based on cyclic coordinate descent. 
}

\value{
  \item{X}{
  The design matrix.
  }
  \item{cases}{
  The number of cases
  }
  \item{predictors}{
  The number of predictors
  }  
  \item{lambda}{
  The value of penalization parameter \code{lambda} used.
  }
  \item{residual}{
  A vector of length \code{p} listing the residuals
  }
  \item{estimate}{
  The estimate of the coefficients
  }
  \item{nonzeros}{
  The number "selected" variables included in the model.
  }
  \item{selected}{
  The name of the "selected" variables included in the model.
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
\code{\link{print.logit.reg}}

\code{\link{summary.logit.reg}}

\code{\link{cv.logit.reg}}

\code{\link{plot.cv.logit.reg}}

\code{\link{l1.reg}}
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
out2

}