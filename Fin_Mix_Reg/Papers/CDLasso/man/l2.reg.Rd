\name{l2.reg}
\alias{l2.reg}
\alias{l2.reg.default}

\title{
Cyclic Coordinate Descent for L2 regression
}
\description{
Cyclic Coordinate Descent for L2 regression with \emph{p} predictors and \emph{n} cases
}
\usage{
l2.reg(X, Y, lambda = 1)
}
\arguments{
  \item{X}{
\code{p} x \code{n} design matrix - Note that the rows of X correspond to predictors and the columns to cases.
}
  \item{Y}{
Outcome of length \code{n}
}
  \item{lambda}{
Penalization Parameter. For optimal \code{lambda}, use \link{cv.l2.reg}.
}
}

\details{
\code{l2.reg} performs an algorithm for estimating regression coefficients in a penalized L2 regression model. The algorithm is based on cyclic coordinate descent. For the new L1 algorithm that is faster, see (\link{l1.reg}).
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
  \item{L2}{
  The sum of the residuals
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
	Wu, T.T. and Lange, K. (2008). Coordinate Descent Algorithms for Lasso Penalized Regression. Annals of Applied Statistics, Volume 2, No 1, 224-244.

}

\author{
Edward Grant, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
}

\seealso{
\code{\link{print.l2.reg}}

\code{\link{summary.l2.reg}}

\code{\link{cv.l2.reg}}

\code{\link{plot.cv.l2.reg}}

\code{\link{l1.reg}}
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
out2
}