\name{l1.reg}
\alias{l1.reg}
\alias{l1.reg.default}

\title{
Greedy Coordinate Descent for L1 regression
}
\description{
Greedy Coordinate Descent for L1 regression with \emph{p} predictors and \emph{n} cases
}
\usage{
	l1.reg(X, Y, lambda = 1)
}
\arguments{
  \item{X}{
\code{p} x \code{n} design matrix - Note that the rows of X correspond to predictors and the columns to cases.
}
  \item{Y}{
Outcome of length \code{n}
}
  \item{lambda}{
Penalization Parameter. To find optimal \code{lambda}, use \link{cv.l1.reg}.
}
}

\details{
\code{l1.reg} performs a new algorithm for estimating regression coefficients with a lasso penalty. The algorithm is based on greedy coordinate descent and Edgeworth's algorithm for ordinary L1 regression. This L1 algorithm is faster than the cyclic coordinate descent in L2 regression (\link{l2.reg}).
}

\value{
  \item{X}{
  The design matrix.
  }
  \item{Y}{
  The outcome variable for cases.
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
  \item{objective}{
  The value of the objective function
  }
  \item{residual}{
  A vector of length \code{p} listing the residuals
  }
  \item{L1}{
  The sum of the residuals
  }
  \item{estimate}{
  The estimate of the coefficients
  }
  \item{nonzeros}{
  The name of "selected" variables included in the model.
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
\code{\link{print.l1.reg}}

\code{\link{summary.l1.reg}}

\code{\link{cv.l1.reg}}

\code{\link{plot.cv.l1.reg}}

\code{\link{l2.reg}}
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
out2
}