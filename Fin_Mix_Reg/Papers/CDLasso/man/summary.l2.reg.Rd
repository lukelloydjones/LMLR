\name{summary.l2.reg}
\alias{summary.l2.reg}

\title{
Print all results of Cyclic Coordinate Descent for L2 Regression
}
\description{
Print full summary of results of Cyclic Coordinate Descent for L2 Regression.  
}
\usage{
\method{summary}{l2.reg}(object, ...)
}

\arguments{
  \item{object}{
Output of l2.reg. Must be of class "l2.reg"
}
  \item{...}{
N/A
}
}
\details{
\code{summary.l2.reg} produces full output from \link{l2.reg}. For selected output, see \link{print.l2.reg}.
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
  \item{L2}{
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
\code{\link{summary.l2.reg}}

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

#Lasso penalized L2 regression
out<-l2.reg(x,y,lambda=2)

#Re-estimate parameters without penalization
out2<-l2.reg(x[out$selected,],y,lambda=0)
summary(out2)
}