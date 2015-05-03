\name{summary.l1.reg}
\alias{summary.l1.reg}

\title{
Print all results of Greedy Coordinate Descent for L1 Regression
}
\description{
Print full summary of results of Greedy Coordinate Descent for L1 Regression.  
}
\usage{
\method{summary}{l1.reg}(object, ...)
}

\arguments{
  \item{object}{
Output of l1.reg. Must be of class "l1.reg"
}
  \item{...}{
N/A
}
}
\details{
\code{summary.l1.reg} produces full output from \link{l1.reg}. For selected output, see \link{print.l1.reg}.
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
summary(out2)
}